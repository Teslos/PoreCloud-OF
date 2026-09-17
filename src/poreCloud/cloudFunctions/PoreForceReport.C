/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see PoreForceReport.H
\*---------------------------------------------------------------------------*/

#include "PoreForceReport.H"
#include "poreCloudSphereAverage.H"
#include "mathematicalConstants.H"
#include "gravityMeshObject.H"
#include "volFields.H"
#include "Pstream.H"
#include <iomanip>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PoreForceReport<CloudType>::PoreForceReport
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    writeControl_
    (
        this->coeffDict().template getOrDefault<word>("writeControl", "writeTime")
    ),
    writeInterval_
    (
        this->coeffDict().template getOrDefault<label>("writeInterval", 1)
    ),
    exclusionCoeff_
    (
        this->coeffDict().template getOrDefault<scalar>("exclusionCoeff", -1.5)
    ),
    JName_(this->coeffDict().template getOrDefault<word>("J", "J_MHD")),
    useSphereAverage_
    (
        this->coeffDict().template getOrDefault<bool>("useSphereAverage", true)
    ),
    epsilonName_
    (
        this->coeffDict().template getOrDefault<word>("epsilonName", "epsilon1")
    ),
    csvName_
    (
        this->coeffDict().template getOrDefault<fileName>
        (
            "file",
            owner.name() + "_pores.csv"
        )
    ),
    BFieldPtr_(nullptr),
    csvPtr_(nullptr),
    stepCounter_(0),
    birthPos_()
{
    if (writeControl_ != "writeTime" && writeControl_ != "timeStep")
    {
        FatalErrorInFunction
            << "Unknown writeControl " << writeControl_
            << "; expected writeTime or timeStep"
            << exit(FatalError);
    }

    // Cross-check against the force model actually integrated, so a mismatched
    // coefficient cannot silently make the reported force disagree with the
    // one that moved the bubble.
    const dictionary& forcesDict =
        owner.subModelProperties().subDict("particleForces");

    if (forcesDict.found("leenovKolin"))
    {
        const dictionary& lk = forcesDict.subDict("leenovKolin");
        const scalar cForce = lk.getOrDefault<scalar>("exclusionCoeff", -1.5);

        if (mag(cForce - exclusionCoeff_) > SMALL)
        {
            WarningInFunction
                << "exclusionCoeff in poreForceReport (" << exclusionCoeff_
                << ") differs from the leenovKolin force model ("
                << cForce << "). The reported Fx_N will not be the force "
                << "that was actually integrated." << nl;
        }

        const bool sForce = lk.getOrDefault<bool>("useSphereAverage", true);

        if (sForce != useSphereAverage_)
        {
            WarningInFunction
                << "useSphereAverage in poreForceReport ("
                << Switch(useSphereAverage_)
                << ") differs from the leenovKolin force model ("
                << Switch(sForce) << "). The reported Fx_N will not be the "
                << "force that was actually integrated." << nl;
        }
    }

    // The stream is opened lazily on first write, not here: cloud function
    // objects are cloned into the cloud's list, and a clone constructed from
    // the copy constructor would otherwise hold a null stream while the
    // original owned the file.
}


template<class CloudType>
Foam::PoreForceReport<CloudType>::PoreForceReport
(
    const PoreForceReport<CloudType>& pfr
)
:
    CloudFunctionObject<CloudType>(pfr),
    writeControl_(pfr.writeControl_),
    writeInterval_(pfr.writeInterval_),
    exclusionCoeff_(pfr.exclusionCoeff_),
    JName_(pfr.JName_),
    useSphereAverage_(pfr.useSphereAverage_),
    epsilonName_(pfr.epsilonName_),
    csvName_(pfr.csvName_),
    BFieldPtr_(nullptr),
    csvPtr_(nullptr),
    stepCounter_(pfr.stepCounter_),
    birthPos_(pfr.birthPos_)
{}


// * * * * * * * * * * * * * * * * Private  * * * * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PoreForceReport<CloudType>::openCSV()
{
    if (!Pstream::master())
    {
        return;
    }

    // globalPath(), not path(): under MPI, Time::path() resolves to the
    // processor-local case directory, which is why PoreTracker's CSV lands in
    // processor0/ instead of the case root.
    const fileName outFile =
        this->owner().mesh().time().globalPath()/csvName_;

    const bool exists = Foam::isFile(outFile);

    // Append rather than truncate, so a restart continues the series instead
    // of erasing everything written before it.
    csvPtr_.reset
    (
        new std::ofstream
        (
            outFile,
            exists ? std::ios::app : std::ios::out
        )
    );

    if (!exists)
    {
        *csvPtr_
            << "Time,PoreID,IsKeyhole,Volume_m3,Cx_m,Cy_m,Cz_m,"
            << "Fx_N,Fy_N,Fz_N,"
            << "BirthTime,BirthCx,BirthCy,BirthCz,NCells,"
            << "Ux,Uy,Uz,d_m,dOverDx,Active,"
            << "Uslipx,Uslipy,Uslipz,"
            << "Fbuoyx,Fbuoyy,Fbuoyz,"
            << "T_K,epsilon1,Coverage"
            << std::endl;
    }

    Info<< "    poreForceReport: writing " << outFile << nl;
}


template<class CloudType>
bool Foam::PoreForceReport<CloudType>::shouldWrite() const
{
    const Time& runTime = this->owner().mesh().time();

    if (writeControl_ == "writeTime")
    {
        return runTime.writeTime();
    }

    return (stepCounter_ % max(label(1), writeInterval_)) == 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PoreForceReport<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    ++stepCounter_;

    if (!shouldWrite())
    {
        return;
    }

    if (!csvPtr_ && Pstream::master())
    {
        openCSV();
    }

    const fvMesh& mesh = this->owner().mesh();
    const Time& runTime = mesh.time();
    const scalar t = runTime.value();

    if (!BFieldPtr_)
    {
        BFieldPtr_.reset(new poreCloud::magneticField(mesh));
    }
    BFieldPtr_->update(t);

    const volVectorField* JPtr = mesh.template findObject<volVectorField>(JName_);
    const volScalarField* epsPtr =
        mesh.template findObject<volScalarField>(epsilonName_);
    const volScalarField* TPtr = mesh.template findObject<volScalarField>("T");
    const volScalarField* rhoPtr =
        mesh.template findObject<volScalarField>("rho");
    const volVectorField* UPtr = mesh.template findObject<volVectorField>("U");

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(runTime);

    // Pack every local parcel into a flat scalar list for the gather.  A POD
    // struct would need IO operators; a flat list is what PoreTracker uses for
    // the same reason.
    DynamicList<scalar> local(nFields*this->owner().size());

    for (const parcelType& p : this->owner())
    {
        const label celli = p.cell();
        const scalar d = p.d();
        const scalar a = 0.5*d;
        const scalar V = constant::mathematical::pi/6.0*d*d*d;
        const point pos = p.position();

        // Electromagnetic exclusion force over the bubble sphere - the
        // Lagrangian analogue of PoreTracker's volume sum.
        vector Fem(Zero);
        label nCells = 1;
        scalar cover = 1.0;

        if (JPtr && BFieldPtr_->active() && celli >= 0)
        {
            label nJ = 1;
            label nB = 1;
            scalar sumVJ = 0;
            scalar sumVB = 0;

            // Radius zero collapses sphereAverage to the host-cell value, so
            // the two sampling modes share one code path.
            const scalar aSample = useSphereAverage_ ? a : 0.0;

            const vector J0 = poreCloud::sphereAverage
            (
                mesh, JPtr->primitiveField(), celli, pos, aSample, nJ, sumVJ
            );
            const vector B = poreCloud::sphereAverage
            (
                mesh, BFieldPtr_->B().primitiveField(), celli, pos, aSample,
                nB, sumVB
            );

            Fem = exclusionCoeff_*V*(J0 ^ B);
            nCells = nJ;

            if (useSphereAverage_)
            {
                cover = poreCloud::coverage(sumVJ, a);
            }
        }

        const scalar rhoc =
            (rhoPtr && celli >= 0) ? rhoPtr->primitiveField()[celli] : 0.0;

        const vector Uc =
            (UPtr && celli >= 0) ? UPtr->primitiveField()[celli] : vector::zero;

        // Net buoyancy including the bubble's own weight
        const vector Fbuoy = V*(rhoc - p.rho())*g.value();

        const scalar Tcell =
            (TPtr && celli >= 0) ? TPtr->primitiveField()[celli] : 0.0;

        const scalar eps =
            (epsPtr && celli >= 0) ? epsPtr->primitiveField()[celli] : 0.0;

        const scalar dx = poreCloud::cellLength(mesh, celli);

        const label poreID = p.origProc()*1000000 + p.origId();

        local.append(scalar(poreID));
        local.append(V);
        local.append(pos.x());
        local.append(pos.y());
        local.append(pos.z());
        local.append(Fem.x());
        local.append(Fem.y());
        local.append(Fem.z());
        local.append(t - p.age());          // exact birth time
        local.append(p.U().x());
        local.append(p.U().y());
        local.append(p.U().z());
        local.append(d);
        local.append(d/max(dx, SMALL));
        local.append(p.active() ? 1.0 : 0.0);
        local.append(p.U().x() - Uc.x());
        local.append(p.U().y() - Uc.y());
        local.append(p.U().z() - Uc.z());
        local.append(Fbuoy.x());
        local.append(Fbuoy.y());
        local.append(Fbuoy.z());
        local.append(Tcell);
        local.append(eps);
        local.append(scalar(nCells));
        local.append(cover);
        local.append(0.0);                  // reserved
    }

    // Gather to master
    List<List<scalar>> allData(Pstream::nProcs());
    allData[Pstream::myProcNo()] = local;
    Pstream::gatherList(allData);

    if (!Pstream::master() || !csvPtr_)
    {
        return;
    }

    std::ofstream& os = *csvPtr_;
    os << std::setprecision(10);

    for (const List<scalar>& procData : allData)
    {
        const label n = procData.size()/nFields;

        for (label i = 0; i < n; ++i)
        {
            const scalar* r = &procData[i*nFields];

            const label poreID = label(r[0]);
            const point c(r[2], r[3], r[4]);

            // Record the birth position the first time this bubble is seen.
            auto it = birthPos_.find(poreID);
            if (it == birthPos_.end())
            {
                it = birthPos_.emplace(poreID, c).first;
            }
            const point& b = it->second;

            os  << t << ','
                << poreID << ','
                << 0 << ','                      // IsKeyhole
                << r[1] << ','                   // Volume_m3
                << r[2] << ',' << r[3] << ',' << r[4] << ','
                << r[5] << ',' << r[6] << ',' << r[7] << ','
                << r[8] << ','                   // BirthTime
                << b.x() << ',' << b.y() << ',' << b.z() << ','
                << label(r[23]) << ','           // NCells
                << r[9] << ',' << r[10] << ',' << r[11] << ','
                << r[12] << ','                  // d_m
                << r[13] << ','                  // dOverDx
                << label(r[14]) << ','           // Active
                << r[15] << ',' << r[16] << ',' << r[17] << ','
                << r[18] << ',' << r[19] << ',' << r[20] << ','
                << r[21] << ',' << r[22] << ','
                << r[24]                         // Coverage
                << '\n';
        }
    }

    os.flush();
}


// ************************************************************************* //
