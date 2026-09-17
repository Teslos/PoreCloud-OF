/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see MeltPoolInjection.H
\*---------------------------------------------------------------------------*/

#include "MeltPoolInjection.H"
#include "mathematicalConstants.H"
#include "polyMeshTetDecomposition.H"
#include "globalIndex.H"
#include "Pstream.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::MeltPoolInjection<CloudType>::generatePositions(const scalar dt)
{
    positions_.clear();
    injectorCells_.clear();
    injectorTetFaces_.clear();
    injectorTetPts_.clear();
    diameters_.clear();
    velocities_.clear();

    const fvMesh& mesh = this->owner().mesh();
    const scalarField& V = mesh.V();
    Random& rnd = this->owner().rndGen();

    const volScalarField* alphaPtr =
        mesh.template findObject<volScalarField>(alphaName_);
    const volScalarField* epsPtr =
        mesh.template findObject<volScalarField>(epsilonName_);

    if (!alphaPtr)
    {
        WarningInFunction
            << "Field " << alphaName_ << " not found; no bubbles injected."
            << nl;
        return;
    }

    const scalarField& alpha = alphaPtr->primitiveField();

    // Local velocity field, used when U0 was not specified
    const volVectorField* UPtr = mesh.template findObject<volVectorField>("U");

    DynamicList<vector> positions(256);
    DynamicList<label> injectorCells(256);
    DynamicList<label> injectorTetFaces(256);
    DynamicList<label> injectorTetPts(256);
    DynamicList<vector> velocities(256);

    // Fractional carry-over, so a rate that yields far less than one bubble
    // per cell per step still produces the correct long-run count instead of
    // truncating to zero everywhere.
    scalar newParticlesTotal = fractionalCarry_;
    label addParticlesTotal = 0;

    forAll(alpha, celli)
    {
        // Liquid metal: the solver's own definition (alpha.metal and epsilon1
        // both above 0.5; laserbeamFoam writes this as its "condition" field).
        if (alpha[celli] <= alphaMin_)
        {
            continue;
        }

        if (epsPtr && epsPtr->primitiveField()[celli] <= epsMin_)
        {
            continue;
        }

        newParticlesTotal += rate_*V[celli]*dt;

        label addParticles = 0;
        const scalar diff = newParticlesTotal - addParticlesTotal;
        if (diff >= 1)
        {
            addParticles = floor(diff);
            addParticlesTotal += addParticles;
        }

        if (addParticles == 0)
        {
            continue;
        }

        // Decompose the cell into tets and pick one weighted by volume, so
        // positions are uniform in space rather than biased to the centroid.
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, celli);

        scalarList cTetVFrac(cellTetIs.size(), Zero);
        for (label tetI = 1; tetI < cellTetIs.size() - 1; tetI++)
        {
            cTetVFrac[tetI] =
                cTetVFrac[tetI-1] + cellTetIs[tetI].tet(mesh).mag()/V[celli];
        }
        cTetVFrac.last() = 1.0;

        for (label pI = 0; pI < addParticles; pI++)
        {
            const scalar volFrac = rnd.sample01<scalar>();
            label tetI = 0;
            forAll(cTetVFrac, vfI)
            {
                if (cTetVFrac[vfI] > volFrac)
                {
                    tetI = vfI;
                    break;
                }
            }

            positions.append(cellTetIs[tetI].tet(mesh).randomPoint(rnd));
            injectorCells.append(celli);
            injectorTetFaces.append(cellTetIs[tetI].face());
            injectorTetPts.append(cellTetIs[tetI].tetPt());

            // Born moving with the melt unless told otherwise: a bubble
            // nucleating in flowing metal is already being advected.
            velocities.append
            (
                haveU0_ || !UPtr ? U0_ : UPtr->primitiveField()[celli]
            );
        }
    }

    fractionalCarry_ = newParticlesTotal - addParticlesTotal;

    // Parallel: gather every processor's positions onto every processor, so
    // the parcelI indexing in inject() is globally consistent.  Cell/tet ids
    // stay -1 on processors that do not own the position, which is what
    // validInjection() keys on.
    const label myProci = UPstream::myProcNo();
    globalIndex globalPositions(positions.size());

    const label total = globalPositions.totalSize();

    positions_.setSize(total, point::max);
    injectorCells_.setSize(total, -1);
    injectorTetFaces_.setSize(total, -1);
    injectorTetPts_.setSize(total, -1);
    velocities_.setSize(total, Zero);

    SubList<vector>(positions_, globalPositions.range(myProci)) = positions;
    SubList<label>(injectorCells_, globalPositions.range(myProci)) =
        injectorCells;
    SubList<label>(injectorTetFaces_, globalPositions.range(myProci)) =
        injectorTetFaces;
    SubList<label>(injectorTetPts_, globalPositions.range(myProci)) =
        injectorTetPts;
    SubList<vector>(velocities_, globalPositions.range(myProci)) = velocities;

    if (UPstream::parRun())
    {
        Pstream::listReduce(positions_, minOp<point>());
        Pstream::listReduce(injectorCells_, maxOp<label>());
        Pstream::listReduce(injectorTetFaces_, maxOp<label>());
        Pstream::listReduce(injectorTetPts_, maxOp<label>());
        Pstream::listReduce(velocities_, maxOp<vector>());
    }

    // Diameters must be identical on every processor, and are drawn after the
    // gather so that the RNG sequence does not have to match across ranks.
    diameters_.setSize(total, 0.0);
    forAll(diameters_, i)
    {
        diameters_[i] = sizeDistribution_->sample();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MeltPoolInjection<CloudType>::MeltPoolInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    rate_(this->coeffDict().getScalar("rate")),
    alphaName_
    (
        this->coeffDict().template getOrDefault<word>("alphaName", "alpha.metal")
    ),
    epsilonName_
    (
        this->coeffDict().template getOrDefault<word>("epsilonName", "epsilon1")
    ),
    alphaMin_(this->coeffDict().template getOrDefault<scalar>("alphaMin", 0.5)),
    epsMin_(this->coeffDict().template getOrDefault<scalar>("epsMin", 0.5)),
    duration_(this->coeffDict().getScalar("duration")),
    maxParcels_
    (
        this->coeffDict().template getOrDefault<label>("maxParcels", 50000)
    ),
    haveU0_(this->coeffDict().found("U0")),
    U0_
    (
        this->coeffDict().template getOrDefault<vector>("U0", Zero)
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    nInjectedTotal_(0),
    fractionalCarry_(0),
    positions_(),
    injectorCells_(),
    injectorTetFaces_(),
    injectorTetPts_(),
    diameters_(),
    velocities_()
{
    // Only used as the normaliser in the base class's injection trigger
    // (newVolumeFraction = volumeToInject/volumeTotal_ must exceed zero).
    // parcelBasisType is expected to be 'fixed', so it plays no part in
    // setting the number of particles per parcel.
    this->volumeTotal_ = 1.0;

    Info<< "    Melt-pool bubble injection:" << nl
        << "        rate            : " << rate_ << " 1/m^3/s" << nl
        << "        liquid test     : " << alphaName_ << " > " << alphaMin_
        << " and " << epsilonName_ << " > " << epsMin_ << nl
        << "        duration        : " << duration_ << " s" << nl
        << "        maxParcels      : " << maxParcels_ << nl;
}


template<class CloudType>
Foam::MeltPoolInjection<CloudType>::MeltPoolInjection
(
    const MeltPoolInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    rate_(im.rate_),
    alphaName_(im.alphaName_),
    epsilonName_(im.epsilonName_),
    alphaMin_(im.alphaMin_),
    epsMin_(im.epsMin_),
    duration_(im.duration_),
    maxParcels_(im.maxParcels_),
    haveU0_(im.haveU0_),
    U0_(im.U0_),
    sizeDistribution_(im.sizeDistribution_.clone()),
    nInjectedTotal_(im.nInjectedTotal_),
    fractionalCarry_(im.fractionalCarry_),
    positions_(im.positions_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    diameters_(im.diameters_),
    velocities_(im.velocities_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MeltPoolInjection<CloudType>::updateMesh()
{
    // Lists are rebuilt from scratch every timestep in parcelsToInject(), so
    // a topology change needs no remapping here - only invalidation.
    positions_.clear();
    injectorCells_.clear();
    injectorTetFaces_.clear();
    injectorTetPts_.clear();
    diameters_.clear();
    velocities_.clear();
}


template<class CloudType>
Foam::scalar Foam::MeltPoolInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::MeltPoolInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    // time0/time1 are relative to SOI
    if (time1 <= 0 || time0 >= duration_)
    {
        return 0;
    }

    if (nInjectedTotal_ >= maxParcels_)
    {
        return 0;
    }

    const scalar dt = max(scalar(0), min(time1, duration_) - max(time0, scalar(0)));

    if (dt <= 0)
    {
        return 0;
    }

    generatePositions(dt);

    label n = positions_.size();

    // Respect the safety cap
    if (nInjectedTotal_ + n > maxParcels_)
    {
        n = max(label(0), maxParcels_ - nInjectedTotal_);
        positions_.setSize(n);
        injectorCells_.setSize(n);
        injectorTetFaces_.setSize(n);
        injectorTetPts_.setSize(n);
        diameters_.setSize(n);
        velocities_.setSize(n);
    }

    nInjectedTotal_ += n;

    return n;
}


template<class CloudType>
Foam::scalar Foam::MeltPoolInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar vol = 0;

    forAll(diameters_, i)
    {
        const scalar d = diameters_[i];
        vol += constant::mathematical::pi/6.0*d*d*d;
    }

    return vol;
}


template<class CloudType>
void Foam::MeltPoolInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
    tetFacei = injectorTetFaces_[parcelI];
    tetPti = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::MeltPoolInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    parcel.U() = velocities_[parcelI];
    parcel.d() = diameters_[parcelI];
}


template<class CloudType>
bool Foam::MeltPoolInjection<CloudType>::fullyDescribed() const
{
    // Diameter and velocity are set above; density comes from
    // constantProperties via the base class.
    return false;
}


template<class CloudType>
bool Foam::MeltPoolInjection<CloudType>::validInjection(const label parcelI)
{
    return injectorCells_[parcelI] > -1;
}


// ************************************************************************* //
