/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see LeenovKolinForce.H
\*---------------------------------------------------------------------------*/

#include "LeenovKolinForce.H"
#include "poreCloudSphereAverage.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LeenovKolinForce<CloudType>::LeenovKolinForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    JName_(this->coeffs().template getOrDefault<word>("J", "J_MHD")),
    exclusionCoeff_
    (
        this->coeffs().template getOrDefault<scalar>("exclusionCoeff", -1.5)
    ),
    useSphereAverage_
    (
        this->coeffs().template getOrDefault<bool>("useSphereAverage", true)
    ),
    BFieldPtr_(nullptr),
    JPtr_(nullptr),
    JInterpPtr_(nullptr),
    BInterpPtr_(nullptr)
{
    Info<< "    Leenov-Kolin exclusion force:" << nl
        << "        J field         : " << JName_ << nl
        << "        exclusionCoeff  : " << exclusionCoeff_ << nl
        << "        sphere-averaged : " << Switch(useSphereAverage_) << nl;
}


template<class CloudType>
Foam::LeenovKolinForce<CloudType>::LeenovKolinForce
(
    const LeenovKolinForce& lkf
)
:
    ParticleForce<CloudType>(lkf),
    JName_(lkf.JName_),
    exclusionCoeff_(lkf.exclusionCoeff_),
    useSphereAverage_(lkf.useSphereAverage_),
    BFieldPtr_(nullptr),      // rebuilt on the next cacheFields
    JPtr_(nullptr),
    JInterpPtr_(nullptr),
    BInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LeenovKolinForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const fvMesh& mesh = this->mesh();

        if (!BFieldPtr_)
        {
            BFieldPtr_.reset(new poreCloud::magneticField(mesh));
            BFieldPtr_->info();
        }

        // laserbeamFoam's functionObjects fire at the top of the time loop,
        // before ++runTime, so runTime.value() is the time of the completed
        // step - which is the time the cached J belongs to.
        BFieldPtr_->update(mesh.time().value());

        JPtr_ = mesh.findObject<volVectorField>(JName_);

        if (!JPtr_)
        {
            // MHD off, or the solver did not register J: force stays zero
            // rather than aborting the run.
            WarningInFunction
                << "Current density field " << JName_
                << " not found in the registry; Leenov-Kolin force is zero."
                << nl;
            return;
        }

        if (!useSphereAverage_)
        {
            JInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    *JPtr_
                ).ptr()
            );

            BInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    BFieldPtr_->B()
                ).ptr()
            );
        }
    }
    else
    {
        JInterpPtr_.clear();
        BInterpPtr_.clear();
        JPtr_ = nullptr;
    }
}


template<class CloudType>
Foam::forceSuSp Foam::LeenovKolinForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    if (!JPtr_ || !BFieldPtr_ || !BFieldPtr_->active())
    {
        return value;
    }

    const scalar d = p.d();
    const scalar V = constant::mathematical::pi/6.0*d*d*d;

    vector J0(Zero);
    vector B(Zero);

    if (useSphereAverage_)
    {
        const fvMesh& mesh = this->mesh();
        const point pos = p.position();
        const scalar a = 0.5*d;

        label nJ = 0;
        label nB = 0;

        J0 = poreCloud::sphereAverage
        (
            mesh, JPtr_->primitiveField(), p.cell(), pos, a, nJ
        );

        B = poreCloud::sphereAverage
        (
            mesh, BFieldPtr_->B().primitiveField(), p.cell(), pos, a, nB
        );
    }
    else
    {
        J0 = JInterpPtr_().interpolate(p.coordinates(), p.currentTetIndices());
        B  = BInterpPtr_().interpolate(p.coordinates(), p.currentTetIndices());
    }

    // Body force: explicit only.  KinematicParcel::calcVelocity discards
    // Fncp.Sp() (the implicit non-coupled splitting is commented out there),
    // so anything returned in Sp() would be silently dropped.
    value.Su() = exclusionCoeff_*V*(J0 ^ B);

    return value;
}


// ************************************************************************* //
