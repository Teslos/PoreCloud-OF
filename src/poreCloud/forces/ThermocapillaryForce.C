/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see ThermocapillaryForce.H
\*---------------------------------------------------------------------------*/

#include "ThermocapillaryForce.H"
#include "poreCloudSphereAverage.H"
#include "mathematicalConstants.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermocapillaryForce<CloudType>::ThermocapillaryForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    gradTName_(this->coeffs().template getOrDefault<word>("gradT", "gradT")),
    epsilonName_
    (
        this->coeffs().template getOrDefault<word>("epsilon1", "epsilon1")
    ),
    dSigmadT_(0),
    viscosityRatio_
    (
        this->coeffs().template getOrDefault<scalar>("viscosityRatio", 0)
    ),
    conductivityRatio_
    (
        this->coeffs().template getOrDefault<scalar>("conductivityRatio", 0)
    ),
    epsMin_(this->coeffs().template getOrDefault<scalar>("epsMin", 0.5)),
    Cmig_(this->coeffs().template getOrDefault<scalar>("Cmig", 1.0)),
    useSphereAverage_
    (
        this->coeffs().template getOrDefault<bool>("useSphereAverage", true)
    ),
    gradTPtr_(nullptr),
    epsilonPtr_(nullptr),
    gradTInterpPtr_(nullptr)
{
    // dsigma/dT: prefer the dictionary, else take the solver's own value so
    // the cloud cannot silently disagree with the Eulerian Marangoni term.
    if (!this->coeffs().readIfPresent("dSigmadT", dSigmadT_))
    {
        const IOdictionary* tpPtr =
            mesh.findObject<IOdictionary>("transportProperties");

        if (tpPtr)
        {
            // Marangoni_Constant is a dimensionedScalar in transportProperties
            dimensionedScalar mc("Marangoni_Constant", dimless, 0);
            if (tpPtr->readIfPresent("Marangoni_Constant", mc))
            {
                dSigmadT_ = mc.value();
            }
        }

        if (mag(dSigmadT_) < SMALL)
        {
            WarningInFunction
                << "dSigmadT is zero: neither the force dictionary nor "
                << "transportProperties/Marangoni_Constant supplied a value. "
                << "Thermocapillary migration will be inactive." << nl;
        }
    }

    Info<< "    Thermocapillary migration force:" << nl
        << "        dSigmadT        : " << dSigmadT_ << " N/m/K" << nl
        << "        mu*, k*         : " << viscosityRatio_
        << ", " << conductivityRatio_ << nl
        << "        Cmig            : " << Cmig_ << nl
        << "        sphere-averaged : " << Switch(useSphereAverage_) << nl;
}


template<class CloudType>
Foam::ThermocapillaryForce<CloudType>::ThermocapillaryForce
(
    const ThermocapillaryForce& tcf
)
:
    ParticleForce<CloudType>(tcf),
    gradTName_(tcf.gradTName_),
    epsilonName_(tcf.epsilonName_),
    dSigmadT_(tcf.dSigmadT_),
    viscosityRatio_(tcf.viscosityRatio_),
    conductivityRatio_(tcf.conductivityRatio_),
    epsMin_(tcf.epsMin_),
    Cmig_(tcf.Cmig_),
    useSphereAverage_(tcf.useSphereAverage_),
    gradTPtr_(nullptr),
    epsilonPtr_(nullptr),
    gradTInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermocapillaryForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const fvMesh& mesh = this->mesh();

        gradTPtr_ = mesh.findObject<volVectorField>(gradTName_);
        epsilonPtr_ = mesh.findObject<volScalarField>(epsilonName_);

        if (!gradTPtr_)
        {
            WarningInFunction
                << "Field " << gradTName_ << " not found in the registry; "
                << "thermocapillary force is zero." << nl;
            return;
        }

        if (!useSphereAverage_)
        {
            gradTInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    *gradTPtr_
                ).ptr()
            );
        }
    }
    else
    {
        gradTInterpPtr_.clear();
        gradTPtr_ = nullptr;
        epsilonPtr_ = nullptr;
    }
}


template<class CloudType>
Foam::forceSuSp Foam::ThermocapillaryForce<CloudType>::calcNonCoupled
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

    if (!gradTPtr_ || mag(dSigmadT_) < SMALL)
    {
        return value;
    }

    const label celli = p.cell();

    // No free surface to support a surface-tension gradient in solid metal.
    if (epsilonPtr_ && celli >= 0)
    {
        if (epsilonPtr_->primitiveField()[celli] < epsMin_)
        {
            return value;
        }
    }

    const scalar a = 0.5*p.d();

    vector gradT(Zero);

    if (useSphereAverage_)
    {
        label nT = 0;
        gradT = poreCloud::sphereAverage
        (
            this->mesh(),
            gradTPtr_->primitiveField(),
            celli,
            p.position(),
            a,
            nT
        );
    }
    else
    {
        gradT =
            gradTInterpPtr_().interpolate(p.coordinates(), p.currentTetIndices());
    }

    // F = 6 pi mu a U_YGB, with mu cancelling:
    //     F = -Cmig * 12 pi a^2 dSigmadT grad(T) / ((2+3mu*)(2+k*))
    const scalar denom =
        (2.0 + 3.0*viscosityRatio_)*(2.0 + conductivityRatio_);

    const scalar coeff =
       -Cmig_*12.0*constant::mathematical::pi*a*a*dSigmadT_/denom;

    // Explicit only - see the Sp caveat in LeenovKolinForce.
    value.Su() = coeff*gradT;

    return value;
}


// ************************************************************************* //
