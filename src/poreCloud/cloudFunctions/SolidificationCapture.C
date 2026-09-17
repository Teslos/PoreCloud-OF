/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see SolidificationCapture.H
\*---------------------------------------------------------------------------*/

#include "SolidificationCapture.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SolidificationCapture<CloudType>::SolidificationCapture
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    epsilonName_
    (
        this->coeffDict().template getOrDefault<word>("epsilonName", "epsilon1")
    ),
    captureThreshold_
    (
        this->coeffDict().template getOrDefault<scalar>("captureThreshold", 0.5)
    ),
    allowRemelt_
    (
        this->coeffDict().template getOrDefault<bool>("allowRemelt", false)
    ),
    nCaptured_(0)
{}


template<class CloudType>
Foam::SolidificationCapture<CloudType>::SolidificationCapture
(
    const SolidificationCapture<CloudType>& sc
)
:
    CloudFunctionObject<CloudType>(sc),
    epsilonName_(sc.epsilonName_),
    captureThreshold_(sc.captureThreshold_),
    allowRemelt_(sc.allowRemelt_),
    nCaptured_(sc.nCaptured_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SolidificationCapture<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    nCaptured_ = 0;
}


template<class CloudType>
void Foam::SolidificationCapture<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    label n = nCaptured_;
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        Info<< "    Bubbles entrapped by solidification front: " << n << nl;
    }
}


template<class CloudType>
bool Foam::SolidificationCapture<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    const typename parcelType::trackingData& td
)
{
    const fvMesh& mesh = this->owner().mesh();

    const volScalarField* epsPtr =
        mesh.template findObject<volScalarField>(epsilonName_);

    if (!epsPtr)
    {
        // Without a liquid fraction there is no front to detect; leave the
        // parcel alone rather than freezing the whole cloud.
        return true;
    }

    const label celli = p.cell();

    if (celli < 0 || celli >= mesh.nCells())
    {
        return true;
    }

    const scalar eps = epsPtr->primitiveField()[celli];

    if (p.active())
    {
        if (eps < captureThreshold_)
        {
            p.active(false);
            ++nCaptured_;
        }
    }
    else if (allowRemelt_ && eps >= captureThreshold_)
    {
        p.active(true);
    }

    // Keep the parcel: an entrapped pore is the result of interest, not
    // something to delete from the cloud.
    return true;
}


// ************************************************************************* //
