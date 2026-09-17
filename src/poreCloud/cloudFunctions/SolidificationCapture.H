/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    PoreCloud-OF
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later, as OpenFOAM.

Class
    Foam::SolidificationCapture

Description
    Freezes a bubble in place once the solidification front has passed it.

    A pore that the advancing solid front overtakes is entrapped: it stops
    moving and becomes the porosity that survives in the finished part.  That
    is the event the whole study is about - whether the magnetic field can get
    a bubble out of the melt before this happens.

    The mechanism is the one FreezeParticles uses:  setting active(false) makes
    KinematicParcel::move take the

    \verbatim
        else { p.stepFraction() += f; }
    \endverbatim

    branch instead of trackToFace(), so the parcel holds its position and cell
    exactly while continuing to be written to disk and to run function-object
    hooks.  That is a physically exact pinned pore, and it mirrors the Eulerian
    porePinning drag the solver already applies to gas in the mushy zone.

    Capture is tested on the liquid fraction at the bubble's cell:

    \verbatim
        epsilon1 < captureThreshold   =>   captured
    \endverbatim

    Capture is irreversible by default: remelting a captured bubble would
    require re-deciding where the front is relative to a body that is no longer
    being advected, and the study's interest is in final entrapment.  Set
    allowRemelt true to release bubbles when the metal melts again.

    No per-parcel capture time is stored.  The Active column written by
    PoreForceReport makes the transition visible, and the capture time is the
    first row where it drops to 0 - the same way PoreTracker's Python side
    infers pore death from the last row present.

Usage
    \verbatim
    cloudFunctions
    {
        solidificationCapture
        {
            epsilonName      epsilon1;
            captureThreshold 0.5;
            allowRemelt      false;
        }
    }
    \endverbatim

SourceFiles
    SolidificationCapture.C

\*---------------------------------------------------------------------------*/

#ifndef SolidificationCapture_H
#define SolidificationCapture_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SolidificationCapture Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SolidificationCapture
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Name of the liquid fraction field
        const word epsilonName_;

        //- Liquid fraction below which a bubble is entrapped
        const scalar captureThreshold_;

        //- Release captured bubbles if the metal remelts
        const bool allowRemelt_;

        //- Number of bubbles captured this evolve, for reporting
        label nCaptured_;


public:

    //- Runtime type information
    TypeName("solidificationCapture");


    // Constructors

        //- Construct from dictionary
        SolidificationCapture
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        SolidificationCapture(const SolidificationCapture<CloudType>& sc);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new SolidificationCapture<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SolidificationCapture() = default;


    // Member Functions

        //- Reset the per-evolve counter
        virtual void preEvolve(const typename parcelType::trackingData& td);

        //- Report captures
        virtual void postEvolve(const typename parcelType::trackingData& td);

        //- Freeze the parcel if the front has passed it
        virtual bool postMove
        (
            parcelType& p,
            const scalar dt,
            const point& position0,
            const typename parcelType::trackingData& td
        );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidificationCapture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
