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
    Foam::LeenovKolinForce

Description
    Electromagnetic exclusion force on a non-conducting inclusion (a gas pore)
    suspended in a current-carrying conducting melt.

    A gas void cannot carry current.  The current it displaces concentrates in
    the surrounding metal and, interacting with B, exerts a net force on the
    void - the Leenov-Kolin force.  For a sphere of volume V in an undisturbed
    current density J0:

    \verbatim
        F = C * V * (J0 x B)
    \endverbatim

    \heading The coefficient C

    C is dictionary-settable and intentionally so.  The parent solver contains
    two independently-derived implementations that disagree by a factor of two:

      - laserbeamFoam branch azimuthal-field, UEqn.H:185-202
            F/V = -1.5 * (1-alpha1) * (J0 x B)
        with J0 = fvc::average(fvc::interpolate(J_MHD))

      - laserbeamFoam branch fixes-solver, UEqn.H:61-76
            F/V = -0.75 * smoothed(LorentzForce) * (1-alpha1) * nearLiquid
        cited to Takahashi & Taniguchi (2003)

    Commit 0340af85 kept -3/2 and rejected -3/4; both coefficients appear in
    the inclusion-removal literature.  Since C scales the drift velocity
    linearly, and hence every pore-expulsion conclusion drawn from this solver,
    this class does not adjudicate: it defaults to -1.5 to match the live
    Eulerian implementation and lets the dictionary select otherwise, so the
    question can be settled by running both.

    \heading Why no double counting

    The Eulerian exclusionForce in UEqn.H is gated by (1 - alpha1).  Lagrangian
    bubbles never modify alpha.metal, so alpha1 = 1 at a bubble's location and
    that term is identically zero there.  Better still, J_MHD sampled at a
    bubble IS the undisturbed J0, because no void has been carved out of the
    conductivity field - which is cleaner than the Eulerian version, which must
    reconstruct J0 by smoothing J out of the void it does carve.

    \heading Finite-size sampling

    With useSphereAverage (default on), J0 and B are volume-averaged over the
    bubble sphere via poreCloud::sphereAverage, reproducing the integral the
    companion Eulerian tool PoreTracker performs over a resolved pore.  For
    35 um pores on a 6-9 um mesh this matters: the parcel spans ~10 cells and
    centroid interpolation is not representative.

Usage
    \verbatim
    particleForces
    {
        leenovKolin
        {
            J                J_MHD;     // registered by laserbeamFoam
            exclusionCoeff   -1.5;      // see discussion above
            useSphereAverage true;
        }
    }
    \endverbatim

SourceFiles
    LeenovKolinForce.C

\*---------------------------------------------------------------------------*/

#ifndef LeenovKolinForce_H
#define LeenovKolinForce_H

#include "ParticleForce.H"
#include "interpolation.H"
#include "poreCloudMagneticField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LeenovKolinForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LeenovKolinForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Name of the current density field
        const word JName_;

        //- Leenov-Kolin coefficient
        const scalar exclusionCoeff_;

        //- Volume-average the fields over the bubble sphere
        const bool useSphereAverage_;

        //- Reconstructed magnetic induction field (B is not registered)
        mutable autoPtr<poreCloud::magneticField> BFieldPtr_;

        //- Pointer to the registered current density, valid within a timestep
        const volVectorField* JPtr_;

        //- Centroid interpolators, used when useSphereAverage_ is false
        autoPtr<interpolation<vector>> JInterpPtr_;
        autoPtr<interpolation<vector>> BInterpPtr_;


public:

    //- Runtime type information
    TypeName("leenovKolin");


    // Constructors

        //- Construct from mesh
        LeenovKolinForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        LeenovKolinForce(const LeenovKolinForce& lkf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new LeenovKolinForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LeenovKolinForce() = default;


    // Member Functions

        // Access

            //- The Leenov-Kolin coefficient in use
            scalar exclusionCoeff() const noexcept { return exclusionCoeff_; }


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LeenovKolinForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
