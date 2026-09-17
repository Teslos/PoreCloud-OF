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
    Foam::ThermocapillaryForce

Description
    Thermocapillary (Marangoni) migration of a gas bubble in a temperature
    gradient.

    A temperature gradient along a bubble surface produces a surface-tension
    gradient, which drives interfacial flow and propels the bubble.  The
    Young-Goldstein-Block (YGB) creeping-flow migration velocity is

    \verbatim
                        2                 a   dsigma
        U_YGB = - -------------------- * --- * ------ * grad(T)
                  (2 + 3 mu*)(2 + k*)    mu_m    dT
    \endverbatim

    with mu* = mu_bubble/mu_melt and k* = k_bubble/k_melt.  For a gas bubble in
    liquid metal both ratios are ~0, giving U_YGB = -(a/2 mu_m)(dsigma/dT) grad T.

    Since dsigma/dT < 0 for this alloy (Marangoni_Constant = -3.1e-4 N/m/K),
    U_YGB is parallel to grad(T): bubbles migrate toward hot metal, i.e. back
    toward the laser - which competes with the electromagnetic exclusion force
    that this study is trying to use to expel them.  Including it matters for
    that reason; it is not a small correction bolted on for completeness.

    \heading Expressing migration as a force

    The cloud integrates a force balance, so the migration velocity is supplied
    as the equivalent Stokes force that would produce it:

    \verbatim
        F = 6 pi mu_m a U_YGB
                                12 pi a^2      dsigma
          = - --------------------------- * ------ * grad(T)
              (2 + 3 mu*)(2 + k*)               dT
    \endverbatim

    Note mu_m cancels.  The 6 pi mu a factor assumes the rigid-sphere Stokes
    drag law; if the cloud's drag model differs (TomiyamaDrag does), the
    realised velocity will differ somewhat from U_YGB.  This is the usual
    engineering compromise and is scaled by the dictionary entry Cmig, which
    defaults to 1.

    \heading Gating

    The force is zero where the metal is not molten (epsilon1 below threshold):
    there is no free surface to support a surface-tension gradient in solid.

Usage
    \verbatim
    particleForces
    {
        thermocapillary
        {
            gradT           gradT;      // registered by laserbeamFoam TEqn.H
            epsilon1        epsilon1;   // liquid fraction
            dSigmadT        -3.1e-4;    // default: read from transportProperties
            viscosityRatio  0;          // mu_bubble / mu_melt
            conductivityRatio 0;        // k_bubble / k_melt
            epsMin          0.5;
            Cmig            1;
            useSphereAverage true;
        }
    }
    \endverbatim

SourceFiles
    ThermocapillaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef ThermocapillaryForce_H
#define ThermocapillaryForce_H

#include "ParticleForce.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ThermocapillaryForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermocapillaryForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Name of the temperature gradient field
        const word gradTName_;

        //- Name of the liquid fraction field
        const word epsilonName_;

        //- dsigma/dT [N/m/K]
        scalar dSigmadT_;

        //- mu_bubble / mu_melt
        const scalar viscosityRatio_;

        //- k_bubble / k_melt
        const scalar conductivityRatio_;

        //- Liquid-fraction threshold below which the force is switched off
        const scalar epsMin_;

        //- Empirical scale factor on the migration force
        const scalar Cmig_;

        //- Volume-average over the bubble sphere
        const bool useSphereAverage_;

        //- Pointers to registered fields, valid within a timestep
        const volVectorField* gradTPtr_;
        const volScalarField* epsilonPtr_;

        //- Centroid interpolator, used when useSphereAverage_ is false
        autoPtr<interpolation<vector>> gradTInterpPtr_;


public:

    //- Runtime type information
    TypeName("thermocapillary");


    // Constructors

        //- Construct from mesh
        ThermocapillaryForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ThermocapillaryForce(const ThermocapillaryForce& tcf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ThermocapillaryForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThermocapillaryForce() = default;


    // Member Functions

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
    #include "ThermocapillaryForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
