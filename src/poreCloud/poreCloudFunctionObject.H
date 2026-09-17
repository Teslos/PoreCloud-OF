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
    Foam::functionObjects::poreCloudFunctionObject

Description
    Drives a one-way-coupled Lagrangian bubble cloud inside an unmodified
    laserbeamFoam run.

    Structurally this follows OpenFOAM's own
    functionObjects::icoUncoupledKinematicCloud, which demonstrates that a
    kinematic cloud can be owned and evolved by a functionObject rather than by
    the solver.  That is what keeps this project out of the solver: the cloud
    only reads fields, so no fork of laserbeamFoam and no change to UEqn.H is
    needed, and the validated azimuthal-field binary stays exactly as built.

    The stock class cannot be used directly because it does

    \verbatim
        mesh_.lookupObject<singlePhaseTransportModel>("transportProperties")
    \endverbatim

    and laserbeamFoam uses incompressibleTwoPhaseMixture, so that lookup aborts
    the run.  The carrier properties are therefore taken straight from the
    registry instead.  Note that mu is NOT registered by the solver -
    incompressibleTwoPhaseMixture::mu() returns a NO_REGISTER temporary - so it
    is formed here as rho*nu and refreshed every step, because nu is re-blended
    by the solver's Arrhenius viscosity model on every PIMPLE iteration.

    \heading Timing

    laserbeamFoam contains no explicit functionObjects().execute() call; the
    hook lives inside Time::run(), which is the while condition of the solver's
    time loop.  execute() therefore fires once per timestep at the TOP of the
    loop, before ++runTime, so runTime.value() is still the time of the
    just-completed step and every field - U, T, epsilon1, alpha.metal, J_MHD,
    gradT, rho, nu - holds converged end-of-PIMPLE values for that step.  This
    is the correct and self-consistent field set to advect bubbles through.

    Sub-timestep coupling is not available to a functionObject; the cloud
    sub-cycles internally against maxCo in its own solution dictionary, which
    is the appropriate granularity here.

Usage
    \verbatim
    functions
    {
        poreCloud
        {
            type            poreCloud;
            libs            ("libPoreCloud.so");
            executeControl  timeStep;
            executeInterval 1;
            writeControl    writeTime;
        }
    }
    \endverbatim

    Cloud sub-models are read from constant/<cloudName>Properties, i.e.
    constant/poreCloudProperties by default.

    \table
        Property   | Description                  | Required | Default
        type       | Type name: poreCloud         | yes      |
        U          | Velocity field name          | no       | U
        rho        | Density field name           | no       | rho
        nu         | Kinematic viscosity name     | no       | nu
        cloudName  | Name of the cloud            | no       | poreCloud
    \endtable

SourceFiles
    poreCloudFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_poreCloudFunctionObject_H
#define functionObjects_poreCloudFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "basicKinematicCloud.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                 Class poreCloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class poreCloudFunctionObject
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Acceleration due to gravity (constant/g, else zero)
        uniformDimensionedVectorField g_;

        //- Carrier density, registered by laserbeamFoam
        const volScalarField& rho_;

        //- Carrier kinematic viscosity, registered by laserbeamFoam
        const volScalarField& nu_;

        //- Carrier dynamic viscosity: rho*nu, owned here because the solver
        //  does not register mu
        volScalarField mu_;

        //- Carrier velocity
        const volVectorField& U_;

        //- Name of the cloud
        word cloudName_;

        //- The bubble cloud
        basicKinematicCloud cloud_;


    // Private Member Functions

        //- No copy construct
        poreCloudFunctionObject(const poreCloudFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const poreCloudFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("poreCloud");


    // Constructors

        //- Construct from Time and dictionary
        poreCloudFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~poreCloudFunctionObject() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Evolve the cloud
        virtual bool execute();

        //- The cloud writes itself via the registry
        virtual bool write();
};


} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
