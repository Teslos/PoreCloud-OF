/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    PoreCloud-OF
-------------------------------------------------------------------------------
Description
    Runtime-selection registration for the PoreCloud sub-models.

    Only the ...Type macros are used here.  The selection TABLES themselves are
    defined exactly once, by liblagrangianIntermediate in
    parcels/derived/basicKinematicParcel/makeBasicKinematicParcelSubmodels.C;
    invoking makeParticleForceModel(CloudType) again in this library would be a
    duplicate definition.  Because the tables are shared, these entries become
    selectable from a cloud properties dictionary with no change to OpenFOAM -
    the library only has to be loaded.

    The macros expand at GLOBAL scope (they open namespace Foam internally), so
    they must not be wrapped in a namespace block here.  Each also emits
    `typedef Foam::basicKinematicCloud::kinematicCloudType kinematicCloudType;`
    - repeating an identical typedef is legal, so the five uses coexist.

\*---------------------------------------------------------------------------*/

#include "basicKinematicCloud.H"

#include "LeenovKolinForce.H"
#include "ThermocapillaryForce.H"
#include "MeltPoolInjection.H"
#include "SolidificationCapture.H"
#include "PoreForceReport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forces
makeParticleForceModelType(LeenovKolinForce, basicKinematicCloud);
makeParticleForceModelType(ThermocapillaryForce, basicKinematicCloud);

// Injection
makeInjectionModelType(MeltPoolInjection, basicKinematicCloud);

// Cloud function objects
makeCloudFunctionObjectType(SolidificationCapture, basicKinematicCloud);
makeCloudFunctionObjectType(PoreForceReport, basicKinematicCloud);


// ************************************************************************* //
