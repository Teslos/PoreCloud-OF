/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see poreCloudFunctionObject.H
\*---------------------------------------------------------------------------*/

#include "poreCloudFunctionObject.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(poreCloudFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        poreCloudFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Local Helpers  * * * * * * * * * * * * * //

namespace
{

//- Look up a registered field, failing with an actionable message rather than
//  OpenFOAM's generic "request for object ... failed".
//
//  Returns a pointer rather than a reference purely so that GCC's
//  -Wdangling-reference heuristic does not fire: it assumes a function
//  returning a reference may alias its arguments, and the field-name arguments
//  here are temporaries.  The referent is a registry object and outlives us.
template<class FieldType>
const FieldType* lookupRequired
(
    const Foam::fvMesh& mesh,
    const Foam::word& name,
    const char* what
)
{
    const FieldType* ptr = mesh.findObject<FieldType>(name);

    if (!ptr)
    {
        FatalErrorInFunction
            << "poreCloud requires the " << what << " field '" << name
            << "' to be registered by the solver, but it was not found." << nl
            << "This functionObject is designed for laserbeamFoam, which "
            << "registers U, rho, nu, T, epsilon1, alpha.metal, gradT and "
            << "J_MHD." << nl
            << "If the solver names it differently, set the corresponding "
            << "entry in the functionObject dictionary." << nl
            << Foam::exit(Foam::FatalError);
    }

    return ptr;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::poreCloudFunctionObject::poreCloudFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    g_(meshObjects::gravity::New(time_)),
    rho_
    (
        *lookupRequired<volScalarField>
        (
            mesh_,
            dict.getOrDefault<word>("rho", "rho"),
            "density"
        )
    ),
    nu_
    (
        *lookupRequired<volScalarField>
        (
            mesh_,
            dict.getOrDefault<word>("nu", "nu"),
            "kinematic viscosity"
        )
    ),
    mu_
    (
        IOobject
        (
            "mu",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*nu_
    ),
    U_
    (
        *lookupRequired<volVectorField>
        (
            mesh_,
            dict.getOrDefault<word>("U", "U"),
            "velocity"
        )
    ),
    cloudName_(dict.getOrDefault<word>("cloudName", "poreCloud")),
    cloud_
    (
        cloudName_,
        rho_,
        U_,
        mu_,
        g_
    )
{
    Info<< nl
        << "poreCloud: Lagrangian bubble cloud '" << cloudName_
        << "' attached to an unmodified solver run." << nl
        << "    sub-models read from constant/" << cloudName_
        << "Properties" << nl
        << "    gravity : " << g_.value() << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::poreCloudFunctionObject::read
(
    const dictionary& dict
)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::poreCloudFunctionObject::execute()
{
    // nu is re-blended by the solver's Arrhenius model every PIMPLE iteration
    // (nu = alpha1*nuA*exp(nuE/T) + (1-alpha1)*nuGas), so mu must be refreshed
    // rather than built once.
    mu_ = rho_*nu_;

    cloud_.evolve();

    return true;
}


bool Foam::functionObjects::poreCloudFunctionObject::write()
{
    // The cloud is an IOobject in the mesh registry, so runTime.write() in the
    // solver writes <time>/lagrangian/<cloudName>/ without help from here.
    return true;
}


// ************************************************************************* //
