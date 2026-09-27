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
    Foam::poreCloud::magneticField

Description
    Reconstructs the laserbeamFoam MHD magnetic induction field B.

    laserbeamFoam builds B as a function-local volVectorField inside UEqn.H
    ("Bfield_MHD", IOobject::NO_READ / NO_WRITE), so it is destroyed the moment
    the momentum predictor returns.  By the time a functionObject's execute()
    runs - at the top of the next time loop iteration - B no longer exists in
    the registry, unlike J_MHD, LorentzForce and exclusionForce which are
    registered and persist.

    This class therefore re-reads the MHD sub-dictionary of the registered
    "transportProperties" IOdictionary and rebuilds B using exactly the
    expressions in
    mhd-laserbeamfoam-solver/applications/solvers/laserbeamFoam/UEqn.H:6-83.

    Two field types are supported, matching the solver:

      fieldType azimuthal
        Static toroidal field about a vertical (y) axis through
        (axisPoint.x, *, axisPoint.z):

            B(x) = B0 * phiHat,
            phiHat = (-(z-z0), 0, (x-x0)) / sqrt(r^2 + rc^2)

        regularised on the axis by coreRadius rc.  Time-independent, so it is
        built once and reused.

      fieldType uniform
        DC when omega <= SMALL:      B = B0 * staticField   (normalised)
        RMF otherwise:               B rotates in rotationPlane at omega,
                                     one of xy, xz, yz.

    Note the solver normalises staticField at read time, so B0 alone sets the
    magnitude; an un-normalised entry such as (0 0 2) must not double |B|.

SourceFiles
    poreCloudMagneticField.C

\*---------------------------------------------------------------------------*/

#ifndef poreCloudMagneticField_H
#define poreCloudMagneticField_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace poreCloud
{

/*---------------------------------------------------------------------------*\
                       Class magneticField Declaration
\*---------------------------------------------------------------------------*/

class magneticField
{
    // Private Data

        const fvMesh& mesh_;

        //- MHD.active
        bool active_;

        //- MHD.fieldType: uniform | azimuthal
        word fieldType_;

        //- MHD.B0 [T]
        scalar B0_;

        //- 2*pi*MHD.frequency [rad/s]
        scalar omega_;

        //- MHD.rotationPlane: xy | xz | yz  (uniform + RMF only)
        word rotationPlane_;

        //- MHD.alternating: B oscillates along staticField_ instead of
        //  rotating in rotationPlane_. Must be read here, not just by the
        //  solver: this class rebuilds B independently, so a flag it does not
        //  know about silently gives the cloud a different field from the one
        //  the momentum equation applied.
        bool alternating_;

        //- MHD.staticField, normalised (uniform + DC only)
        vector staticField_;

        //- MHD.coreRadius [m] (azimuthal only)
        scalar coreRadius_;

        //- MHD.axisPoint [m] (azimuthal only)
        point axisPoint_;

        //- The reconstructed field
        autoPtr<volVectorField> BPtr_;

        //- True once an azimuthal (time-independent) field has been built
        bool builtStatic_;


    // Private Member Functions

        //- Populate B for the azimuthal field type
        void buildAzimuthal();

        //- Populate B for the uniform field type at time t
        void buildUniform(const scalar t);


public:

    // Constructors

        //- Construct from mesh, reading transportProperties/MHD
        explicit magneticField(const fvMesh& mesh);


    // Member Functions

        //- Is MHD active in this case?
        bool active() const noexcept { return active_; }

        //- Field type as read from the dictionary
        const word& fieldType() const noexcept { return fieldType_; }

        //- B0 as read from the dictionary [T]
        scalar B0() const noexcept { return B0_; }

        //- Update B for the given time.  A no-op for the azimuthal type after
        //  the first call, since that field is static.
        void update(const scalar t);

        //- Const access to the reconstructed field
        const volVectorField& B() const;

        //- Report the configuration to the log, once
        void info() const;
};


} // End namespace poreCloud
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
