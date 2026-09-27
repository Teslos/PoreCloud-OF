/*---------------------------------------------------------------------------*\
    PoreCloud-OF   -   see poreCloudMagneticField.H
\*---------------------------------------------------------------------------*/

#include "poreCloudMagneticField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::poreCloud::magneticField::magneticField(const fvMesh& mesh)
:
    mesh_(mesh),
    active_(false),
    fieldType_("uniform"),
    B0_(0),
    omega_(0),
    rotationPlane_("xz"),
    alternating_(false),
    staticField_(1, 0, 0),
    coreRadius_(3.0e-5),
    axisPoint_(mesh.bounds().midpoint()),
    BPtr_(nullptr),
    builtStatic_(false)
{
    // transportProperties is registered by laserbeamFoam (createFields.H:90),
    // so it is available even though the B field itself is not.
    const IOdictionary* tpPtr =
        mesh_.findObject<IOdictionary>("transportProperties");

    if (!tpPtr)
    {
        WarningInFunction
            << "No registered 'transportProperties' dictionary found; "
            << "the magnetic field will be zero and the Leenov-Kolin force "
            << "inactive." << nl;
    }
    else if (tpPtr->found("MHD"))
    {
        const dictionary& mhd = tpPtr->subDict("MHD");

        active_ = mhd.getOrDefault<bool>("active", false);
        B0_ = mhd.getOrDefault<scalar>("B0", 0);

        const scalar freq = mhd.getOrDefault<scalar>("frequency", 0);
        omega_ = 2.0*constant::mathematical::pi*freq;

        fieldType_ = mhd.getOrDefault<word>("fieldType", "uniform");
        if (fieldType_ != "uniform" && fieldType_ != "azimuthal")
        {
            FatalErrorInFunction
                << "Unknown MHD fieldType " << fieldType_
                << "; expected uniform or azimuthal"
                << exit(FatalError);
        }

        rotationPlane_ = mhd.getOrDefault<word>("rotationPlane", "xz");
        if
        (
            rotationPlane_ != "xy"
         && rotationPlane_ != "xz"
         && rotationPlane_ != "yz"
        )
        {
            FatalErrorInFunction
                << "Unknown MHD rotationPlane " << rotationPlane_
                << "; expected xy, xz or yz"
                << exit(FatalError);
        }

        // Normalised at read time, exactly as the solver does, so that B0
        // alone sets the magnitude.
        alternating_ = mhd.getOrDefault<Switch>("alternating", false);
        staticField_ = mhd.getOrDefault<vector>("staticField", vector(1, 0, 0));
        const scalar magStatic = mag(staticField_);
        if (magStatic > SMALL)
        {
            staticField_ /= magStatic;
        }

        coreRadius_ = mhd.getOrDefault<scalar>("coreRadius", 3.0e-5);
        axisPoint_ =
            mhd.getOrDefault<point>("axisPoint", mesh.bounds().midpoint());
    }

    BPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "poreCloud:B",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "B",
                dimensionSet(1, 0, -2, 0, 0, -1, 0),
                Zero
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Private  * * * * * * * * * * * * * * * * //

void Foam::poreCloud::magneticField::buildAzimuthal()
{
    // B(x) = B0 * phiHat,  phiHat = (-(z-z0), 0, (x-x0)) / sqrt(r^2 + rc^2)
    // Matches laserbeamFoam UEqn.H azimuthal branch.
    volVectorField& B = BPtr_();
    vectorField& BI = B.primitiveFieldRef();
    const vectorField& cc = mesh_.C().primitiveField();

    const scalar x0 = axisPoint_.x();
    const scalar z0 = axisPoint_.z();
    const scalar rc2 = coreRadius_*coreRadius_;

    forAll(BI, celli)
    {
        const scalar dx = cc[celli].x() - x0;
        const scalar dz = cc[celli].z() - z0;
        const scalar denom = Foam::sqrt(dx*dx + dz*dz + rc2);

        BI[celli] = (B0_/denom)*vector(-dz, 0.0, dx);
    }

    B.correctBoundaryConditions();
}


void Foam::poreCloud::magneticField::buildUniform(const scalar t)
{
    volVectorField& B = BPtr_();

    vector Bvec(Zero);

    if (omega_ <= SMALL)
    {
        // Static DC field
        Bvec = B0_*staticField_;
    }
    else if (alternating_)
    {
        // AC: B oscillates along the fixed staticField_ axis rather than
        // rotating. rotationPlane_ is not consulted - mirrors UEqn.H.
        //
        // This branch is not cosmetic. Without it an AC case, which sets both
        // frequency and rotationPlane, fell through to the RMF branch below
        // and handed the cloud a rotating field. For staticField (0 1 0) that
        // turned an exactly-zero vertical exclusion force - (J x B)_y vanishes
        // when B is along y - into a spurious one of the same order as F_x,
        // making the run look like it produced vertical force when the
        // momentum equation applied none.
        Bvec = B0_*Foam::cos(omega_*t)*staticField_;
    }
    else
    {
        const scalar c = B0_*Foam::cos(omega_*t);
        const scalar s = B0_*Foam::sin(omega_*t);

        if (rotationPlane_ == "yz")
        {
            Bvec = vector(0.0, c, s);
        }
        else if (rotationPlane_ == "xy")
        {
            Bvec = vector(c, s, 0.0);
        }
        else // "xz"
        {
            Bvec = vector(c, 0.0, s);
        }
    }

    B.primitiveFieldRef() = Bvec;
    B.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::poreCloud::magneticField::update(const scalar t)
{
    if (!active_)
    {
        return;
    }

    if (fieldType_ == "azimuthal")
    {
        // Time-independent: build once.
        if (!builtStatic_)
        {
            buildAzimuthal();
            builtStatic_ = true;
        }
    }
    else
    {
        // DC is also time-independent, but rebuilding a uniform field is a
        // single assignment, so it is not worth special-casing.
        buildUniform(t);
    }
}


const Foam::volVectorField& Foam::poreCloud::magneticField::B() const
{
    return BPtr_();
}


void Foam::poreCloud::magneticField::info() const
{
    if (!active_)
    {
        Info<< "    MHD inactive - Leenov-Kolin force will be zero" << nl;
        return;
    }

    Info<< "    MHD field  : " << fieldType_
        << ", B0 = " << B0_ << " T" << nl;

    if (fieldType_ == "azimuthal")
    {
        Info<< "    axisPoint  : " << axisPoint_
            << ", coreRadius = " << coreRadius_ << " m" << nl;
    }
    else if (omega_ <= SMALL)
    {
        Info<< "    DC direction: " << staticField_ << " (normalised)" << nl;
    }
    else if (alternating_)
    {
        Info<< "    AC axis    : " << staticField_ << " (normalised)"
            << ", omega = " << omega_ << " rad/s" << nl;
    }
    else
    {
        Info<< "    RMF plane  : " << rotationPlane_
            << ", omega = " << omega_ << " rad/s" << nl;
    }
}


// ************************************************************************* //
