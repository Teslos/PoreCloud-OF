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

Application
    poreCloudReplay

Description
    Replays the poreCloud Lagrangian bubble functionObject against the
    already-stored flow fields of a FINISHED laserbeamFoam run, instead of
    re-solving the CFD to try a different cloud configuration.

    Why this exists: the cloud is one-way coupled (it reads the melt fields,
    the melt never sees the parcels) and costs roughly 0.1% of the runtime
    the CFD itself costs. Re-running laserbeamFoam - 2.5 hours for one of the
    cases this was built for - to change a cloud parameter is therefore pure
    waste. This utility turns that 2.5 hours into seconds by driving the same
    poreCloud functionObject (unmodified - see src/poreCloud) against fields
    read back off disk.

    Interpolation scheme (the approximation this tool makes, and the thing
    to validate before trusting its output for anything quantitative):

    A finished case stores flow fields at a coarse cadence (typically every
    1e-5 s, set by controlDict's writeInterval on the original run) while the
    cloud - like the flow solve that produced these fields - wants to
    sub-cycle at something close to the original CFD timestep (typically
    ~1e-7 s, so ~87x finer than storage). Re-reading a ~30 MB timestep of
    fields from disk on every one of those ~87 sub-steps would cost as much
    as the CFD it replaces, defeating the entire point.

    So: this utility advances Time with its OWN timestep, taken from this
    case's controlDict exactly as any solver would (deltaT, adjustTimeStep).
    At every step it locates the two STORED snapshots bracketing the current
    replay time and linearly interpolates U, J_MHD, epsilon1, T and
    alpha.metal between them in place, cell-by-cell (see lerp() below). Only
    when the replay clock crosses into a new bracket is anything re-read from
    disk (see Bracket::moveTo()) - and even then, only the new upper snapshot
    needs reading, since the old upper snapshot becomes the new lower one on
    a monotonic forward replay. rho and nu are not stored on disk (createFields.H
    never writes them) and are re-derived every step from the interpolated
    alpha.metal, using the identical formulas laserbeamFoam itself uses - see
    the comments at the point of use below, and do not change them without
    re-checking applications/solvers/laserbeamFoam/{createFields.H,
    laserbeamFoam.C} in mhd-laserbeamfoam-solver (read-only reference, not
    part of this repo).

    LIMITATION: linear interpolation between snapshots 1e-5 s apart cannot
    reconstruct flow features that vary on the timescale of the original
    ~1e-7 s CFD step - transient turbulent/MHD fluctuations, fast oscillating
    fields (e.g. a rotating-field MHD case with a period comparable to the
    storage cadence), or anything that changes non-monotonically between two
    stored snapshots is smoothed out. Bubble trajectories computed by this
    replay are only as trustworthy as that smoothing is small relative to the
    bubble dynamics being studied - this must be checked case by case (e.g.
    by comparing a replay against a directly-coupled run over a short window)
    before treating replay output as ground truth rather than a fast first
    look.

    Carrier fields are exposed under exactly the names laserbeamFoam
    registers them under (U, rho, nu, T, epsilon1, alpha.metal, J_MHD,
    gradT), because poreCloud's functionObject, force models and injection
    model all find them purely by looking them up in the object registry by
    name (see src/poreCloud/poreCloudFunctionObject.H) - this utility does
    not, and must not, need to touch libPoreCloud at all. libPoreCloud.so is
    loaded the same way it is in a live run: dynamically, via controlDict's
    functions{ poreCloud { libs ("libPoreCloud.so"); ... } } entry, which
    this case must already have (see cases/meltPool-azimuthal or any prior
    laserbeamFoam+poreCloud case).

Usage
    Run like any OpenFOAM application, from a case directory that already
    contains a finished (or partially finished) laserbeamFoam run's time
    directories:

    \verbatim
        poreCloudReplay [-parallel]
    \endverbatim

    controlDict's startTime/endTime/deltaT/adjustTimeStep drive the replay
    clock as normal; startTime must fall within the range of stored
    snapshots (their bracket cannot be extrapolated). -parallel runs against
    a decomposed case exactly like any other OpenFOAM utility. Serial is
    preferred when available: poreCloud's sphere-averaged EM/thermocapillary
    sampling uses mesh::cellCells(), which does not cross a processor
    boundary, so a decomposed run undercounts the sampled volume for bubbles
    straddling a processor interface (see LeenovKolinForce.C and
    docs/specs.md in this repo for the measured error).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include <utility>

using namespace Foam;

// * * * * * * * * * * * * * * * * Local Helpers * * * * * * * * * * * * * //

namespace
{

//- The two on-disk snapshots currently bracketing the replay clock, for one
//  named field. Reloaded from disk only when the bracket changes - see
//  moveTo(). Fields are constructed with registerObject = false: the name
//  (e.g. "T") is already live in the mesh registry as the persistent,
//  interpolated target field that poreCloud looks up, and the registry does
//  not allow two objects under the same name at once.
template<class GeoField>
class Bracket
{
    const fvMesh& mesh_;
    const word fieldName_;
    label loIndex_;
    autoPtr<GeoField> lo_;
    autoPtr<GeoField> hi_;

    autoPtr<GeoField> readAt(const word& instance) const
    {
        return autoPtr<GeoField>::New
        (
            IOobject
            (
                fieldName_,
                instance,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false                   // registerObject = false, see above
            ),
            mesh_
        );
    }

public:

    Bracket(const fvMesh& mesh, const word& fieldName)
    :
        mesh_(mesh),
        fieldName_(fieldName),
        loIndex_(-1)
    {}

    //- Ensure lo()/hi() hold storedTimes[idx] / storedTimes[idx+1].
    void moveTo(const instantList& storedTimes, const label idx)
    {
        if (idx == loIndex_ && lo_ && hi_)
        {
            return;
        }

        if (hi_ && idx == loIndex_ + 1)
        {
            // Monotonic forward replay (the normal case): the old upper
            // snapshot becomes the new lower one - one disk read saved,
            // which is most of the point of caching at all.
            lo_ = std::move(hi_);
        }
        else
        {
            lo_ = readAt(storedTimes[idx].name());
        }

        hi_ = readAt(storedTimes[idx + 1].name());
        loIndex_ = idx;
    }

    const GeoField& lo() const { return *lo_; }
    const GeoField& hi() const { return *hi_; }
};


//- target = (1-w)*lo + w*hi, cell-by-cell (and face-by-face on boundaries).
template<class GeoField>
void lerp(GeoField& target, const GeoField& lo, const GeoField& hi, const scalar w)
{
    target = (1.0 - w)*lo + w*hi;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * Main * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Replay the poreCloud Lagrangian bubble functionObject against "
        "already-stored laserbeamFoam flow fields instead of re-solving "
        "the CFD. See the file header for the interpolation scheme and its "
        "limitations."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // ---------------------------------------------------------------------
    // Discover the stored snapshots (the finished run's field data) - this
    // is completely independent of the replay clock, which is our own and
    // driven by this case's controlDict deltaT/adjustTimeStep below.
    // ---------------------------------------------------------------------
    const instantList allTimes(runTime.times());

    DynamicList<instant> storedList(allTimes.size());
    forAll(allTimes, i)
    {
        // Guard against a directory that doesn't actually carry field data
        // (e.g. a bare 0/ written before the case was ever decomposed/run).
        if (isFile(runTime.path()/allTimes[i].name()/"alpha.metal"))
        {
            storedList.append(allTimes[i]);
        }
    }

    if (storedList.size() < 2)
    {
        FatalErrorInFunction
            << "Found " << storedList.size() << " stored time(s) with "
            << "alpha.metal under " << runTime.path() << "; need at least "
            << "two to interpolate between." << nl
            << "poreCloudReplay replays a FINISHED laserbeamFoam run - point "
            << "it at a case whose time directories already hold field data "
            << "(and, if decomposed, run with -parallel)." << nl
            << exit(FatalError);
    }

    const instantList storedTimes(storedList);

    Info<< "poreCloudReplay: " << storedTimes.size() << " stored snapshots, "
        << storedTimes.first().name() << " to " << storedTimes.last().name()
        << " s (storage cadence ~"
        << storedTimes[1].value() - storedTimes[0].value() << " s)" << nl
        << endl;

    if
    (
        runTime.value() < storedTimes.first().value() - SMALL
     || runTime.value() > storedTimes.last().value() + SMALL
    )
    {
        FatalErrorInFunction
            << "controlDict startTime (" << runTime.value() << ") is "
            << "outside the stored snapshot range ["
            << storedTimes.first().value() << ", "
            << storedTimes.last().value() << "]. Replay interpolates "
            << "between stored snapshots and cannot extrapolate beyond "
            << "them - set startTime inside that range." << nl
            << exit(FatalError);
    }

    // ---------------------------------------------------------------------
    // Mixture properties, read exactly as laserbeamFoam's createFields.H
    // does (immiscibleIncompressibleTwoPhaseMixture reads "rho" the same
    // way from the same subdicts) - see mhd-laserbeamfoam-solver
    // applications/solvers/laserbeamFoam/createFields.H (read-only
    // reference, not part of this repo, do not assume it stays in sync;
    // re-check it if this ever looks wrong).
    // ---------------------------------------------------------------------
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dictionary& metalDict = transportProperties.subDict("metal");
    const dictionary& gasDict = transportProperties.subDict("gas");

    const dimensionedScalar rho1("rho", dimDensity, metalDict);
    const dimensionedScalar rho2("rho", dimDensity, gasDict);

    // Temperature-dependent (Arrhenius) viscosity override, mirroring
    // laserbeamFoam.C's per-PIMPLE-iteration block:
    //   nu(T) = alpha1 * nuA*exp(nuE/max(T,300)) + (1-alpha1) * nuGas
    // active only when both nuA and nuE are given for the metal phase.
    const bool nuArrheniusActive =
        metalDict.found("nuA") && metalDict.found("nuE");

    dimensionedScalar nuA("nuA", dimViscosity, 0.0);
    dimensionedScalar nuE("nuE", dimTemperature, 0.0);
    dimensionedScalar nuGas("nu", dimViscosity, 0.0);
    dimensionedScalar nu1("nu", dimViscosity, 0.0);
    dimensionedScalar nu2("nu", dimViscosity, 0.0);

    if (nuArrheniusActive)
    {
        nuA = dimensionedScalar("nuA", dimViscosity, metalDict);
        nuE = dimensionedScalar("nuE", dimTemperature, metalDict);
        nuGas = dimensionedScalar("nu", dimViscosity, gasDict);

        Info<< "poreCloudReplay: Arrhenius viscosity active - nu(T) = "
            << nuA.value() << " * exp(" << nuE.value() << "/T)" << nl
            << endl;
    }
    else
    {
        // No Arrhenius entries: fall back to
        // incompressibleTwoPhaseMixture::calcNu()'s mass-weighted blend of
        // two constant (Newtonian) per-phase viscosities -
        //   nu = (a1*rho1*nu1 + (1-a1)*rho2*nu2) / (a1*rho1 + (1-a1)*rho2)
        // See src/transportModels/incompressible/incompressibleTwoPhase
        // Mixture/incompressibleTwoPhaseMixture.C::calcNu() in
        // mhd-laserbeamfoam-solver. Only covers the Newtonian
        // viscosityModel; a case using a different per-phase
        // viscosityModel (BirdCarreau, Casson, ...) is not reproduced here.
        nu1 = dimensionedScalar("nu", dimViscosity, metalDict);
        nu2 = dimensionedScalar("nu", dimViscosity, gasDict);

        Info<< "poreCloudReplay: no metal.nuA/nuE - using the mass-weighted "
            << "Newtonian nu blend instead of the Arrhenius model." << nl
            << endl;
    }

    // ---------------------------------------------------------------------
    // Registered carrier fields. Names and dimensions match createFields.H
    // exactly: poreCloud's functionObject, force models and injection model
    // all find these purely by looking them up in the registry by name.
    // ---------------------------------------------------------------------
    volScalarField alphaMetal
    (
        IOobject
        (
            "alpha.metal", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha.metal", dimless, 0.0)
    );

    volScalarField epsilon1
    (
        IOobject
        (
            "epsilon1", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("epsilon1", dimless, 0.0)
    );

    volScalarField T
    (
        IOobject
        (
            "T", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0.0)
    );

    volVectorField U
    (
        IOobject
        (
            "U", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimVelocity, vector::zero)
    );

    // Current density [A/m^2] - dimensions match createFields.H's J_MHD.
    volVectorField J_MHD
    (
        IOobject
        (
            "J_MHD", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "J_MHD", dimensionSet(0, -2, 0, 0, 0, 1, 0), vector::zero
        )
    );

    // Derived every step, never read from disk: laserbeamFoam never writes
    // rho (createFields.H: READ_IF_PRESENT, no AUTO_WRITE) or nu (owned by
    // the transport mixture, likewise unwritten).
    volScalarField rho
    (
        IOobject
        (
            "rho", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0.0)
    );

    volScalarField nu
    (
        IOobject
        (
            "nu", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("nu", dimViscosity, 0.0)
    );

    // gradT dimensions match createFields.H's gradT (K/m).
    volVectorField gradT
    (
        IOobject
        (
            "gradT", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
    );

    // ---------------------------------------------------------------------
    // Disk-snapshot brackets, one per field actually stored on disk.
    // ---------------------------------------------------------------------
    Bracket<volScalarField> alphaBracket(mesh, "alpha.metal");
    Bracket<volScalarField> epsilonBracket(mesh, "epsilon1");
    Bracket<volScalarField> TBracket(mesh, "T");
    Bracket<volVectorField> UBracket(mesh, "U");
    Bracket<volVectorField> JBracket(mesh, "J_MHD");

    label loIndex = 0;

    // Interpolate the on-disk carrier fields to time t, and re-derive rho,
    // nu and gradT from the result. Disk is touched only when t crosses
    // into a new bracket (Bracket::moveTo()).
    auto updateCarrierFields = [&](const scalar t)
    {
        const scalar tc =
            Foam::max
            (
                Foam::min(t, storedTimes.last().value()),
                storedTimes.first().value()
            );

        while
        (
            loIndex < storedTimes.size() - 2
         && storedTimes[loIndex + 1].value() < tc - SMALL
        )
        {
            ++loIndex;
        }

        alphaBracket.moveTo(storedTimes, loIndex);
        epsilonBracket.moveTo(storedTimes, loIndex);
        TBracket.moveTo(storedTimes, loIndex);
        UBracket.moveTo(storedTimes, loIndex);
        JBracket.moveTo(storedTimes, loIndex);

        const scalar tLo = storedTimes[loIndex].value();
        const scalar tHi = storedTimes[loIndex + 1].value();
        const scalar w = (tHi > tLo) ? (tc - tLo)/(tHi - tLo) : 0.0;

        lerp(alphaMetal, alphaBracket.lo(), alphaBracket.hi(), w);
        lerp(epsilon1, epsilonBracket.lo(), epsilonBracket.hi(), w);
        lerp(T, TBracket.lo(), TBracket.hi(), w);
        lerp(U, UBracket.lo(), UBracket.hi(), w);
        lerp(J_MHD, JBracket.lo(), JBracket.hi(), w);

        // rho = alpha1*rho1 + (1-alpha1)*rho2 - createFields.H, verbatim.
        rho = alphaMetal*rho1 + (scalar(1) - alphaMetal)*rho2;

        if (nuArrheniusActive)
        {
            // laserbeamFoam.C's Arrhenius override, verbatim (including the
            // alpha clamp to [0,1] and the T floor at 300 K).
            const volScalarField alphaClipped(min(max(alphaMetal, 0.0), 1.0));
            const volScalarField Tfloored
            (
                max(T, dimensionedScalar("Tfloor", dimTemperature, 300.0))
            );

            nu = alphaClipped*nuA*exp(nuE/Tfloored)
               + (scalar(1) - alphaClipped)*nuGas;
        }
        else
        {
            // incompressibleTwoPhaseMixture::calcNu()'s mass-weighted
            // blend, verbatim.
            const volScalarField alphaClipped(min(max(alphaMetal, 0.0), 1.0));

            nu =
            (
                alphaClipped*rho1*nu1
              + (scalar(1) - alphaClipped)*rho2*nu2
            )
           /(alphaClipped*rho1 + (scalar(1) - alphaClipped)*rho2);
        }

        // TEqn.H: gradT = fvc::grad(T), used by the thermocapillary force.
        gradT = fvc::grad(T);
    };

    // Time::run() fires functionObjects_.start() using the value already on
    // the clock the FIRST time it is called, and functionObjects_.execute()
    // using the (then-current) value on every call after that - in both
    // cases before this loop body runs (see Time::run(), and the ordering
    // note in src/poreCloud/poreCloudFunctionObject.H, which documents the
    // identical ordering in laserbeamFoam's own loop). So carrier fields
    // must already correspond to runTime.value() before the very first
    // "while (runTime.run())" check, which is why this call sits here
    // rather than as the first line inside the loop.
    updateCarrierFields(runTime.value());

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << endl;

        updateCarrierFields(runTime.value());

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
