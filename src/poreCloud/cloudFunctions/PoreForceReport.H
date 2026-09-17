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
    Foam::PoreForceReport

Description
    Per-bubble CSV output, schema-compatible with PoreTracker.

    The first fifteen columns are exactly PoreTracker's, in the same order and
    with the same meaning, so the existing analysis chain
    (PoreTracker-OF/scripts/PoreStats.py, mhd-openfoam/visualize_poretracker.py,
    plot_pore_trajectories_csv.py, compare_poretracker.py) reads this file
    without modification:

    \verbatim
    Time,PoreID,IsKeyhole,Volume_m3,Cx_m,Cy_m,Cz_m,Fx_N,Fy_N,Fz_N,
    BirthTime,BirthCx,BirthCy,BirthCz,NCells
    \endverbatim

    Fx_N..Fz_N carry the electromagnetic exclusion force, which is what
    PoreTracker puts there too (it accumulates
    -exclusionScale * sum(V_c * fLorentz_c) over the resolved pore).  The
    Lagrangian equivalent is recomputed here over the bubble sphere, so the two
    columns are directly comparable between a resolved-VoF pore and an injected
    bubble - which is the cross-validation the project depends on.

    Additional columns follow:

    \verbatim
    Ux,Uy,Uz            bubble velocity [m/s]
    d_m                 bubble diameter [m]
    dOverDx             d / cell length - the finite-size ratio
    Active              1 while free, 0 once entrapped by the front
    Uslipx,y,z          U_bubble - U_melt [m/s]
    Fbuoyx,y,z          buoyancy + weight, V (rho_melt - rho_gas) g [N]
    T_K                 local temperature [K]
    epsilon1            local liquid fraction
    \endverbatim

    \heading What is deliberately not reported

    Drag, added mass and lift are not broken out.  Those are evaluated inside
    KinematicParcel::calcVelocity from trackingData that is only valid mid-move
    and is stale by the time a cloud function object's postEvolve runs;
    recomputing them here would mean duplicating each drag correlation and
    risking silent divergence from what the cloud actually integrated.  Uslip
    is reported instead, from which any drag law can be evaluated in
    post-processing against the same slip velocity the solver used.

    The EM force IS recomputed here rather than intercepted, so exclusionCoeff
    below must match the value given to the leenovKolin force model.  A
    mismatch is reported at construction.

    \heading dOverDx

    The point-particle assumption behind Lagrangian tracking requires the
    bubble to be smaller than a cell.  At 35 um on a 6-9 um melt-pool mesh it
    is not, by roughly a factor of ten.  The sphere-average sampling is what
    makes that defensible, but the ratio is written out per bubble per step so
    the regime is visible in the data rather than buried in an assumption.

    \heading PoreID

    Encoded as origProc*1000000 + origId, which is unique and stable for the
    life of a parcel including across processor migration.  IsKeyhole is always
    0: an injected bubble is never the keyhole.

    BirthTime is exact (t - age).  BirthCx..z are recorded when a bubble is
    first seen by this object, so after a restart they refer to the restart
    position rather than the true nucleation site.

Usage
    \verbatim
    cloudFunctions
    {
        poreForceReport
        {
            writeControl    writeTime;    // or timeStep
            writeInterval   1;
            exclusionCoeff  -1.5;         // must match leenovKolin
            useSphereAverage true;        // must match leenovKolin
            J               J_MHD;
        }
    }
    \endverbatim

SourceFiles
    PoreForceReport.C

\*---------------------------------------------------------------------------*/

#ifndef PoreForceReport_H
#define PoreForceReport_H

#include "CloudFunctionObject.H"
#include "poreCloudMagneticField.H"
#include <fstream>
#include <memory>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PoreForceReport Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PoreForceReport
:
    public CloudFunctionObject<CloudType>
{
public:

    //- Number of scalars packed per bubble for the parallel gather
    static constexpr Foam::label nFields = 26;


private:

    // Private Data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Write cadence: writeTime | timeStep
        word writeControl_;

        //- Interval for writeControl timeStep
        label writeInterval_;

        //- Leenov-Kolin coefficient, must match the force model
        scalar exclusionCoeff_;

        //- Name of the current density field
        word JName_;

        //- Volume-average over the bubble sphere.  Must mirror the leenovKolin
        //  force model, else the reported force is not the integrated one.
        bool useSphereAverage_;

        //- Name of the liquid fraction field
        word epsilonName_;

        //- Output file name
        fileName csvName_;

        //- Reconstructed magnetic field
        mutable autoPtr<poreCloud::magneticField> BFieldPtr_;

        //- Output stream, master only
        mutable std::unique_ptr<std::ofstream> csvPtr_;

        //- Step counter for the timeStep cadence
        label stepCounter_;

        //- Birth positions, master only, keyed on PoreID
        std::map<label, point> birthPos_;


    // Private Member Functions

        //- Open the CSV and write its header (master only)
        void openCSV();

        //- Should a row be written at this call?
        bool shouldWrite() const;


public:

    //- Runtime type information
    TypeName("poreForceReport");


    // Constructors

        //- Construct from dictionary
        PoreForceReport
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PoreForceReport(const PoreForceReport<CloudType>& pfr);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new PoreForceReport<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PoreForceReport() = default;


    // Member Functions

        //- Gather and write the per-bubble rows
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PoreForceReport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
