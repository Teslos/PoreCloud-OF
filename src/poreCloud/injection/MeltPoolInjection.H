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
    Foam::MeltPoolInjection

Description
    Continuous stochastic injection of gas bubbles into the liquid region of a
    melt pool.

    Unlike CellZoneInjection - from which the random-point-in-cell and parallel
    position-gathering machinery here is adapted - the injection region is not
    a static cellZone fixed at construction.  The melt pool translates with the
    laser and changes shape every timestep, so the candidate cell set is
    recomputed on each call to parcelsToInject() from a field predicate:

    \verbatim
        alpha.metal > alphaMin   (metal, not gas)
        epsilon1    > epsMin     (molten, not solid)
    \endverbatim

    which is exactly the solver's own definition of liquid metal (laserbeamFoam
    writes the same test as its "condition" field).

    Injection is a rate, not a one-shot seeding: the expected number of bubbles
    nucleated in a cell over a timestep is

    \verbatim
        n = rate [1/m^3/s] * V_cell * dt
    \endverbatim

    with the fractional remainder carried across cells so that low rates still
    produce the correct long-run count instead of rounding to zero everywhere.

    \heading Why this replaces the existing workflow

    Porosity was previously introduced by stopping the run, executing
    findPoreSites.py to write a setFieldsDict, carving resolved VoF spheres
    into alpha.metal with setFields, and restarting (see mhd-openfoam/injectPore
    and Allrun.repore).  That limits a study to a handful of hand-placed pores
    of one radius, all injected at a single instant.  This model injects a
    population with a size distribution, continuously, with no restart.

Usage
    \verbatim
    injectionModels
    {
        melt
        {
            type            meltPoolInjection;
            massTotal       0;          // unused with parcelBasisType fixed
            SOI             5e-4;
            duration        2e-3;
            parcelBasisType fixed;
            nParticle       1;          // one real bubble per parcel

            rate            1e12;       // nucleation rate [1/m^3/s]
            alphaName       alpha.metal;
            epsilonName     epsilon1;
            alphaMin        0.5;
            epsMin          0.5;
            maxParcels      20000;      // safety cap on cloud size

            // U0 omitted => bubbles are born with the local melt velocity
            sizeDistribution
            {
                type            normal;
                normalDistribution
                {
                    mu          35e-6;
                    sigma       10e-6;
                    minValue    10e-6;
                    maxValue    70e-6;
                }
            }
        }
    }
    \endverbatim

SourceFiles
    MeltPoolInjection.C

\*---------------------------------------------------------------------------*/

#ifndef MeltPoolInjection_H
#define MeltPoolInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class MeltPoolInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MeltPoolInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        //- Nucleation rate per unit liquid volume [1/m^3/s]
        const scalar rate_;

        //- Name of the metal phase fraction field
        const word alphaName_;

        //- Name of the liquid fraction field
        const word epsilonName_;

        //- Metal phase fraction threshold
        const scalar alphaMin_;

        //- Liquid fraction threshold
        const scalar epsMin_;

        //- Radial seeding: bias injection toward a moving axis so that the
        //  distribution of birth radii matches a measured one, instead of the
        //  default uniform-per-unit-liquid-volume.
        //
        //  This is not cosmetic. Uniform seeding places bubbles wherever there
        //  is liquid, which puts ~46% of them beyond 150 um of the keyhole -
        //  whereas resolved-VoF runs show 92-100% of real pores forming within
        //  100 um of it. Outcomes depend strongly on birth position (the
        //  up-fraction varies by a factor of ~4 across radius, and the sign of
        //  the trend differs between field configurations), so the seeding
        //  choice moves population statistics by tens of percentage points and
        //  can reorder configurations. See docs/specs.md.
        //
        //  Empty radialProfile_ leaves the original uniform behaviour.
        vector axisPoint_;
        vector axisVelocity_;

        //- Bin edges [m] and the target probability mass in each bin.
        scalarList radialEdges_;
        scalarList radialWeights_;

        //- Injection duration [s]
        const scalar duration_;

        //- Safety cap on the total number of parcels injected
        const label maxParcels_;

        //- Fixed initial velocity, if supplied; otherwise the local melt
        //  velocity is used
        const bool haveU0_;
        const vector U0_;

        //- Parcel size distribution
        const autoPtr<distributionModel> sizeDistribution_;

        //- Running total of parcels injected, for maxParcels_
        label nInjectedTotal_;

        //- Carry-over of the fractional parcel count between timesteps
        scalar fractionalCarry_;

        // Per-timestep injection lists, regenerated in parcelsToInject()

            List<vector> positions_;
            labelList injectorCells_;
            labelList injectorTetFaces_;
            labelList injectorTetPts_;
            scalarList diameters_;
            List<vector> velocities_;


    // Private Member Functions

        //- Rebuild the per-timestep injection lists from the current fields
        void generatePositions(const scalar dt);


public:

    //- Runtime type information
    TypeName("meltPoolInjection");


    // Constructors

        //- Construct from dictionary
        MeltPoolInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        MeltPoolInjection(const MeltPoolInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new MeltPoolInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MeltPoolInjection() = default;


    // Member Functions

        //- Injection lists are rebuilt every step, so nothing to remap
        virtual void updateMesh();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether injection of parcelI is allowed
            virtual bool validInjection(const label parcelI);
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeltPoolInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
