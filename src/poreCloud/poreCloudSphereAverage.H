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

Namespace
    Foam::poreCloud

Description
    Finite-size field sampling for Lagrangian bubbles that span several cells.

    Melt-pool meshes in this study have ~6-9 um cells while the pores of
    interest are 10-70 um in diameter, i.e. up to ~10 cells across.  The stock
    ParticleForce models interpolate the carrier field at the parcel centroid,
    which is both inaccurate at that size and inconsistent with how the
    companion Eulerian tool (PoreTracker) measures the same quantity - it sums
    cell-by-cell over the whole resolved pore:

    \verbatim
        F = -k * sum_cells( V_c * f_c )            // PoreTracker/PoreExtract.C
    \endverbatim

    sphereAverage() reproduces that integral by breadth-first walking the
    cell-cell adjacency from the parcel's host cell and returning the
    VOLUME-WEIGHTED MEAN of the field over the cells whose centres lie within
    the bubble radius:

    \verbatim
        <f> = sum(V_c * f_c) / sum(V_c)
    \endverbatim

    A mean rather than PoreTracker's raw sum is returned deliberately:

      - The caller multiplies by the bubble's own analytic volume
        V = 4/3 pi a^3, so the product equals the integral while remaining
        independent of how many cells happened to be collected.
      - Normalising by the volume actually collected keeps the result an
        average over the cells that WERE seen, rather than silently losing a
        fraction of the force outright.

    \heading Parallel limitation - measured, not hypothetical

    fvMesh::cellCells() does not cross a processorPolyPatch, so in a decomposed
    run the walk stops at a subdomain boundary and a bubble straddling one
    samples only its local side.  Normalising bounds the damage but does not
    remove it: the sampled region is genuinely smaller and differently centred.

    Measured on the 0.2 T azimuthal melt-pool case: 181 bubbles at identical
    fixed positions (d/dx 3.1-5.4, ~50 cells per stencil), serial vs 18 ranks
    (scotch), matched one-to-one to within 1e-11 m.

    \verbatim
      coverage           serial   mean 0.980   min 0.893
                       18 ranks   mean 0.893   min 0.317

      |dF|/F by parallel coverage      n      mean      max
        >= 0.89  (serial floor)      130     1.79%    9.14%
        0.70 - 0.89                   24    13.05%   71.76%
        0.50 - 0.70                   21     9.45%   29.90%
        <  0.50                        6    31.20%   75.08%
      all bubbles                    181     5.15%   75.08%

      control: centroid sampling, no stencil    0.54%   15.11%
    \endverbatim

    The control isolates the cause: with no stencil at all the serial/parallel
    difference is 0.5% (the phiE Poisson solve converging differently under
    decomposition), so the rest is truncation.

    The filter that follows from this: in serial no bubble fell below coverage
    0.89 - that is the floor set by cell centres not tiling a sphere exactly,
    at this resolution - so a parallel bubble materially below it has been cut
    by a processor boundary rather than by discretisation.  Discarding those
    keeps 72% of the sample and brings the error back to serial-like levels,
    mean 5.15% -> 1.79%, max 75% -> 9.1%.  Note the floor is resolution
    dependent: take it from a serial run of the same case, do not assume 0.89.

    A proper fix needs remote stencil values, via mapDistribute (as
    extendedCentredCellToCellStencil does) or a layered halo exchange (as
    PoreExtract.C:322-382 does for its CCL).
      - When the radius falls below one cell it reduces exactly to the host
        cell value, so the same code path is valid for sub-grid bubbles.

    The number of cells collected is reported so the caller can emit it as a
    diagnostic (PoreTracker writes the equivalent as its NCells column).

SourceFiles
    poreCloudSphereAverage.H  (header-only template)

\*---------------------------------------------------------------------------*/

#ifndef poreCloudSphereAverage_H
#define poreCloudSphereAverage_H

#include "fvMesh.H"
#include "DynamicList.H"
#include "bitSet.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace poreCloud
{

//- Hard cap on the BFS, so a pathological radius cannot walk the whole mesh.
//  70 um / 6 um ~ 12 cells across => O(1e3) cells; 20000 is generous.
static constexpr label maxSphereCells = 20000;


//- Volume-weighted mean of an internal field over the cells of a sphere.
//
//  \param mesh     the mesh
//  \param fld      internal (primitive) field to average
//  \param celli    host cell of the parcel (always included)
//  \param centre   sphere centre, i.e. the parcel position
//  \param radius   sphere radius, i.e. the bubble radius
//  \param nCells   [out] number of cells collected
//  \param sumV     [out] total volume of the cells collected.  Compare with
//                  the analytic sphere volume to detect a truncated stencil -
//                  see coverage() below.
//
//  \return volume-weighted mean; the host cell value if radius is sub-cell.
template<class Type>
Type sphereAverage
(
    const fvMesh& mesh,
    const Field<Type>& fld,
    const label celli,
    const point& centre,
    const scalar radius,
    label& nCells,
    scalar& sumV
)
{
    nCells = 0;
    sumV = 0;

    if (celli < 0 || celli >= mesh.nCells())
    {
        return Zero;
    }

    const scalarField& V = mesh.V();
    const vectorField& C = mesh.C().primitiveField();

    // Sub-cell bubble: nothing to average, the host cell is the answer.
    // Cheap early-out that also covers radius <= 0.
    if (radius <= SMALL)
    {
        nCells = 1;
        sumV = mesh.V()[celli];
        return fld[celli];
    }

    const scalar r2 = radius*radius;

    const labelListList& cellCells = mesh.cellCells();

    bitSet visited(mesh.nCells());
    DynamicList<label> queue(256);

    // The host cell is unconditionally included: for a sub-cell bubble its
    // centre may lie further than radius from the parcel position, and
    // returning an empty average there would be wrong.
    visited.set(celli);
    queue.append(celli);

    Type sumVf = Zero;

    // Index-based BFS; queue grows in place (same idiom as PoreExtract::localCCL)
    for (label qi = 0; qi < queue.size(); ++qi)
    {
        const label c = queue[qi];

        sumVf += V[c]*fld[c];
        sumV  += V[c];
        ++nCells;

        if (queue.size() >= maxSphereCells)
        {
            break;
        }

        // Only expand through cells that are inside the sphere: the region is
        // convex, so this cannot miss an interior cell once the host is in.
        for (const label nbr : cellCells[c])
        {
            if (!visited.test(nbr) && magSqr(C[nbr] - centre) <= r2)
            {
                visited.set(nbr);
                queue.append(nbr);
            }
        }
    }

    if (sumV <= VSMALL)
    {
        nCells = 1;
        sumV = V[celli];
        return fld[celli];
    }

    return sumVf/sumV;
}


//- Convenience overload for callers that do not need the collected volume.
template<class Type>
Type sphereAverage
(
    const fvMesh& mesh,
    const Field<Type>& fld,
    const label celli,
    const point& centre,
    const scalar radius,
    label& nCells
)
{
    scalar sumV = 0;
    return sphereAverage(mesh, fld, celli, centre, radius, nCells, sumV);
}


//- Fraction of the bubble's analytic volume actually covered by the stencil.
//
//  Below 1 the walk was cut short - by a domain boundary, or (the case that
//  matters) by a processor boundary, since fvMesh::cellCells() does not cross a
//  processorPolyPatch.  A decomposed run therefore samples a smaller region
//  than the same bubble would sample in serial, and the sampled average is
//  biased accordingly.  Reported per bubble so affected samples can be
//  identified and filtered rather than silently trusted.
inline scalar coverage(const scalar sumV, const scalar radius)
{
    const scalar Vsphere =
        4.0/3.0*constant::mathematical::pi*radius*radius*radius;

    return Vsphere > VSMALL ? min(scalar(1), sumV/Vsphere) : scalar(1);
}


//- Representative cell length at a cell, used to report the finite-size ratio
//  d/dx so that violation of the point-particle assumption stays visible in
//  the output rather than implicit.
inline scalar cellLength(const fvMesh& mesh, const label celli)
{
    if (celli < 0 || celli >= mesh.nCells())
    {
        return SMALL;
    }

    return Foam::cbrt(mesh.V()[celli]);
}


} // End namespace poreCloud
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
