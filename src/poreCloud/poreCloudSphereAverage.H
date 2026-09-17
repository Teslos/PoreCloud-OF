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
      - It degrades gracefully across processor boundaries.  fvMesh::cellCells()
        does not cross a processorPolyPatch, so the walk truncates at the local
        mesh; normalising by the volume actually collected keeps the result a
        correct average over the visible part instead of silently losing a
        fraction of the force.
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
    label& nCells
)
{
    nCells = 0;

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
    scalar sumV = 0;

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
        return fld[celli];
    }

    return sumVf/sumV;
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
