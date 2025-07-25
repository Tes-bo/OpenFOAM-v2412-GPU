/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "isoSurfaceTopo.H"
#include "polyMesh.H"
#include "volFields.H"
#include "edgeHashes.H"
#include "tetCell.H"
#include "DynamicField.H"
#include "syncTools.H"
#include "indirectPrimitivePatch.H"
#include "polyMeshTetDecomposition.H"
#include "foamVtkInternalMeshWriter.H"
#include "foamVtkLineWriter.H"
#include "foamVtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isoSurfaceBaseMethods.H"
defineIsoSurfaceInterpolateMethods(Foam::isoSurfaceTopo);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceTopo, 0);
}


// Get/set snapIndex (0, 1 or 2) at given position
// 0 = no snap
// 1 = snap to first edge end
// 2 = snap to second edge end
// NB: 4 lower bits left free for regular tet-cut information

#undef  SNAP_END_VALUE
#undef  SNAP_END_ENCODE
#define SNAP_END_ENCODE(pos, val)   (((val) << (4 + 2 * pos)))
#define SNAP_END_VALUE(pos, val)    (((val) >> (4 + 2 * pos)) & 0x3)


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check for tet values above/below given (iso) value
// Result encoded as an integer, with possible snapping information too
inline static int getTetCutIndex
(
    scalar p0,
    scalar p1,
    scalar p2,
    scalar p3,
    const scalar val,
    const bool doSnap
) noexcept
{
    int cutIndex
    (
        (p0 < val ? 1 : 0)  // point 0
      | (p1 < val ? 2 : 0)  // point 1
      | (p2 < val ? 4 : 0)  // point 2
      | (p3 < val ? 8 : 0)  // point 3
    );

    if (doSnap && cutIndex && cutIndex != 0xF)
    {
        // Not all below or all

        // Calculate distances (for snapping)
        p0 -= val; if (cutIndex & 1) p0 *= -1;
        p1 -= val; if (cutIndex & 2) p1 *= -1;
        p2 -= val; if (cutIndex & 4) p2 *= -1;
        p3 -= val; if (cutIndex & 8) p3 *= -1;

        // Add snap index into regular edge cut index
        // Snap to end if less than approx 1% of the distance.
        // - only valid if there is also a corresponding sign change
        #undef  ADD_SNAP_INDEX
        #define ADD_SNAP_INDEX(pos, d1, d2, idx1, idx2)                \
        switch (cutIndex & (idx1 | idx2))                              \
        {                                                              \
            case idx1 : /* first below, second above */                \
            case idx2 : /* first above, second below */                \
                cutIndex |= SNAP_END_ENCODE                            \
                (                                                      \
                    pos,                                               \
                    ((d1 * 100 < d2) ? 1 : (d2 * 100 < d1) ? 2 : 0)    \
                );                                                     \
                break;                                                 \
        }

        ADD_SNAP_INDEX(0, p0, p1, 1, 2);    // Edge 0: 0 -> 1
        ADD_SNAP_INDEX(1, p0, p2, 1, 4);    // Edge 1: 0 -> 2
        ADD_SNAP_INDEX(2, p0, p3, 1, 8);    // Edge 2: 0 -> 3
        ADD_SNAP_INDEX(3, p3, p1, 8, 2);    // Edge 3: 3 -> 1
        ADD_SNAP_INDEX(4, p1, p2, 2, 4);    // Edge 4: 1 -> 2
        ADD_SNAP_INDEX(5, p3, p2, 8, 4);    // Edge 5: 3 -> 2
        #undef ADD_SNAP_INDEX
    }

    return cutIndex;
}


// Append three labels to list.
// Filter out degenerate (eg, snapped) tris. Flip face as requested
inline static void appendTriLabels
(
    DynamicList<label>& verts,
    const label a,
    const label b,
    const label c,
    const bool flip  // Flip normals
)
{
    if (a != b && b != c && c != a)
    {
        verts.append(a);
        if (flip)
        {
            verts.append(c);
            verts.append(b);
        }
        else
        {
            verts.append(b);
            verts.append(c);
        }
    }
}


// Return point reference to mesh points or cell-centres
inline static const point& getMeshPointRef
(
    const polyMesh& mesh,
    const label pointi
)
{
    return
    (
        pointi < mesh.nPoints()
      ? mesh.points()[pointi]
      : mesh.cellCentres()[pointi - mesh.nPoints()]
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceTopo::tetCutAddressing::tetCutAddressing
(
    const label nCutCells,
    const bool useSnap,
    const bool useDebugCuts
)
:
    vertsToPointLookup_(12*nCutCells),
    snapVertsLookup_(0),

    pointToFace_(10*nCutCells),
    pointFromDiag_(10*nCutCells),

    pointToVerts_(10*nCutCells),
    cutPoints_(12*nCutCells),

    debugCutTets_(),
    debugCutTetsOn_(useDebugCuts)
{
    // Per cell: 5 pyramids cut, each generating 2 triangles

    // Per cell: number of intersected edges:
    //          - four faces cut so 4 mesh edges + 4 face-diagonal edges
    //          - 4 of the pyramid edges

    if (useSnap)
    {
        // Some, but not all, cells may have point snapping
        snapVertsLookup_.reserve(2*nCutCells);
    }
    if (debugCutTetsOn_)
    {
        debugCutTets_.reserve(6*nCutCells);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoSurfaceTopo::tetCutAddressing::clearDebug()
{
    debugCutTets_.clearStorage();
}


void Foam::isoSurfaceTopo::tetCutAddressing::clearDiagonal()
{
    pointToFace_.clearStorage();
    pointFromDiag_.clearStorage();
}


void Foam::isoSurfaceTopo::tetCutAddressing::clearHashes()
{
    vertsToPointLookup_.clear();
    snapVertsLookup_.clear();
}


Foam::label Foam::isoSurfaceTopo::tetCutAddressing::generatePoint
(
    label facei,
    bool edgeIsDiagonal,
    const int snapEnd,
    const edge& vertices
)
{
    // Generate new point, unless it already exists for edge
    // or corresponds to a snapped point (from another edge)

    label pointi = vertsToPointLookup_.lookup(vertices, -1);
    if (pointi == -1)
    {
        bool addNewPoint(true);

        const label snapPointi =
        (
            (snapEnd == 1) ? vertices.first()
          : (snapEnd == 2) ? vertices.second()
          : -1
        );

        if (snapPointi == -1)
        {
            // No snapped point
            pointi = pointToVerts_.size();
            pointToVerts_.append(vertices);
        }
        else
        {
            // Snapped point. No corresponding face or diagonal
            facei = -1;
            edgeIsDiagonal = false;

            pointi = snapVertsLookup_.lookup(snapPointi, -1);
            addNewPoint = (pointi == -1);
            if (addNewPoint)
            {
                pointi = pointToVerts_.size();
                snapVertsLookup_.insert(snapPointi, pointi);
                pointToVerts_.append(edge(snapPointi, snapPointi));
            }
        }

        if (addNewPoint)
        {
            pointToFace_.append(facei);
            pointFromDiag_.append(edgeIsDiagonal);
        }

        vertsToPointLookup_.insert(vertices, pointi);
    }

    return pointi;
}


bool Foam::isoSurfaceTopo::tetCutAddressing::generatePoints
(
    const label facei,
    const int tetCutIndex,
    const tetCell& tetLabels,

    // Per tet edge whether is face diag etc
    const FixedList<bool, 6>& edgeIsDiagonal
)
{
    bool flip(false);
    const label nCutPointsOld(cutPoints_.size());

    // Form the vertices of the triangles for each case
    switch (tetCutIndex & 0x0F)
    {
        case 0x00:
        case 0x0F:
        break;

        // Cut point 0
        case 0x0E: flip = true; [[fallthrough]];    // Point 0 above cut
        case 0x01:                                  // Point 0 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[0],
                    SNAP_END_VALUE(0, tetCutIndex),
                    tetLabels.edge(0)  // 0 -> 1
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[1],
                    SNAP_END_VALUE(1, tetCutIndex),
                    tetLabels.edge(1)  // 0 -> 2
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[2],
                    SNAP_END_VALUE(2, tetCutIndex),
                    tetLabels.edge(2)  // 0 -> 3
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
        }
        break;

        // Cut point 1
        case 0x0D: flip = true; [[fallthrough]];    // Point 1 above cut
        case 0x02:                                  // Point 1 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[0],
                    SNAP_END_VALUE(0, tetCutIndex),
                    tetLabels.edge(0)  // 0 -> 1
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[3],
                    SNAP_END_VALUE(3, tetCutIndex),
                    tetLabels.edge(3)  // 3 -> 1
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[4],
                    SNAP_END_VALUE(4, tetCutIndex),
                    tetLabels.edge(4)  // 1 -> 2
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
        }
        break;

        // Cut point 0/1 | 2/3
        case 0x0C: flip = true; [[fallthrough]];    // Point 0/1 above cut
        case 0x03:                                  // Point 0/1 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[1],
                    SNAP_END_VALUE(1, tetCutIndex),
                    tetLabels.edge(1)  // 0 -> 2
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[2],
                    SNAP_END_VALUE(2, tetCutIndex),
                    tetLabels.edge(2)  // 0 -> 3
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[3],
                    SNAP_END_VALUE(3, tetCutIndex),
                    tetLabels.edge(3)  // 3 -> 1
                )
            );
            const label cutD
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[4],
                    SNAP_END_VALUE(4, tetCutIndex),
                    tetLabels.edge(4)  // 1 -> 2
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
            appendTriLabels(cutPoints_, cutA, cutC, cutD, flip);
        }
        break;

        // Cut point 2
        case 0x0B: flip = true; [[fallthrough]];    // Point 2 above cut
        case 0x04:                                  // Point 2 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[1],
                    SNAP_END_VALUE(1, tetCutIndex),
                    tetLabels.edge(1)  // 0 -> 2
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[4],
                    SNAP_END_VALUE(4, tetCutIndex),
                    tetLabels.edge(4)  // 1 -> 2
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[5],
                    SNAP_END_VALUE(5, tetCutIndex),
                    tetLabels.edge(5)  // 3 -> 2
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
        }
        break;

        // Cut point 0/2 | 1/3
        case 0x0A: flip = true; [[fallthrough]];    // Point 0/2 above cut
        case 0x05:                                  // Point 0/2 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[0],
                    SNAP_END_VALUE(0, tetCutIndex),
                    tetLabels.edge(0)  // 0 -> 1
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[4],
                    SNAP_END_VALUE(4, tetCutIndex),
                    tetLabels.edge(4)  // 1 -> 2
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[5],
                    SNAP_END_VALUE(5, tetCutIndex),
                    tetLabels.edge(5)  // 3 -> 2
                )
            );
            const label cutD
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[2],
                    SNAP_END_VALUE(2, tetCutIndex),
                    tetLabels.edge(2)  // 0 -> 3
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
            appendTriLabels(cutPoints_, cutA, cutC, cutD, flip);
        }
        break;

        // Cut point 1/2 | 0/3
        case 0x09: flip = true; [[fallthrough]];    // Point 1/2 above cut
        case 0x06:                                  // Point 1/2 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[0],
                    SNAP_END_VALUE(0, tetCutIndex),
                    tetLabels.edge(0)  // 0 -> 1
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[3],
                    SNAP_END_VALUE(3, tetCutIndex),
                    tetLabels.edge(3)  // 3 -> 1
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[5],
                    SNAP_END_VALUE(5, tetCutIndex),
                    tetLabels.edge(5)  // 3 -> 2
                )
            );
            const label cutD
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[1],
                    SNAP_END_VALUE(1, tetCutIndex),
                    tetLabels.edge(1)  // 0 -> 2
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
            appendTriLabels(cutPoints_, cutA, cutC, cutD, flip);
        }
        break;

        // Cut point 3
        case 0x07: flip = true; [[fallthrough]];    // Point 3 above cut
        case 0x08:                                  // Point 3 below cut
        {
            const label cutA
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[2],
                    SNAP_END_VALUE(2, tetCutIndex),
                    tetLabels.edge(2)  // 0 -> 3
                )
            );
            const label cutB
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[5],
                    SNAP_END_VALUE(5, tetCutIndex),
                    tetLabels.edge(5)  // 3 -> 2
                )
            );
            const label cutC
            (
                generatePoint
                (
                    facei,
                    edgeIsDiagonal[3],
                    SNAP_END_VALUE(3, tetCutIndex),
                    tetLabels.edge(3)  // 3 -> 1
                )
            );

            appendTriLabels(cutPoints_, cutA, cutB, cutC, flip);
        }
        break;
    }

    const bool added(nCutPointsOld != cutPoints_.size());

    if (added && debugCutTetsOn_)
    {
        debugCutTets_.append(tetLabels.shape());
    }

    return added;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Requires mesh_, tetBasePtIs
void Foam::isoSurfaceTopo::generateTriPoints
(
    const label celli,
    const bool isTet,
    const labelList& tetBasePtIs,
    tetCutAddressing& tetCutAddr
) const
{
    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const cell& cFaces = mesh_.cells()[celli];
    const bool doSnap = this->snap();

    if (isTet)
    {
        // For tets don't do cell-centre decomposition, just use the
        // tet points and values

        const label facei = cFaces[0];
        const face& f0 = faces[facei];

        // Get the other point from f1. Tbd: check if not duplicate face
        // (ACMI / ignoreBoundaryFaces_).
        const face& f1 = faces[cFaces[1]];
        label apexi = -1;
        forAll(f1, fp)
        {
            apexi = f1[fp];
            if (!f0.found(apexi))
            {
                break;
            }
        }

        const label p0 = f0[0];
        label p1 = f0[1];
        label p2 = f0[2];

        if (faceOwner[facei] == celli)
        {
            std::swap(p1, p2);
        }

        const tetCell tetLabels(p0, p1, p2, apexi);
        const int tetCutIndex
        (
            getTetCutIndex
            (
                pVals_[p0],
                pVals_[p1],
                pVals_[p2],
                pVals_[apexi],
                iso_,
                doSnap
            )
        );

        tetCutAddr.generatePoints
        (
            facei,
            tetCutIndex,
            tetLabels,
            FixedList<bool, 6>(false)  // Not face diagonal
        );
    }
    else
    {
        for (const label facei : cFaces)
        {
            if
            (
               !mesh_.isInternalFace(facei)
             && ignoreBoundaryFaces_.test(facei-mesh_.nInternalFaces())
            )
            {
                continue;
            }

            const face& f = faces[facei];

            label fp0 = tetBasePtIs[facei];

            // Fallback
            if (fp0 < 0)
            {
                fp0 = 0;
            }

            const label p0 = f[fp0];
            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); ++i)
            {
                label p1 = f[fp];
                fp = f.fcIndex(fp);
                label p2 = f[fp];

                FixedList<bool, 6> edgeIsDiagonal(false);
                if (faceOwner[facei] == celli)
                {
                    std::swap(p1, p2);
                    if (i != 2) edgeIsDiagonal[1] = true;
                    if (i != f.size()-1) edgeIsDiagonal[0] = true;
                }
                else
                {
                    if (i != 2) edgeIsDiagonal[0] = true;
                    if (i != f.size()-1) edgeIsDiagonal[1] = true;
                }

                const tetCell tetLabels(p0, p1, p2, mesh_.nPoints()+celli);
                const int tetCutIndex
                (
                    getTetCutIndex
                    (
                        pVals_[p0],
                        pVals_[p1],
                        pVals_[p2],
                        cVals_[celli],
                        iso_,
                        doSnap
                    )
                );

                tetCutAddr.generatePoints
                (
                    facei,
                    tetCutIndex,
                    tetLabels,
                    edgeIsDiagonal
                );
            }
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::isoSurfaceTopo::triangulateOutside
(
    const bool filterDiag,
    const primitivePatch& pp,
    const boolUList& pointFromDiag,
    const labelUList& pointToFace,
    const label cellID,

    // outputs
    DynamicList<face>& compactFaces,
    DynamicList<label>& compactCellIDs
)
{
    // We can form pockets:
    // - 1. triangle on face
    // - 2. multiple triangles on interior (from diag edges)
    // - the edge loop will be pocket since it is only the diag
    //   edges that give it volume?

    // Retriangulate the exterior loops
    const labelListList& edgeLoops = pp.edgeLoops();
    const labelList& mp = pp.meshPoints();

    for (const labelList& loop : edgeLoops)
    {
        if (loop.size() > 2)
        {
            face& f = compactFaces.emplace_back(loop.size());

            label fpi = 0;
            forAll(f, i)
            {
                const label pointi = mp[loop[i]];
                if (filterDiag && pointFromDiag[pointi])
                {
                    const label prevPointi = mp[loop[loop.fcIndex(i)]];
                    if
                    (
                        pointFromDiag[prevPointi]
                     && (pointToFace[pointi] != pointToFace[prevPointi])
                    )
                    {
                        f[fpi++] = pointi;
                    }
                    else
                    {
                        // Filter out diagonal point
                    }
                }
                else
                {
                    f[fpi++] = pointi;
                }
            }

            if (fpi > 2)
            {
                f.resize(fpi);
            }
            else
            {
                // Keep original face
                forAll(f, i)
                {
                    const label pointi = mp[loop[i]];
                    f[i] = pointi;
                }
            }
            compactCellIDs.append(cellID);
        }
    }
}


void Foam::isoSurfaceTopo::removeInsidePoints
(
    Mesh& s,
    const bool filterDiag,

    // inputs
    const boolUList& pointFromDiag,
    const labelUList& pointToFace,
    const labelUList& start,                // Per cell the starting triangle

    // outputs
    DynamicList<label>& compactCellIDs      // Per returned tri the cellID
)
{
    const pointField& points = s.points();

    compactCellIDs.clear();
    compactCellIDs.reserve(s.size()/4);

    DynamicList<face> compactFaces(s.size()/4);

    for (label celli = 0; celli < start.size()-1; ++celli)
    {
        // Triangles for the current cell

        const label nTris = start[celli+1]-start[celli];

        if (nTris)
        {
            const primitivePatch pp
            (
                SubList<face>(s, nTris, start[celli]),
                points
            );

            triangulateOutside
            (
                filterDiag,
                pp,
                pointFromDiag,
                pointToFace,
                celli,

                compactFaces,
                compactCellIDs
            );
        }
    }

    s.swapFaces(compactFaces);  // Use new faces
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceTopo::isoSurfaceTopo
(
    const polyMesh& mesh,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const isoSurfaceParams& params,
    const bitSet& ignoreCells
)
:
    isoSurfaceBase(mesh, cellValues, pointValues, iso, params)
{
    // The cell cut type
    List<cutType> cellCutType_(mesh.nCells(), cutType::UNVISITED);

    // Time description (for debug output)
    const word timeDesc(word::printf("%08d", mesh_.time().timeIndex()));

    if (debug)
    {
        Pout<< "isoSurfaceTopo:" << nl
            << "    cell min/max  : " << minMax(cVals_) << nl
            << "    point min/max : " << minMax(pVals_) << nl
            << "    isoValue      : " << iso << nl
            << "    filter        : "
            << isoSurfaceParams::filterNames[params.filter()] << nl
            << "    mesh span     : " << mesh.bounds().mag() << nl
            << "    ignoreCells   : " << ignoreCells.count()
            << " / " << cVals_.size() << nl
            << endl;
    }

    this->ignoreCyclics();

    label nBlockedCells = 0;

    // Mark ignoreCells as blocked
    nBlockedCells += blockCells(cellCutType_, ignoreCells);

    // Mark cells outside bounding box as blocked
    nBlockedCells +=
        blockCells(cellCutType_, params.getClipBounds(), volumeType::OUTSIDE);

    // Adjusted tet base points to improve tet quality
    labelList tetBasePtIs
    (
        polyMeshTetDecomposition::adjustTetBasePtIs(mesh_, debug)
    );


    // Determine cell cuts
    const label nCutCells = calcCellCuts(cellCutType_);

    if (debug)
    {
        Pout<< "isoSurfaceTopo : candidate cells cut "
            << nCutCells
            << " blocked " << nBlockedCells
            << " total " << mesh_.nCells() << endl;
    }

    if (debug && isA<fvMesh>(mesh))
    {
        const auto& fvmesh = refCast<const fvMesh>(mesh);

        volScalarField debugField
        (
            IOobject
            (
                "isoSurfaceTopo.cutType",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvmesh,
            dimensionedScalar(dimless)
        );

        auto& debugFld = debugField.primitiveFieldRef();

        forAll(cellCutType_, celli)
        {
            debugFld[celli] = cellCutType_[celli];
        }

        Info<< "Writing cut types: " << debugField.objectRelPath() << nl;
        debugField.write();
    }

    // Additional debugging
    if (debug & 8)
    {
        // Write debug cuts cells in VTK format
        {
            constexpr uint8_t realCut(cutType::CUT | cutType::TETCUT);
            labelList debugCutCells(nCutCells, Zero);

            label nout = 0;
            forAll(cellCutType_, celli)
            {
                if ((cellCutType_[celli] & realCut) != 0)
                {
                    debugCutCells[nout] = celli;
                    ++nout;
                    if (nout >= nCutCells) break;
                }
            }

            // The mesh subset cut
            vtk::vtuCells vtuCells;
            vtuCells.reset(mesh_, debugCutCells);

            vtk::internalMeshWriter writer
            (
                mesh_,
                vtuCells,
                fileName
                (
                    mesh_.time().globalPath()
                  / ("isoSurfaceTopo." + timeDesc + "-cutCells")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();
            writer.writeCellData("cutField", cVals_);

            // PointData
            writer.beginPointData();
            writer.writePointData("cutField", pVals_);

            Info<< "isoSurfaceTopo : (debug) wrote "
                << returnReduce(nCutCells, sumOp<label>())
                << " cut cells: "
                << writer.output().name() << nl;
        }
    }


    tetCutAddressing tetCutAddr
    (
        nCutCells,
        this->snap(),
        (debug & 8)  // Enable debug tets
    );

    labelList startTri(mesh_.nCells()+1, Zero);
    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        startTri[celli] = tetCutAddr.nFaces();
        if ((cellCutType_[celli] & cutType::ANYCUT) != 0)
        {
            generateTriPoints
            (
                celli,
                // Same as tetMatcher::test(mesh_, celli),
                bool(cellCutType_[celli] & cutType::TETCUT),

                tetBasePtIs,
                tetCutAddr
            );
        }
    }
    startTri.back() = tetCutAddr.nFaces();

    // Information not needed anymore:
    tetBasePtIs.clear();
    tetCutAddr.clearHashes();


    // From list of vertices -> triangular faces
    faceList allTriFaces(startTri.back());
    {
        auto& verts = tetCutAddr.cutPoints();

        label verti = 0;
        for (face& f : allTriFaces)
        {
            f.resize(3);
            f[0] = verts[verti++];
            f[1] = verts[verti++];
            f[2] = verts[verti++];
        }
        verts.clearStorage();  // Not needed anymore
    }


    // The cells cut by the triangular faces
    meshCells_.resize(startTri.back());
    for (label celli = 0; celli < startTri.size()-1; ++celli)
    {
        // All triangles for the current cell
        labelList::subList
        (
            meshCells_,
            (startTri[celli+1] - startTri[celli]),
            startTri[celli]
        ) = celli;
    }


    pointToVerts_.transfer(tetCutAddr.pointToVerts());

    pointField allTriPoints
    (
        this->interpolateTemplate
        (
            mesh_.cellCentres(),
            mesh_.points()
        )
    );


    // Assign to MeshedSurface
    static_cast<Mesh&>(*this) = Mesh
    (
        std::move(allTriPoints),
        std::move(allTriFaces),
        surfZoneList()  // zones not required (one zone)
    );

    if (debug)
    {
        Pout<< "isoSurfaceTopo : generated "
            << Mesh::size() << " triangles "
            << Mesh::points().size() << " points" << endl;
    }

    // Write debug triangulated surface
    if ((debug & 8) && (params.filter() != filterType::NONE))
    {
        const Mesh& s = *this;

        vtk::surfaceWriter writer
        (
            s.points(),
            s,
            fileName
            (
                mesh_.time().globalPath()
              / ("isoSurfaceTopo." + timeDesc + "-triangles")
            )
        );

        writer.writeGeometry();

        // CellData
        writer.beginCellData();
        writer.writeProcIDs();
        writer.write("cellID", meshCells_);

        // PointData
        writer.beginPointData();
        {
            // NB: may have non-compact surface points
            // --> use points().size() not nPoints()!

            labelList pointStatus(s.points().size(), Zero);

            forAll(pointToVerts_, i)
            {
                const edge& verts = pointToVerts_[i];
                if (verts.first() == verts.second())
                {
                    // Duplicate index (ie, snapped)
                    pointStatus[i] = 1;
                }
                if (tetCutAddr.pointFromDiag().test(i))
                {
                    // Point on triangulation diagonal
                    pointStatus[i] = -1;
                }
            }

            writer.write("point-status", pointStatus);
        }

        Info<< "isoSurfaceTopo : (debug) wrote "
            << returnReduce(s.size(), sumOp<label>())
            << " triangles : "
            << writer.output().name() << nl;
    }


    // Now:
    // - generated faces and points are assigned to *this
    // - per point we know:
    //  - pointOnDiag: whether it is on a face-diagonal edge
    //  - pointToFace: from what pyramid (cell+face) it was produced
    //    (note that the pyramid faces are shared between multiple mesh faces)
    //  - pointToVerts_ : originating mesh vertex or cell centre

    if (params.filter() == filterType::NONE)
    {
        // Compact out unused (snapped) points
        if (this->snap())
        {
            Mesh& s = *this;

            labelList pointMap;          // Back to original point
            s.compactPoints(pointMap);   // Compact out unused points
            pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointMap)();
        }
    }
    else
    {
        // Initial filtering

        Mesh& s = *this;

        // Triangulate outside
        // (filter edges to cell centres and optionally face diagonals)
        DynamicList<label> compactCellIDs;  // Per tri the cell

        removeInsidePoints
        (
            *this,
            // Filter face diagonal
            (
                params.filter() == filterType::DIAGCELL
             || params.filter() == filterType::NONMANIFOLD
            ),
            tetCutAddr.pointFromDiag(),
            tetCutAddr.pointToFace(),
            startTri,
            compactCellIDs
        );

        labelList pointMap;          // Back to original point
        s.compactPoints(pointMap);   // Compact out unused points

        pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointMap)();
        meshCells_.transfer(compactCellIDs);

        if (debug)
        {
            Pout<< "isoSurfaceTopo :"
                " after removing cell centre and face-diag triangles: "
                << Mesh::size() << " faces "
                << Mesh::points().size() << " points"
                << endl;
        }
    }

    // Diagonal filter information not needed anymore
    tetCutAddr.clearDiagonal();


    // For more advanced filtering (eg, removal of open edges)
    // need the boundary and other 'protected' points

    bitSet isProtectedPoint;
    if
    (
        (params.filter() == filterType::NONMANIFOLD)
     || tetCutAddr.debugCutTetsOn()
    )
    {
        // Mark points on mesh outside as 'protected'
        // - never erode these edges

        isProtectedPoint.resize(mesh_.nPoints());

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            ++facei
        )
        {
            isProtectedPoint.set(mesh_.faces()[facei]);
        }

        // Include faces that would be exposed from mesh subset
        if (nBlockedCells)
        {
            const labelList& faceOwn = mesh_.faceOwner();
            const labelList& faceNei = mesh_.faceNeighbour();

            for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
            {
                // If only one cell is blocked, the face corresponds
                // to an exposed subMesh face
                if
                (
                    (cellCutType_[faceOwn[facei]] == cutType::BLOCKED)
                 != (cellCutType_[faceNei[facei]] == cutType::BLOCKED)
                )
                {
                    isProtectedPoint.set(mesh_.faces()[facei]);
                }
            }
        }
    }

    // Initial cell cut information not needed anymore
    cellCutType_.clear();


    // Additional debugging
    if (tetCutAddr.debugCutTetsOn())
    {
        // Write debug cut tets in VTK format
        {
            const auto& debugCuts = tetCutAddr.debugCutTets();

            // The TET shapes, using the mesh_ points information
            vtk::vtuCells vtuCells;
            vtuCells.resetShapes(debugCuts);

            // Use all points and all cell-centres
            vtuCells.setNumPoints(mesh_.nPoints());
            vtuCells.addPointCellLabels(identity(mesh_.nCells()));

            vtk::internalMeshWriter writer
            (
                mesh_,
                vtuCells,
                fileName
                (
                    mesh_.time().globalPath()
                  / ("isoSurfaceTopo." + timeDesc + "-cutTets")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();

            // Quality of the cut tets
            {
                Field<scalar> cutTetQuality(debugCuts.size());
                forAll(cutTetQuality, teti)
                {
                    cutTetQuality[teti] = tetPointRef
                    (
                        getMeshPointRef(mesh_, debugCuts[teti][0]),
                        getMeshPointRef(mesh_, debugCuts[teti][1]),
                        getMeshPointRef(mesh_, debugCuts[teti][2]),
                        getMeshPointRef(mesh_, debugCuts[teti][3])
                    ).quality();
                }
                writer.writeCellData("tetQuality", cutTetQuality);
            }

            // PointData
            if (this->snap())
            {
                writer.beginPointData();

                labelList pointStatus(vtuCells.nFieldPoints(), Zero);

                for (const edge& verts : pointToVerts_)
                {
                    if (verts.first() == verts.second())
                    {
                        // Duplicate index (ie, snapped)
                        pointStatus[verts.first()] = 1;
                    }
                }

                writer.writePointData("point-status", pointStatus);
            }

            Info<< "isoSurfaceTopo : (debug) wrote "
                << returnReduce(debugCuts.size(), sumOp<label>())
                << " cut tets: "
                << writer.output().name() << nl;
        }

        // Determining open edges. Same logic as used later...

        labelHashSet openEdgeIds(0);

        {
            const Mesh& s = *this;

            const labelList& mp = s.meshPoints();
            const edgeList& surfEdges = s.edges();
            const labelListList& edgeFaces = s.edgeFaces();
            openEdgeIds.reserve(s.size());

            forAll(edgeFaces, edgei)
            {
                const labelList& eFaces = edgeFaces[edgei];
                if (eFaces.size() == 1)
                {
                    // Open edge (not originating from a boundary face)

                    const edge& e = surfEdges[edgei];
                    const edge& verts0 = pointToVerts_[mp[e.first()]];
                    const edge& verts1 = pointToVerts_[mp[e.second()]];

                    if
                    (
                        isProtectedPoint.test(verts0.first())
                     && isProtectedPoint.test(verts0.second())
                     && isProtectedPoint.test(verts1.first())
                     && isProtectedPoint.test(verts1.second())
                    )
                    {
                        // Open edge on boundary face. Keep
                    }
                    else
                    {
                        // Open edge
                        openEdgeIds.insert(edgei);
                    }
                }
            }

            const label nOpenEdges
            (
                returnReduce(openEdgeIds.size(), sumOp<label>())
            );

            if (nOpenEdges)
            {
                const edgeList debugEdges
                (
                    surfEdges,
                    openEdgeIds.sortedToc()
                );

                vtk::lineWriter writer
                (
                    s.points(),
                    debugEdges,
                    fileName
                    (
                        mesh_.time().globalPath()
                      / ("isoSurfaceTopo." + timeDesc + "-openEdges")
                    )
                );

                writer.writeGeometry();

                // CellData
                writer.beginCellData();
                writer.writeProcIDs();

                Info<< "isoSurfaceTopo : (debug) wrote "
                    << nOpenEdges << " open edges: "
                    << writer.output().name() << nl;
            }
            else
            {
                Info<< "isoSurfaceTopo : no open edges" << nl;
            }
        }

        // Write debug surface with snaps
        if (this->snap())
        {
            const Mesh& s = *this;

            vtk::surfaceWriter writer
            (
                s.points(),
                s,
                fileName
                (
                    mesh_.time().globalPath()
                  / ("isoSurfaceTopo." + timeDesc + "-surface")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();
            writer.write("cellID", meshCells_);

            // PointData
            writer.beginPointData();
            {
                // NB: may have non-compact surface points
                // --> use points().size() not nPoints()!

                labelList pointStatus(s.points().size(), Zero);

                forAll(pointToVerts_, i)
                {
                    const edge& verts = pointToVerts_[i];
                    if (verts.first() == verts.second())
                    {
                        // Duplicate index (ie, snapped)
                        pointStatus[i] = 1;
                    }
                }

                writer.write("point-status", pointStatus);
            }

            Info<< "isoSurfaceTopo : (debug) wrote "
                << returnReduce(s.size(), sumOp<label>())
                << " faces : "
                << writer.output().name() << nl;
        }
    }
    tetCutAddr.clearDebug();


    if (params.filter() == filterType::NONMANIFOLD)
    {
        // We remove verts on face diagonals. This is in fact just
        // straightening the edges of the face through the cell. This can
        // close off 'pockets' of triangles and create open or
        // multiply-connected triangles

        // Solved by eroding open-edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // The list of surface faces that should be retained after erosion
        Mesh& surf = *this;
        labelList faceAddr(identity(surf.size()));

        bitSet faceSelection;

        while (true)
        {
            // Shadow the surface for the purposes of erosion
            uindirectPrimitivePatch erosion
            (
                UIndirectList<face>(surf, faceAddr),
                surf.points()
            );

            faceSelection.clear();
            faceSelection.resize(erosion.size());

            const labelList& mp = erosion.meshPoints();
            const edgeList& surfEdges = erosion.edges();
            const labelListList& edgeFaces = erosion.edgeFaces();

            label nEdgeRemove = 0;

            forAll(edgeFaces, edgei)
            {
                const labelList& eFaces = edgeFaces[edgei];
                if (eFaces.size() == 1)
                {
                    // Open edge (not originating from a boundary face)

                    const edge& e = surfEdges[edgei];
                    const edge& verts0 = pointToVerts_[mp[e.first()]];
                    const edge& verts1 = pointToVerts_[mp[e.second()]];

                    if
                    (
                        isProtectedPoint.test(verts0.first())
                     && isProtectedPoint.test(verts0.second())
                     && isProtectedPoint.test(verts1.first())
                     && isProtectedPoint.test(verts1.second())
                    )
                    {
                        // Open edge on boundary face. Keep
                    }
                    else
                    {
                        // Open edge. Mark for erosion
                        faceSelection.set(eFaces[0]);
                        ++nEdgeRemove;
                    }
                }
            }

            if (debug)
            {
                Pout<< "isoSurfaceTopo :"
                    << " removing " << faceSelection.count()
                    << " / " << faceSelection.size()
                    << " faces on " << nEdgeRemove << " open edges" << endl;
            }

            if (returnReduceAnd(faceSelection.none()))
            {
                break;
            }

            // Remove the faces from the addressing
            inplaceSubset(faceSelection, faceAddr, true);  // True = remove
        }


        // Finished erosion (if any)
        // - retain the faces listed in the updated addressing

        if (surf.size() != faceAddr.size())
        {
            faceSelection.clear();
            faceSelection.resize(surf.size());
            faceSelection.set(faceAddr);

            inplaceSubsetMesh(faceSelection);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoSurfaceTopo::inplaceSubsetMesh(const bitSet& include)
{
    labelList pointMap;
    labelList faceMap;
    Mesh filtered
    (
        Mesh::subsetMesh(include, pointMap, faceMap)
    );
    Mesh::transfer(filtered);

    meshCells_ = UIndirectList<label>(meshCells_, faceMap)();

    pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointMap)();
}


// ************************************************************************* //
