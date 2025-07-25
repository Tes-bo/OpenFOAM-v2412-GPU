/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "isoSurfaceBase.H"
#include "polyMesh.H"
#include "tetMatcher.H"
#include "cyclicACMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isoSurfaceBaseMethods.H"
defineIsoSurfaceInterpolateMethods(Foam::isoSurfaceBase);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Test face for edge cuts
inline static bool isFaceCut
(
    const scalar isoval,
    const scalarField& pointValues,
    const labelUList& indices
)
{
    auto iter = indices.cbegin();
    const auto last = indices.cend();

    // if (iter == last) return false;  // ie, indices.empty()

    const bool aLower = (pointValues[*iter] < isoval);
    ++iter;

    while (iter != last)
    {
        if (aLower != (pointValues[*iter] < isoval))
        {
            return true;
        }
        ++iter;
    }

    return false;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceBase::isoSurfaceBase
(
    const polyMesh& mesh,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const isoSurfaceParams& params
)
:
    meshedSurface(),
    isoSurfaceParams(params),
    mesh_(mesh),
    cVals_(cellValues),
    pVals_(pointValues),
    iso_(iso),
    ignoreBoundaryFaces_(),
    meshCells_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoSurfaceBase::ignoreCyclics()
{
    // Determine boundary pyramids to ignore (originating from ACMI faces)
    // Maybe needs to be optional argument or more general detect colocated
    // faces.

    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (isA<cyclicACMIPolyPatch>(pp))
        {
            ignoreBoundaryFaces_.set(labelRange(pp.offset(), pp.size()));
        }
    }
}


Foam::label Foam::isoSurfaceBase::countCutType
(
    const UList<cutType>& cuts,
    const uint8_t maskValue
)
{
    label count = 0;

    for (const cutType cut : cuts)
    {
        if (maskValue ? (cut & maskValue) != 0 : !cut)
        {
            ++count;
        }
    }

    return count;
}


void Foam::isoSurfaceBase::resetCuts(UList<cutType>& cuts)
{
    for (cutType& cut : cuts)
    {
        if (cut != cutType::BLOCKED)
        {
            cut = cutType::UNVISITED;
        }
    }
}


Foam::label Foam::isoSurfaceBase::blockCells
(
    UList<cutType>& cuts,
    const bitSet& ignoreCells
) const
{
    label count = 0;

    for (const label celli : ignoreCells)
    {
        if (celli >= cuts.size())
        {
            break;
        }

        cuts[celli] = cutType::BLOCKED;
        ++count;
    }

    return count;
}


Foam::label Foam::isoSurfaceBase::blockCells
(
    UList<cutType>& cuts,
    const boundBox& bb,
    const volumeType::type volType
) const
{
    label count = 0;

    // Mark cells inside/outside bounding box as blocked
    const bool keepInside = (volType == volumeType::INSIDE);

    if (!keepInside && volType != volumeType::OUTSIDE)
    {
        // Could warn about invalid...
    }
    else if (bb.good())
    {
        const pointField& cc = mesh_.cellCentres();

        forAll(cuts, celli)
        {
            if
            (
                cuts[celli] == cutType::UNVISITED
             && (bb.contains(cc[celli]) ? keepInside : !keepInside)
            )
            {
                cuts[celli] = cutType::BLOCKED;
                ++count;
            }
        }
    }

    return count;
}


Foam::label Foam::isoSurfaceBase::calcCellCuts(List<cutType>& cuts) const
{
    // Don't consider SPHERE cuts in the total number of cells cut
    constexpr uint8_t realCut(cutType::CUT | cutType::TETCUT);

    cuts.resize(mesh_.nCells(), cutType::UNVISITED);

    label nCuts = 0;
    forAll(cuts, celli)
    {
        if (cuts[celli] == cutType::UNVISITED)
        {
            cuts[celli] = getCellCutType(celli);

            if ((cuts[celli] & realCut) != 0)
            {
                ++nCuts;
            }
        }
    }

    return nCuts;
}


Foam::isoSurfaceBase::cutType
Foam::isoSurfaceBase::getFaceCutType(const label facei) const
{
    return
    (
        (
            mesh_.isInternalFace(facei)
         || !ignoreBoundaryFaces_.test(facei-mesh_.nInternalFaces())
        )
     && isFaceCut(iso_, pVals_, mesh_.faces()[facei])
    ) ? cutType::CUT : cutType::NOTCUT;
}


Foam::isoSurfaceBase::cutType
Foam::isoSurfaceBase::getCellCutType(const label celli) const
{
    // Tet version
    if (tetMatcher::test(mesh_, celli))
    {
        for (const label facei : mesh_.cells()[celli])
        {
            if
            (
                !mesh_.isInternalFace(facei)
             && ignoreBoundaryFaces_.test(facei-mesh_.nInternalFaces())
            )
            {
                continue;
            }

            if (isFaceCut(iso_, pVals_, mesh_.faces()[facei]))
            {
                return cutType::TETCUT;
            }
        }

        return cutType::NOTCUT;
    }


    // Regular cell
    label nPyrEdges = 0;
    label nPyrCuts = 0;

    const bool cellLower = (cVals_[celli] < iso_);

    for (const label facei : mesh_.cells()[celli])
    {
        if
        (
           !mesh_.isInternalFace(facei)
         && ignoreBoundaryFaces_.test(facei-mesh_.nInternalFaces())
        )
        {
            continue;
        }

        const face& f = mesh_.faces()[facei];

        // Count pyramid edges cut
        for (const label pointi : f)
        {
            ++nPyrEdges;

            if (cellLower != (pVals_[pointi] < iso_))
            {
                ++nPyrCuts;
            }
        }
    }

    if (nPyrCuts == 0)
    {
        return cutType::NOTCUT;
    }

    // If all pyramid edges are cut (no outside faces),
    // it is a sphere cut

    return (nPyrCuts == nPyrEdges) ? cutType::SPHERE : cutType::CUT;
}


// ************************************************************************* //
