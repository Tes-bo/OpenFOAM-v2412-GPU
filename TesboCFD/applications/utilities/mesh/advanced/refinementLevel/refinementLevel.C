/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    refinementLevel

Group
    grpMeshAdvancedUtilities

Description
    Attempt to determine the refinement levels of a refined cartesian mesh.
    Run BEFORE snapping.

    Writes
    - volScalarField 'refinementLevel' with current refinement level.
    - cellSet 'refCells' which are the cells that need to be refined to satisfy
      2:1 refinement.

    Works by dividing cells into volume bins.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "SortList.H"
#include "labelIOList.H"
#include "fvMesh.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return true if any cells had to be split to keep a difference between
// neighbouring refinement levels < limitDiff. Puts cells into refCells and
// update refLevel to account for refinement.
bool limitRefinementLevel
(
    const primitiveMesh& mesh,
    labelList& refLevel,
    cellSet& refCells
)
{
    const labelListList& cellCells = mesh.cellCells();

    label oldNCells = refCells.size();

    forAll(cellCells, celli)
    {
        const labelList& cCells = cellCells[celli];

        forAll(cCells, i)
        {
            if (refLevel[cCells[i]] > (refLevel[celli]+1))
            {
                // Found neighbour with >=2 difference in refLevel.
                refCells.insert(celli);
                refLevel[celli]++;
                break;
            }
        }
    }

    if (refCells.size() > oldNCells)
    {
        Info<< "Added an additional " << refCells.size() - oldNCells
            << " cells to satisfy 1:2 refinement level"
            << endl;

        return true;
    }

    return false;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Attempt to determine refinement levels of a refined cartesian mesh.\n"
        "Run BEFORE snapping!"
    );

    argList::addBoolOption
    (
        "readLevel",
        "Read level from refinementLevel file"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    Info<< "Dividing cells into bins depending on cell volume.\nThis will"
        << " correspond to refinement levels for a mesh with only 2x2x2"
        << " refinement\n"
        << "The upper range for every bin is always 1.1 times the lower range"
        << " to allow for some truncation error."
        << nl << endl;

    const bool readLevel = args.found("readLevel");

    const scalarField& vols = mesh.cellVolumes();

    SortList<scalar> sortedVols(vols);

    // All cell labels, sorted per bin.
    DynamicList<DynamicList<label>> bins;

    // Lower/upper limits
    DynamicList<scalarMinMax> limits;

    // Create bin0. Have upperlimit as factor times lowerlimit.
    bins.emplace_back();
    limits.emplace_back(sortedVols[0], 1.1*sortedVols[0]);

    forAll(sortedVols, i)
    {
        if (sortedVols[i] > limits.back().max())
        {
            // New value outside of current bin
            Info<< "Collected " << bins.back() << " elements in bin "
                << limits.back().min() << " .. "
                << limits.back().max() << endl;

            // Create new bin.
            bins.emplace_back();
            limits.emplace_back(sortedVols[i], 1.1 * sortedVols[i]);

            Info<< "Creating new bin "
                << limits.back().min() << " .. "
                << limits.back().max() << endl;
        }

        // Add to current bin.
        bins.back().push_back(sortedVols.indices()[i]);
    }
    Info<< endl;

    //
    // Write to cellSets.
    //

    Info<< "Volume bins:" << nl;
    forAll(bins, bini)
    {
        const auto& bin = bins[bini];

        cellSet cells(mesh, "vol" + Foam::name(bini), bin.size());
        cells.insert(bin);

        Info<< "    " << limits[bini].min() << " .. " << limits[bini].max()
            << "  : writing " << bin.size() << " cells to cellSet "
            << cells.name() << endl;

        cells.write();
    }



    //
    // Convert bins into refinement level.
    //


    // Construct fvMesh to be able to construct volScalarField

    fvMesh fMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointField(mesh.points()),  // Could we safely re-use the data?
        faceList(mesh.faces()),
        cellList(mesh.cells())
    );

    // Add the boundary patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    polyPatchList newPatches(patches.size());

    forAll(newPatches, patchi)
    {
        newPatches.set
        (
            patchi,
            patches[patchi].clone(fMesh.boundaryMesh())
        );
    }

    fMesh.addFvPatches(newPatches);


    // Refinement level
    IOobject refHeader
    (
        "refinementLevel",
        runTime.timeName(),
        polyMesh::defaultRegion,
        runTime
    );

    if (!readLevel && refHeader.typeHeaderOk<labelIOList>(true))
    {
        WarningInFunction
            << "Detected " << refHeader.name() << " file in "
            << polyMesh::defaultRegion <<  " directory. Please remove to"
            << " recreate it or use the -readLevel option to use it"
            << endl;
        return 1;
    }


    labelIOList refLevel
    (
        IOobject
        (
            "refinementLevel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nCells(), Zero)
    );

    if (readLevel)
    {
        refLevel = labelIOList(refHeader);
    }

    // Construct volScalarField with same info for post processing
    volScalarField postRefLevel
    (
        IOobject
        (
            "refinementLevel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fMesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    // Set cell values
    forAll(bins, bini)
    {
        const auto& bin = bins[bini];

        forAll(bin, i)
        {
            refLevel[bin[i]] = bins.size() - bini - 1;
            postRefLevel[bin[i]] = refLevel[bin[i]];
        }
    }

    volScalarField::Boundary& postRefLevelBf =
        postRefLevel.boundaryFieldRef();

    // For volScalarField: set boundary values to same as cell.
    // Note: could also put
    // zeroGradient b.c. on postRefLevel and do evaluate.
    forAll(postRefLevel.boundaryField(), patchi)
    {
        const polyPatch& pp = patches[patchi];

        fvPatchScalarField& bField = postRefLevelBf[patchi];

        Info<< "Setting field for patch "<< endl;

        forAll(bField, facei)
        {
            label own = mesh.faceOwner()[pp.start() + facei];

            bField[facei] = postRefLevel[own];
        }
    }

    Info<< "Determined current refinement level and writing to "
        << postRefLevel.name() << " (as volScalarField; for post processing)"
        << nl
        << polyMesh::defaultRegion/refLevel.name()
        << " (as labelIOList; for meshing)" << nl
        << endl;

    refLevel.write();
    postRefLevel.write();


    // Find out cells to refine to keep to 2:1 refinement level restriction

    // Cells to refine
    cellSet refCells(mesh, "refCells", 100);

    while
    (
        limitRefinementLevel
        (
            mesh,
            refLevel,       // current refinement level
            refCells        // cells to refine
        )
    )
    {}

    if (refCells.size())
    {
        Info<< "Collected " << refCells.size() << " cells that need to be"
            << " refined to get closer to overall 2:1 refinement level limit"
            << nl
            << "Written cells to be refined to cellSet " << refCells.name()
            << nl << endl;

        refCells.write();

        Info<< "After refinement this tool can be run again to see if the 2:1"
            << " limit is observed all over the mesh" << nl << endl;
    }
    else
    {
        Info<< "All cells in the mesh observe the 2:1 refinement level limit"
            << nl << endl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
