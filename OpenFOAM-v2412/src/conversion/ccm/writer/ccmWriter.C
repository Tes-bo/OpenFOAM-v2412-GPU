/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ccmWriter.H"
#include "cellModel.H"

#include "ccmInternal.H" // include last to avoid any strange interactions

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Name for the topology file reference
Foam::string Foam::ccm::writer::defaultMeshName = "meshExport";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Create a filled linear map with 'size' from 'start + 1'
void Foam::ccm::writer::addLinearMap
(
    const string& mapName,
    ccmID& mapId,
    label size,
    label start
) const
{
    if (globalState_->hasError() || size <= 0)
    {
        return;
    }

    CCMIONewEntity
    (
        &(globalState_->error),
        (globalState_->root),
        kCCMIOMap,
        mapName.c_str(),
        &mapId
    );
    assertNoError("creating linearMap '" + mapName + "' node");

    List<int> data(size);
    forAll(data, i)
    {
        data[i] = start + 1 + i;
    }

    CCMIOWriteMap
    (
        &(globalState_->error),
        mapId,
        size,
        (start + size),
        data.cdata(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("writing linearMap '" + mapName + "' data");
}


Foam::label Foam::ccm::writer::findDefaultBoundary() const
{
    return mesh_.boundaryMesh().findPatchID(defaultBoundaryName);
}


void Foam::ccm::writer::writeBoundaryRegion
(
    const ccmID& probNode
) const
{
    // Create dictionary lookup for constant/boundaryRegion
    dictionary typeDict;

    forAllConstIters(boundaryRegion_, iter)
    {
        const dictionary& dict = iter();
        word nameEntry;
        word typeEntry;

        if
        (
            dict.readIfPresent("Label", nameEntry)
         && dict.readIfPresent("BoundaryType", typeEntry)
         && !typeDict.found(nameEntry)
        )
        {
            typeDict.add(nameEntry, typeEntry);
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label defaultId = findDefaultBoundary();

    // Force write Default_Boundary_Region as BoundaryRegion-0
    {
        ccmID nodeId;
        CCMIONewIndexedEntity
        (
            &(globalState_->error),
            probNode,
            kCCMIOBoundaryRegion,
            0,
            nullptr,
            &nodeId
        );
        CCMIOWriteOptstr
        (
            nullptr,
            nodeId,
            "Label",
            defaultBoundaryName
        );
        CCMIOWriteOptstr
        (
            nullptr,
            nodeId,
            "BoundaryType",
            "wall"
        );
    }

    forAll(patches, patchI)
    {
        word patchName = patches[patchI].name();
        word patchType = patches[patchI].type();

        label regionId = patchI;
        if (regionId == defaultId)
        {
            continue;  // Skip - already written
        }
        else if (defaultId == -1 || regionId < defaultId)
        {
            ++regionId;
        }

        // Use BoundaryType from constant/boundaryRegion
        typeDict.readIfPresent(patchName, patchType);

        ccmID nodeId;
        CCMIONewIndexedEntity
        (
            &(globalState_->error),
            probNode,
            kCCMIOBoundaryRegion,
            regionId,
            nullptr,
            &nodeId
        );
        CCMIOWriteOptstr
        (
            nullptr,
            nodeId,
            "Label",
            patchName.c_str()
        );
        CCMIOWriteOptstr
        (
            nullptr,
            nodeId,
            "BoundaryType",
            patchType.c_str()
        );
    }
}


void Foam::ccm::writer::writeCellTable
(
    const ccmID& probNode
) const
{
    if (!cellTable_.size())
    {
        Info<< "No cellTable" << endl;
        return;
    }

    ccmID nodeId;

    forAllConstIters(cellTable_, iter)
    {
        label intVal = iter.key();
        const dictionary& dict = iter();

        CCMIONewIndexedEntity
        (
            &(globalState_->error),
            probNode,
            kCCMIOCellType,
            intVal,
            nullptr,
            &nodeId
        );

        wordList toc = dict.toc();
        for (const word& keyword : toc)
        {
            int pos = keyword.find("Id");

            // Tags containing 'Id' are integers
            if (pos > 0)
            {
                dict.readEntry(keyword, intVal);

                CCMIOWriteOpti
                (
                    nullptr,
                    nodeId,
                    keyword.c_str(),
                    intVal
                );
            }
            else if (pos < 0)
            {
                const word strVal(dict.get<word>(keyword));

                CCMIOWriteOptstr
                (
                    nullptr,
                    nodeId,
                    keyword.c_str(),
                    strVal.c_str()
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write out problem description
void Foam::ccm::writer::writeProblem
(
    const ccmID& stateNode
) const
{
    ccmID probNode;
    CCMIONewEntity
    (
        &(globalState_->error),
        (globalState_->root),
        kCCMIOProblemDescription,
        nullptr,
        &probNode
    );

    writeCellTable(probNode);
    writeBoundaryRegion(probNode);

    // let the state know about our problem
    CCMIOWriteState
    (
        &(globalState_->error),
        stateNode,
        probNode,
        nullptr
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ccm::writer::writer
(
    const fileName& file,
    const polyMesh& mesh,
    const bool backup
)
:
    base(),
    maps_(new ccmMaps),
    mesh_(mesh),
    // Mapping between OpenFOAM and PROSTAR primitives
    prostarShapeLookup_
    {
        { cellModel::ref(cellModel::HEX).index(), STARCDCore::starcdHex },
        { cellModel::ref(cellModel::PRISM).index(), STARCDCore::starcdPrism },
        { cellModel::ref(cellModel::TET).index(), STARCDCore::starcdTet },
        { cellModel::ref(cellModel::PYR).index(), STARCDCore::starcdPyr }
    },
    boundaryRegion_(mesh),
    cellTable_(mesh)
{
    // Writing/re-writing existing files is too annoying - start anew
    if (backup)
    {
        if (Foam::mvBak(file))
        {
            Info<< "moved existing file -> " << fileName(file + ".bak*") << nl;
        }
    }
    else if (exists(file, false))
    {
        Foam::rm(file);
        Info<< "removed existing file: " << file << nl;
    }

    // Reinitialize error state from the return value
    globalState_->error = CCMIOOpenFile
    (
        nullptr,
        file.c_str(),
        kCCMIOWrite,
        &(globalState_->root)
    );
    assertNoError("Error opening file for writing");

    // We always need the cell map
    addLinearMap
    (
        "cell Map",
        maps_->cells,
        mesh_.nCells()
    );

    // We often need the internalFaces map
    addLinearMap
    (
        "internalFaces Map",
        maps_->internalFaces,
        mesh_.nInternalFaces()
    );

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    maps_->boundary.setSize(patches.size());

    // We always need maps for the boundary regions
    forAll(patches, patchI)
    {
        if (patches[patchI].size() > 0)
        {
            string mapName = "boundaryMap-" + Foam::name(patchI);

            addLinearMap
            (
                mapName,
                maps_->boundary[patchI],
                patches[patchI].size(),
                patches[patchI].start()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ccm::writer::~writer()
{
    close();
}


// ************************************************************************* //
