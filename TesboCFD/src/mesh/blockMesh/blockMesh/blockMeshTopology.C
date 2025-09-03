/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "blockMesh.H"
#include "blockMeshTools.H"
#include "Time.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Replace (<block> <face>) face description
// with the corresponding block face
// - otherwise check point labels are in range

template<class Source>
static void rewritePatchLabels
(
    const Source& source,
    const word& patchName,
    const PtrList<block>& blocks,
    const label nPoints,
    faceList& patchFaces
)
{
    const label nBlocks = blocks.size();

    forAll(patchFaces, facei)
    {
        face& f = patchFaces[facei];

        // Replace (<block> <face>) face description
        // with the corresponding block face
        if (f.size() == 2)
        {
            const label bi = f[0];
            const label fi = f[1];

            if (bi < 0 || bi >= nBlocks)
            {
                FatalIOErrorInFunction(source)
                    << "Block index out of range."
                    << " Patch face (" << bi << ' ' << fi << ")\n"
                    << "    Number of blocks = " << nBlocks
                    << ", block index = " << bi << nl
                    << "    on patch " << patchName << ", face " << facei
                    << exit(FatalIOError);
            }
            else if (fi >= blocks[bi].blockShape().nFaces())
            {
                FatalIOErrorInFunction(source)
                    << "Block face index out of range."
                    << " Patch face (" << bi << ' ' << fi << ")\n"
                    << "    Number of block faces = "
                    << blocks[bi].blockShape().nFaces()
                    << ", face index = " << fi << nl
                    << "    on patch " << patchName << ", face " << facei
                    << exit(FatalIOError);
            }
            else
            {
                f = blocks[bi].blockShape().face(fi);
            }
        }
        else
        {
            for (const label pointi : f)
            {
                if (pointi < 0 || pointi >= nPoints)
                {
                    FatalIOErrorInFunction(source)
                        << "Point label " << pointi
                        << " out of range 0.." << (nPoints - 1) << nl
                        << "    on patch " << patchName << ", face " << facei
                        << exit(FatalIOError);
                }
            }
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMesh::readPatches
(
    const dictionary& meshDescription,
    faceListList& tmpBlocksPatches,
    wordList& patchNames,
    wordList& patchTypes,
    wordList& nbrPatchNames
)
{
    // Collect all variables
    dictionary varDict(meshDescription.subOrEmptyDict("namedVertices"));
    varDict.merge(meshDescription.subOrEmptyDict("namedBlocks"));


    ITstream& patchStream = meshDescription.lookup("patches");

    // Read number of patches in mesh
    label nPatches = 0;

    if (patchStream.peek().isLabel())
    {
        patchStream >> nPatches;

        tmpBlocksPatches.setSize(nPatches);
        patchNames.setSize(nPatches);
        patchTypes.setSize(nPatches);
        nbrPatchNames.setSize(nPatches);
    }

    // Read beginning of blocks
    patchStream.readBegin("patches");

    nPatches = 0;

    token lastToken(patchStream);
    while (!lastToken.isPunctuation(token::END_LIST))
    {
        if (tmpBlocksPatches.size() <= nPatches)
        {
            tmpBlocksPatches.setSize(nPatches + 1);
            patchNames.setSize(nPatches + 1);
            patchTypes.setSize(nPatches + 1);
            nbrPatchNames.setSize(nPatches + 1);
        }

        patchStream.putBack(lastToken);

        patchStream
            >> patchTypes[nPatches]
            >> patchNames[nPatches];

        // Read patch faces
        tmpBlocksPatches[nPatches] = blockMeshTools::read<face>
        (
            patchStream,
            varDict
        );


        // Check for multiple patches
        for (label i = 0; i < nPatches; i++)
        {
            if (patchNames[nPatches] == patchNames[i])
            {
                FatalIOErrorInFunction(patchStream)
                    << "Duplicate patch " << patchNames[nPatches]
                    << " at line " << patchStream.lineNumber()
                    << exit(FatalIOError);
            }
        }

        rewritePatchLabels
        (
            patchStream,
            patchNames[nPatches],
            *this,
            vertices_.size(),
            tmpBlocksPatches[nPatches]
        );

        ++nPatches;


        // Split old style cyclics

        if (patchTypes[nPatches-1] == cyclicPolyPatch::typeName)
        {
            word halfA = patchNames[nPatches-1] + "_half0";
            word halfB = patchNames[nPatches-1] + "_half1";

            FatalIOErrorInFunction(patchStream)
                << "Old-style cyclic definition."
                << " Splitting patch "
                << patchNames[nPatches-1] << " into two halves "
                << halfA << " and " << halfB << endl
                << "    Alternatively use new 'boundary' dictionary syntax."
                << exit(FatalIOError);

            // Add extra patch
            if (tmpBlocksPatches.size() <= nPatches)
            {
                tmpBlocksPatches.setSize(nPatches + 1);
                patchNames.setSize(nPatches + 1);
                patchTypes.setSize(nPatches + 1);
                nbrPatchNames.setSize(nPatches + 1);
            }

            // Update halfA info
            patchNames[nPatches-1] = halfA;
            nbrPatchNames[nPatches-1] = halfB;

            // Update halfB info
            patchTypes[nPatches] = patchTypes[nPatches-1];
            patchNames[nPatches] = halfB;
            nbrPatchNames[nPatches] = halfA;

            // Split faces
            if ((tmpBlocksPatches[nPatches-1].size() % 2) != 0)
            {
                FatalIOErrorInFunction(patchStream)
                    << "Size of cyclic faces is not a multiple of 2 :"
                    << tmpBlocksPatches[nPatches-1]
                    << exit(FatalIOError);
            }
            label sz = tmpBlocksPatches[nPatches-1].size()/2;
            faceList unsplitFaces(tmpBlocksPatches[nPatches-1], true);
            tmpBlocksPatches[nPatches-1] = faceList
            (
                SubList<face>(unsplitFaces, sz)
            );
            tmpBlocksPatches[nPatches] = faceList
            (
                SubList<face>(unsplitFaces, sz, sz)
            );

            nPatches++;
        }

        patchStream >> lastToken;
    }
    patchStream.putBack(lastToken);

    // Read end of blocks
    patchStream.readEnd("patches");
}


void Foam::blockMesh::readBoundary
(
    const dictionary& meshDescription,
    wordList& patchNames,
    faceListList& tmpBlocksPatches,
    PtrList<dictionary>& patchDicts
)
{
    // Collect all variables
    dictionary varDict(meshDescription.subOrEmptyDict("namedVertices"));
    varDict.merge(meshDescription.subOrEmptyDict("namedBlocks"));


    // Read like boundary file
    const PtrList<entry> patchesInfo
    (
        meshDescription.lookup("boundary")
    );

    patchNames.setSize(patchesInfo.size());
    tmpBlocksPatches.setSize(patchesInfo.size());
    patchDicts.setSize(patchesInfo.size());

    forAll(tmpBlocksPatches, patchi)
    {
        const entry& patchInfo = patchesInfo[patchi];

        if (!patchInfo.isDict())
        {
            FatalIOErrorInFunction(meshDescription)
                << "Entry " << patchInfo << " in boundary section is not a"
                << " valid dictionary."
                << exit(FatalIOError);
        }

        patchNames[patchi] = patchInfo.keyword();

        // Construct patch dictionary
        patchDicts.set(patchi, new dictionary(patchInfo.dict()));

        // Read block faces
        tmpBlocksPatches[patchi] = blockMeshTools::read<face>
        (
            patchDicts[patchi].lookup("faces"),
            varDict
        );


        rewritePatchLabels
        (
            patchInfo.dict(),
            patchNames[patchi],
            *this,
            vertices_.size(),
            tmpBlocksPatches[patchi]
        );
    }
}


Foam::cellShapeList Foam::blockMesh::getBlockShapes() const
{
    const blockMesh& blocks = *this;

    cellShapeList shapes(blocks.size());

    forAll(blocks, blocki)
    {
        shapes[blocki] = blocks[blocki].blockShape();
    }

    return shapes;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh>
Foam::blockMesh::createTopology
(
    const IOdictionary& meshDescription,
    const word& regionName
)
{
    word defaultPatchName = "defaultFaces";
    word defaultPatchType = emptyPolyPatch::typeName;

    // Read the names/types for the unassigned patch faces
    // this is a bit heavy handed (and ugly), but there is currently
    // no easy way to rename polyMesh patches subsequently
    if (const dictionary* dictptr = meshDescription.findDict("defaultPatch"))
    {
        dictptr->readIfPresent("name", defaultPatchName);
        dictptr->readIfPresent("type", defaultPatchType);
    }


    // Scaling, transformations
    readPointTransforms(meshDescription);

    // Read the block edges
    if (meshDescription.found("edges"))
    {
        if (verbose_)
        {
            Info<< "Creating block edges" << endl;
        }

        blockEdgeList edges
        (
            meshDescription.lookup("edges"),
            blockEdge::iNew(meshDescription, geometry_, vertices_)
        );

        edges_.transfer(edges);
    }
    else
    {
        edges_.clear();
        if (verbose_)
        {
            Info<< "No non-linear block edges defined" << endl;
        }
    }


    // Read the block faces
    if (meshDescription.found("faces"))
    {
        if (verbose_)
        {
            Info<< "Creating block faces" << endl;
        }

        blockFaceList faces
        (
            meshDescription.lookup("faces"),
            blockFace::iNew(meshDescription, geometry_)
        );

        faces_.transfer(faces);
    }
    else
    {
        faces_.clear();
        if (verbose_)
        {
            Info<< "No non-planar block faces defined" << endl;
        }
    }


    // Create the blocks
    if (verbose_)
    {
        Info<< "Creating topology blocks" << endl;
    }
    {
        blockList blocks
        (
            meshDescription.lookup("blocks"),
            block::iNew(meshDescription, vertices_, edges_, faces_)
        );

        transfer(blocks);
    }


    // Create patch information

    faceListList tmpBlocksPatches;
    wordList patchNames;
    PtrList<dictionary> patchDicts;


    if (verbose_)
    {
        Info<< nl << "Creating topology patches - ";
    }

    if (meshDescription.found("patches"))
    {
        if (verbose_)
        {
            Info<< "from patches section" << endl;
        }

        wordList patchTypes;
        wordList nbrPatchNames;

        readPatches
        (
            meshDescription,
            tmpBlocksPatches,
            patchNames,
            patchTypes,
            nbrPatchNames
        );

        if (verbose_)
        {
            Info<< nl
                << "Reading physicalType from existing boundary file" << endl;
        }

        patchDicts.resize(patchNames.size());

        preservePatchTypes
        (
            meshDescription.time(),
            meshDescription.time().constant(),
            polyMesh::meshSubDir,
            patchNames,
            patchDicts,
            defaultPatchName,
            defaultPatchType
        );


        // Add cyclic info (might not be present from older file)
        forAll(patchDicts, patchi)
        {
            if (!patchDicts.set(patchi))
            {
                patchDicts.set(patchi, new dictionary());
            }

            dictionary& dict = patchDicts[patchi];

            // Add but not override type
            if (!dict.found("type"))
            {
                dict.add("type", patchTypes[patchi], false);
            }
            else if (dict.get<word>("type") != patchTypes[patchi])
            {
                FatalIOErrorInFunction(meshDescription)
                    << "For patch " << patchNames[patchi]
                    << " overriding type '" << patchTypes[patchi]
                    << "' with '" << dict.get<word>("type")
                    << "' (read from boundary file)"
                    << exit(FatalIOError);
            }

            // Override neighbour patch name
            if (!nbrPatchNames[patchi].empty())
            {
                dict.set("neighbourPatch", nbrPatchNames[patchi]);
            }
        }
    }
    else if (meshDescription.found("boundary"))
    {
        if (verbose_)
        {
            Info<< "from boundary section" << endl;
        }

        readBoundary
        (
            meshDescription,
            patchNames,
            tmpBlocksPatches,
            patchDicts
        );
    }
    else
    {
        if (verbose_)
        {
            Info<< "with default boundary only!!" << endl;
        }
    }


    if (verbose_)
    {
        Info<< nl << "Creating block mesh topology";
        if (hasPointTransforms())
        {
            Info<< " - scaling/transform applied later";
        }
        Info<< endl;
    }

    auto blockMeshPtr = autoPtr<polyMesh>::New
    (
        IOobject
        (
            regionName,
            meshDescription.time().constant(),
            meshDescription.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointField(vertices_),   // Use a copy of vertices
        getBlockShapes(),
        tmpBlocksPatches,
        patchNames,
        patchDicts,
        defaultPatchName,
        defaultPatchType
    );

    check(*blockMeshPtr, meshDescription);

    return blockMeshPtr;
}


// ************************************************************************* //
