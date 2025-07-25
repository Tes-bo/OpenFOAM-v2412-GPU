/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "ensightOutput.H"
#include "cell.H"
#include "cellShape.H"
#include "face.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "manifoldCellsMeshObject.H"
#include "debug.H"
#include "defineDebugSwitch.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * * * Globals * * * * * * * * * * * * * * * * * //

defineDebugSwitchWithName(Foam::ensightOutput, "ensight", 0);

registerDebugSwitchWithName(Foam::ensightOutput, ensightOutput, "ensight");

int Foam::ensightOutput::maxChunk_
(
    Foam::debug::optimisationSwitch("ensight.maxChunk", 0)
);

registerOptSwitch("ensight.maxChunk", int, Foam::ensightOutput::maxChunk_);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Sizes

Foam::labelList Foam::ensightOutput::Detail::getFaceSizes
(
    const UList<face>& faces
)
{
    labelList list(faces.size());

    auto outIter = list.begin();

    for (const face& f : faces)
    {
        *outIter = f.size();
        ++outIter;
    }

    return list;
}


Foam::labelList Foam::ensightOutput::Detail::getFaceSizes
(
    const UIndirectList<face>& faces
)
{
    labelList list(faces.size());

    auto outIter = list.begin();

    for (const face& f : faces)
    {
        *outIter = f.size();
        ++outIter;
    }

    return list;
}


Foam::labelList Foam::ensightOutput::Detail::getPolysNFaces
(
    const polyMesh& mesh,
    const labelUList& addr
)
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();

    labelList list(addr.size());

    auto outIter = list.begin();

    // The number of faces per element
    for (const label cellId : addr)
    {
        *outIter = meshCells[cellId].size();
        ++outIter;
    }

    return list;
}


Foam::labelList Foam::ensightOutput::Detail::getPolysNPointsPerFace
(
    const polyMesh& mesh,
    const labelUList& addr
)
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    // Count the number of faces per element

    label nTotFaces = 0;
    for (const label cellId : addr)
    {
        nTotFaces += meshCells[cellId].size();
    }

    labelList list(nTotFaces);

    auto outIter = list.begin();

    // The number of points per element face
    for (const label cellId : addr)
    {
        for (const label facei : meshCells[cellId])
        {
            *outIter = meshFaces[facei].size();
            ++outIter;
        }
    }

    return list;
}


void Foam::ensightOutput::Detail::writeLabelListList
(
    ensightGeoFile& os,
    const labelUList& offsets,
    const labelUList& values,
    const label pointOffset
)
{
    const label off = (pointOffset + 1);  // 1-based for Ensight

    const label nLists = (offsets.size() - 1);

    for (label i = 0; i < nLists; ++i)
    {
        const labelUList list
        (
            values.slice(offsets[i], (offsets[i+i] - offsets[i]))
        );
        for (const label pointi : list)
        {
            os.write(pointi + off);
        }
        os.newline();  // One list (cell/faces) per line (ASCII)
    }
}


void Foam::ensightOutput::writeFaceList
(
    ensightGeoFile& os,
    const UList<face>& faces,
    const label pointOffset
)
{
    ensightOutput::Detail::writeLabelListList(os, faces, pointOffset);
}


void Foam::ensightOutput::writeFaceList
(
    ensightGeoFile& os,
    const UIndirectList<face>& faces,
    const label pointOffset
)
{
    ensightOutput::Detail::writeLabelListList(os, faces, pointOffset);
}


void Foam::ensightOutput::writeFaceList
(
    ensightGeoFile& os,
    const CompactListList<label>& faces,
    const label pointOffset
)
{
    ensightOutput::Detail::writeLabelListList(os, faces, pointOffset);
}


void Foam::ensightOutput::writeCellShapes
(
    ensightGeoFile& os,
    const UList<cellShape>& shapes,
    const label pointOffset
)
{
    ensightOutput::Detail::writeLabelListList(os, shapes, pointOffset);
}


Foam::CompactListList<Foam::label>
Foam::ensightOutput::Detail::getPolysFacePoints
(
    const polyMesh& mesh,
    const labelUList& addr,
    const labelList& pointMap
)
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();
    const labelList& owner = mesh.faceOwner();


    // The caller should have already checked for possible overflow,
    // so can skip that here.
    // but still need the sizing for allocations

    label nFaces = 0, nPoints = 0;
    for (const label cellId : addr)
    {
        nFaces += meshCells[cellId].size();

        for (const label faceId : meshCells[cellId])
        {
            nPoints += meshFaces[faceId].size();
        }
    }


    CompactListList<label> compact(nFaces, nPoints);
    labelList& offsets = compact.offsets();
    labelList& verts = compact.values();

    // Restart counts
    nFaces = nPoints = 0;

    for (const label cellId : addr)
    {
        for (const label faceId : meshCells[cellId])
        {
            const face& f = meshFaces[faceId];

            offsets[nFaces++] = nPoints;

            if (faceId < owner.size() && owner[faceId] != cellId)
            {
                // The neighbour of an internal face
                // - handle like face::reverseFace()

                verts[nPoints++] = pointMap[f[0]];
                for (label pti = f.size()-1; pti > 0; --pti)
                {
                    verts[nPoints++] = pointMap[f[pti]];
                }
            }
            else
            {
                for (const label pointi : f)
                {
                    verts[nPoints++] = pointMap[pointi];
                }
            }
        }
    }

    // Finally
    if (nFaces)
    {
        offsets[nFaces] = nPoints;
    }

    return compact;
}


void Foam::ensightOutput::writePolysPoints
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    const labelUList& addr,
    const labelList& pointMap
)
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();
    const labelList& owner = mesh.faceOwner();

    const label off = (1);  // 1-based for Ensight

    for (const label cellId : addr)
    {
        for (const label faceId : meshCells[cellId])
        {
            const face& f = meshFaces[faceId];

            if (faceId < owner.size() && owner[faceId] != cellId)
            {
                // The neighbour of an internal face
                // - write as face::reverseFace()

                os.write(pointMap[f[0]] + off);
                for (label pti = f.size()-1; pti > 0; --pti)
                {
                    os.write(pointMap[f[pti]] + off);
                }
            }
            else
            {
                for (const label pointi : f)
                {
                    os.write(pointMap[pointi] + off);
                }
            }

            os.newline();  // One face per line (ASCII)
        }
    }
}


void Foam::ensightOutput::writePolysPoints
(
    ensightGeoFile& os,
    const cellUList& meshCells,
    const labelUList& addr,
    const faceUList& meshFaces,
    const labelUList& owner
)
{
    const label off = (1);  // 1-based for Ensight

    for (const label cellId : addr)
    {
        for (const label faceId : meshCells[cellId])
        {
            const face& f = meshFaces[faceId];

            if (faceId < owner.size() && owner[faceId] != cellId)
            {
                // The neighbour of an internal face
                // - write as face::reverseFace()

                os.write(f[0] + off);
                for (label pti = f.size()-1; pti > 0; --pti)
                {
                    os.write(f[pti] + off);
                }
            }
            else
            {
                for (const label pointi : f)
                {
                    os.write(pointi + off);
                }
            }

            os.newline();  // One face per line (ASCII)
        }
    }
}


void Foam::ensightOutput::writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const faceUList& faces,
    bool parallel
)
{
    if (!nTotal)
    {
        return;
    }

    parallel = parallel && Pstream::parRun();

    const IntRange<int> senders =
    (
        parallel
      ? Pstream::subProcs()
      : IntRange<int>()
    );

    if (Pstream::master())
    {
        os.writeKeyword(ensightFaces::key(etype));
        os.write(nTotal);
        os.newline();
    }

    if (etype == ensightFaces::elemType::NSIDED)
    {
        // Face sizes (number of points per face)

        labelList send(ensightOutput::Detail::getFaceSizes(faces));

        if (Pstream::master())
        {
            // Main
            os.writeLabels(send);

            // Others
            for (const int proci : senders)
            {
                IPstream fromOther(Pstream::commsTypes::scheduled, proci);
                labelList recv(fromOther);

                os.writeLabels(recv);
            }
        }
        else if (senders)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }


    // List of points id for each face
    if (Pstream::master())
    {
        // Main
        writeFaceList(os, faces);

        // Others
        for (const int proci : senders)
        {
            IPstream fromOther(Pstream::commsTypes::scheduled, proci);
            List<face> recv(fromOther);

            writeFaceList(os, recv);
        }
    }
    else if (senders)
    {
        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster << faces;
    }
}


void Foam::ensightOutput::writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const UIndirectList<face>& faces,
    bool parallel
)
{
    if (!nTotal)
    {
        return;
    }

    parallel = parallel && Pstream::parRun();

    const IntRange<int> senders =
    (
        parallel
      ? Pstream::subProcs()
      : IntRange<int>()
    );


    if (Pstream::master())
    {
        os.writeKeyword(ensightFaces::key(etype));
        os.write(nTotal);
        os.newline();
    }

    if (etype == ensightFaces::elemType::NSIDED)
    {
        // Face sizes (number of points per face)

        labelList send(ensightOutput::Detail::getFaceSizes(faces));

        if (Pstream::master())
        {
            // Main
            os.writeLabels(send);

            // Others
            for (const int proci : senders)
            {
                IPstream fromOther(Pstream::commsTypes::scheduled, proci);
                labelList recv(fromOther);

                os.writeLabels(recv);
            }
        }
        else if (senders)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }


    // List of points id per face

    if (Pstream::master())
    {
        // Main
        writeFaceList(os, faces);

        // Others
        for (const int proci : senders)
        {
            IPstream fromOther(Pstream::commsTypes::scheduled, proci);
            List<face> recv(fromOther);

            writeFaceList(os, recv);
        }
    }
    else if (senders)
    {
        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster << faces;
    }
}


void Foam::ensightOutput::writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    bool parallel
)
{
    for (label typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        writeFaceConnectivity
        (
            os,
            etype,
            part.total(etype),
            UIndirectList<face>(faces, part.faceIds(etype)),
            parallel
        );
    }
}


void Foam::ensightOutput::writeFaceConnectivityPresorted
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    bool parallel
)
{
    for (label typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        writeFaceConnectivity
        (
            os,
            etype,
            part.total(etype),
            SubList<face>(faces, part.range(etype)),
            parallel
        );
    }
}


// ************************************************************************* //
