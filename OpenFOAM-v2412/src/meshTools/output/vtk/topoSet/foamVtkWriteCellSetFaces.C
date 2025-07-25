/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "foamVtkWriteTopoSet.H"
#include "foamVtkIndPatchWriter.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::vtk::writeCellSetFaces
(
    const polyMesh& mesh,
    const cellSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
{
    typedef IndirectList<face> FaceListType;

    indirectPrimitivePatch pp
    (
        FaceListType(mesh.faces(), labelList()),
        mesh.points()
    );
    FaceListType& faces = pp;


    //-------------------------------------------------------------------------

    // External faces of cellset with OpenFOAM cellID as value

    Map<label> cellFaces(2*set.size());

    for (const label celli : static_cast<const labelHashSet&>(set))
    {
        const cell& cFaces = mesh.cells()[celli];

        for (const label facei : cFaces)
        {
            if (mesh.isInternalFace(facei))
            {
                label otherCelli = mesh.faceOwner()[facei];

                if (otherCelli == celli)
                {
                    otherCelli = mesh.faceNeighbour()[facei];
                }

                if (!set.found(otherCelli))
                {
                    cellFaces.insert(facei, celli);
                }
            }
            else
            {
                cellFaces.insert(facei, celli);
            }
        }
    }

    // Use these faces
    faces.addressing() = cellFaces.sortedToc();

    //-------------------------------------------------------------------------

    indirectPatchWriter writer(pp, opts);

    writer.open(file, parallel);

    writer.beginFile(set.name());
    writer.writeGeometry();

    //-------------------------------------------------------------------------

    // CellData - faceID only

    writer.beginCellData(1);
    {
        // For each face, the corresponding cellID

        labelList faceValues(faces.size());

        const labelList& faceIds = faces.addressing();

        // processor-local cellID offset
        const label cellIdOffset =
        (
            writer.parallel() ? globalIndex::calcOffset(mesh.nCells()) : 0
        );

        forAll(faceValues, facei)
        {
            faceValues[facei] = cellFaces[faceIds[facei]] + cellIdOffset;
        }

        writer.write("faceID", faceValues);
    }

    // End CellData/PointData is implicit

    writer.close();

    return true;
}


// ************************************************************************* //
