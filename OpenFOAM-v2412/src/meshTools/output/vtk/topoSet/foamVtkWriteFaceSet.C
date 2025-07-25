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
#include "faceSet.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::vtk::writeFaceSet
(
    const polyMesh& mesh,
    const faceSet& set,
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

    // Use the faces from faceSet
    faces.addressing() = set.sortedToc();

    //-------------------------------------------------------------------------

    indirectPatchWriter writer(pp, opts);

    writer.open(file, parallel);

    writer.beginFile(set.name());
    writer.writeGeometry();


    // CellData - faceID only
    {
        writer.beginCellData(1);

        labelField faceValues(faces.addressing());

        // processor-local faceID offset
        const label faceIdOffset =
        (
            writer.parallel() ? globalIndex::calcOffset(mesh.nFaces()) : 0
        );

        if (faceIdOffset)
        {
            faceValues += faceIdOffset;
        }

        writer.write("faceID", faceValues);

        // End CellData/PointData is implicit
    }

    writer.close();

    return true;
}


// ************************************************************************* //
