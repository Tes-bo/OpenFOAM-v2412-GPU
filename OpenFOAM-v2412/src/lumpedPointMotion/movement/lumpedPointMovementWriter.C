/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "lumpedPointMovement.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "OFstream.H"
#include "foamVtkOutput.H"
#include "foamVtkSurfaceWriter.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::lumpedPointMovement::writeStateVTP
(
    const lumpedPointState& state,
    const fileName& file
) const
{
    if (!Pstream::master())
    {
        // No extra information available from slaves, write on master only.
        return;
    }

    labelListList lines;

    label nLines = controllers_.size();

    if (nLines)
    {
        lines.resize(nLines);
        nLines = 0;

        for (const word& ctrlName : controllers_.sortedToc())
        {
            lines[nLines] = controllers_[ctrlName]->pointLabels();
            ++nLines;
        }
    }
    else
    {
        // Default - global with all points as a single line
        lines.resize(1);
        lines.first() = identity(state.size());
    }

    state.writeVTP(file, lines, originalIds_);
}


void Foam::lumpedPointMovement::writeStateVTP(const fileName& file) const
{
    writeStateVTP(state(), file);
}


void Foam::lumpedPointMovement::writeForcesAndMomentsVTP
(
    const fileName& file,
    const UList<vector>& forces,
    const UList<vector>& moments
) const
{
    if (!Pstream::master())
    {
        // Force, moments already reduced
        return;
    }

    OFstream fos(file);
    std::ostream& os = fos.stdStream();

    autoPtr<vtk::formatter> format = vtk::newFormatter
    (
        os,
        vtk::formatType::INLINE_ASCII
    );

    format().xmlHeader()
        .beginVTKFile<vtk::fileTag::POLY_DATA>();

    //
    // The 'backbone' of lumped mass points
    //
    const label nPoints = state().points().size();

    {
        format()
            .tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, nPoints,
                vtk::fileAttr::NUMBER_OF_VERTS,  nPoints
            );

        // 'points'
        {
            const uint64_t payLoad = vtk::sizeofData<float, 3>(nPoints);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
            vtk::writeList(format(), state().points());
            format().flush();

            format()
                .endDataArray()
                .endTag(vtk::fileTag::POINTS);
        }

        // <Verts>
        format().tag(vtk::fileTag::VERTS);

        //
        // 'connectivity'
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPoints);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            vtk::writeIdentity(format(), nPoints);

            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = linear mapping onto points (with 1 offset)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPoints);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            vtk::writeIdentity(format(), nPoints, 1);

            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::VERTS);
        // </Verts>
    }

    format().beginPointData();

    // forces
    if (forces.size() == nPoints)
    {
        const uint64_t payLoad = vtk::sizeofData<float, 3>(nPoints);

        format().beginDataArray<float, 3>("forces");
        format().writeSize(payLoad);

        vtk::writeList(format(), forces);
        format().flush();

        format().endDataArray();
    }

    // moments
    if (moments.size() == nPoints)
    {
        const uint64_t payLoad = vtk::sizeofData<float, 3>(nPoints);

        format().beginDataArray<float, 3>("moments");
        format().writeSize(payLoad);

        vtk::writeList(format(), moments);
        format().flush();

        format().endDataArray();
    }

    format().endPointData();

    format().endPiece();

    format().endTag(vtk::fileTag::POLY_DATA)
        .endVTKFile();
}


void Foam::lumpedPointMovement::writeZonesVTP
(
    const fileName& file,
    const polyMesh& mesh,
    const pointField& points0
) const
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList patchIds(patchControls_.sortedToc());

    vtk::surfaceWriter writer
    (
        pointField::null(),
        faceList::null(),
        vtk::formatType::INLINE_ASCII,
        file
    );

    for (const label patchi : patchIds)
    {
        const labelList& faceToPoint = patchControls_[patchi].faceToPoint_;

        primitivePatch pp
        (
            SubList<face>(mesh.faces(), patches[patchi].range()),
            points0
        );

        writer.piece(pp.localPoints(), pp.localFaces());

        writer.writeGeometry();

        writer.beginCellData(2);

        writer.writeUniform("patchId", patchi);
        writer.write("lumpedId", faceToPoint);

        writer.endCellData();
    }
}


void Foam::lumpedPointMovement::writeVTP
(
    const fileName& file,
    const polyMesh& mesh,
    const pointField& points0
) const
{
    writeVTP(file, state(), mesh, points0);
}


void Foam::lumpedPointMovement::writeVTP
(
    const fileName& file,
    const lumpedPointState& state,
    const polyMesh& mesh,
    const pointField& points0
) const
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList patchIds(patchControls_.sortedToc());

    pointMesh ptMesh(mesh);

    vtk::surfaceWriter writer
    (
        pointField::null(),
        faceList::null(),
        vtk::formatType::INLINE_ASCII,
        file
    );

    for (const label patchi : patchIds)
    {
        const polyPatch& pp = patches[patchi];

        const pointPatch& ptPatch = ptMesh.boundary()[patchi];

        // Current position (not displacement)
        tmp<pointField> tpts = pointsPosition(state, ptPatch, points0);

        writer.piece(tpts(), pp.localFaces());

        writer.writeGeometry();

        // Face mapping
        const labelList& faceToPoint = patchControls_[patchi].faceToPoint_;

        writer.beginCellData(2);

        writer.writeUniform("patchId", patchi);
        writer.write("lumpedId", faceToPoint);

        writer.endCellData();

        // The interpolator
        const List<lumpedPointInterpolator>& interpList
            = patchControls_[patchi].interp_;

        writer.beginPointData(3);

        // Nearest, Next
        {
            labelList intData(interpList.size());

            forAll(interpList, i)
            {
                intData[i] = interpList[i].nearest();
            }
            writer.write("nearest", intData);

            forAll(interpList, i)
            {
                intData[i] = interpList[i].next1();
            }
            writer.write("next1", intData);


            forAll(interpList, i)
            {
                intData[i] = interpList[i].next2();
            }
            writer.write("next2", intData);
        }

        // Weights
        {
            scalarList floatData(interpList.size());

            forAll(interpList, i)
            {
                floatData[i] = interpList[i].weight0();
            }
            writer.write("weight", floatData);

            forAll(interpList, i)
            {
                floatData[i] = interpList[i].weight1();
            }
            writer.write("weight1", floatData);

            forAll(interpList, i)
            {
                floatData[i] = interpList[i].weight2();
            }
            writer.write("weight2", floatData);
        }

        writer.endPointData();
    }
}


// ************************************************************************* //
