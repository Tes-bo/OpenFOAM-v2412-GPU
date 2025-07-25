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

#include "lumpedPointState.H"
#include "OFstream.H"
#include "sliceRange.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointState::writeVTP
(
    const fileName& outputFile,
    const labelListList& lines,
    const labelList& pointIds
) const
{
    if (!Pstream::master())
    {
        // No extra information available from slaves, write on master only.
        return;
    }

    // local-to-global transformation tensors
    const tensorField& localToGlobal = rotations();

    OFstream fos(outputFile);
    std::ostream& os = fos.stdStream();

    autoPtr<vtk::formatter> format = vtk::newFormatter
    (
        os,
        vtk::formatType::INLINE_ASCII
    );

    format().xmlHeader()
        .beginVTKFile<vtk::fileTag::POLY_DATA>();

    //
    // Lumped mass points and connections,
    // with triangles to visualize location/rotation
    //
    {
        const label nPoints = 3*points_.size();  // 3 points per triangle
        const label nPolys  = points_.size();

        format()
            .tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, nPoints,
                vtk::fileAttr::NUMBER_OF_VERTS,  points_.size(),
                vtk::fileAttr::NUMBER_OF_LINES,  lines.size(),
                vtk::fileAttr::NUMBER_OF_POLYS,  nPolys
            );

        // 'points'
        {
            const uint64_t payLoad = vtk::sizeofData<float, 3>(nPoints);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float, 3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);

            // The lumped points first
            vtk::writeList(format(), points_);

            // Other points (for the triangles) next
            forAll(points_, posi)
            {
                const point& origin = points_[posi];
                const tensor& rotTensor =
                (
                    posi < localToGlobal.size()
                  ? localToGlobal[posi]
                  : pTraits<tensor>::I
                );

                // Local-to-global rotation and translation
                vtk::write(format(), 2*visLength*rotTensor.cx() + origin);
                vtk::write(format(), 1*visLength*rotTensor.cy() + origin);
            }

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
            const uint64_t payLoad = vtk::sizeofData<label>(points_.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            vtk::writeIdentity(format(), points_.size());

            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = linear mapping onto points (with 1 offset)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(points_.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            vtk::writeIdentity(format(), points_.size(), 1);

            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::VERTS);
        // </Verts>


        // <Lines>
        format().tag(vtk::fileTag::LINES);

        label nLinePoints = 0;
        for (const labelList& linePoints : lines)
        {
            nLinePoints += linePoints.size();
        }

        //
        // 'connectivity'
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nLinePoints);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            for (const labelList& linePoints : lines)
            {
                vtk::writeList(format(), linePoints);
            }

            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = N lines
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(lines.size());

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            nLinePoints = 0;
            for (const labelList& linePoints : lines)
            {
                nLinePoints += linePoints.size();
                format().write(nLinePoints);
            }

            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::LINES);
        // </Lines>


        // <Polys>
        format().tag(vtk::fileTag::POLYS);

        //
        // 'connectivity' - 3 points (ie, tri)
        // origins appear first, followed by a point pair for each triangle
        // Eg,
        // - tri 0: (0 N   N+1)
        // - tri 1: (1 N+2 N+3)
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(3*nPolys);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);

            for (label pointi=0, nei=nPolys; pointi < nPolys; ++pointi)
            {
                format().write(pointi);
                format().write(nei); ++nei;
                format().write(nei); ++nei;
            }

            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = single tri
        //
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPolys);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);

            for (const label off : sliceRange(3, nPolys, 3))
            {
                format().write(off);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);
        // </Polys>


        // CELL_DATA
        format().beginCellData();

        // point id
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(points_.size() + lines.size());

            format().beginDataArray<label>("pointId");
            format().writeSize(payLoad);

            // <Verts>
            vtk::writeIdentity(format(), points_.size());

            // <Lines>
            vtk::write(format(), label(-1), lines.size());

            // <Poly>
            vtk::writeIdentity(format(), nPolys);

            format().flush();

            format().endDataArray();
        }

        // original id
        if (pointIds.size() == points_.size())
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(points_.size() + lines.size());

            format().beginDataArray<label>("originalId");
            format().writeSize(payLoad);

            // <Verts>
            vtk::writeList(format(), pointIds);

            // <Lines>
            vtk::write(format(), label(-1), lines.size());

            // <Poly>
            vtk::writeList(format(), pointIds);

            format().flush();

            format().endDataArray();
        }

        // line id
        {
            const uint64_t payLoad =
                vtk::sizeofData<label>(points_.size() + lines.size());

            format().beginDataArray<label>("lineId");
            format().writeSize(payLoad);

            // <Verts>
            vtk::write(format(), label(-1), points_.size());

            // <Lines>
            vtk::writeIdentity(format(), lines.size());

            // <Poly>
            vtk::write(format(), label(-1), nPolys);

            format().flush();

            format().endDataArray();
        }

        format().endCellData();


        // POINT_DATA
        format().beginPointData();

        // point id
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPoints);

            format().beginDataArray<label>("pointId");
            format().writeSize(payLoad);

            // The lumped points first
            vtk::writeIdentity(format(), points_.size());

            // Tag other (triangle) points as -1
            vtk::write(format(), label(-1), 2*points_.size());

            format().flush();

            format().endDataArray();
        }

        // original id
        if (pointIds.size() == points_.size())
        {
            const uint64_t payLoad = vtk::sizeofData<label>(nPoints);

            format().beginDataArray<label>("originalId");
            format().writeSize(payLoad);

            // The lumped points first
            vtk::writeList(format(), pointIds);

            // Tag other (triangle) points as -1
            vtk::write(format(), label(-1), 2*points_.size());

            format().flush();

            format().endDataArray();
        }

        format().endPointData();

        format().endPiece();
    }

    format().endTag(vtk::fileTag::POLY_DATA)
        .endVTKFile();
}


// ************************************************************************* //
