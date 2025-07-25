/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "foamVtkLagrangianWriter.H"
#include "Cloud.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::vtk::lagrangianWriter::cloudDir() const
{
    return (cloud::prefix/cloudName_);
}


Foam::pointField Foam::vtk::lagrangianWriter::positions() const
{
    Cloud<passiveParticle> parcels(mesh_, cloudName_, false);

    pointField pts(parcels.size());

    auto outIter = pts.begin();

    for (const auto& p : parcels)
    {
        *outIter = p.position();
        ++outIter;
    }

    return pts;
}


void Foam::vtk::lagrangianWriter::writeVerts()
{
    // No collectives - can skip on sub-processes
    if (!format_) return;

    const label nVerts = numberOfPoints_;

    // Same payload for connectivity and offsets
    const uint64_t payLoad = vtk::sizeofData<label>(nVerts);


    format().tag(vtk::fileTag::VERTS);

    //
    // 'connectivity'
    // = linear mapping onto points
    //
    {
        format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
        format().writeSize(payLoad);

        vtk::writeIdentity(format(), nVerts);

        format().flush();
        format().endDataArray();
    }

    //
    // 'offsets' (connectivity offsets)
    // = linear mapping onto points (with 1 offset)
    //
    {
        format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
        format().writeSize(payLoad);

        vtk::writeIdentity(format(), nVerts, 1);

        format().flush();
        format().endDataArray();
    }

    format().endTag(vtk::fileTag::VERTS);
}


bool Foam::vtk::lagrangianWriter::beginCellData(label nFields)
{
    return enter_CellData(numberOfPoints_, nFields);
}


bool Foam::vtk::lagrangianWriter::beginPointData(label nFields)
{
    return enter_PointData(numberOfPoints_, nFields);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::lagrangianWriter::lagrangianWriter
(
    const fvMesh& mesh,
    const word& cloudName,
    const vtk::outputOptions opts,
    bool useVerts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts),
    mesh_(mesh),
    cloudName_(cloudName),
    numberOfPoints_(0),
    useVerts_(useVerts)
{
    opts_.append(false); // No append mode (horrible for streaming)
    opts_.legacy(false); // Disallow legacy (see notes)
}


Foam::vtk::lagrangianWriter::lagrangianWriter
(
    const fvMesh& mesh,
    const word& cloudName,
    const fileName& file,
    bool parallel
)
:
    lagrangianWriter(mesh, cloudName)
{
    open(file, parallel);
}


Foam::vtk::lagrangianWriter::lagrangianWriter
(
    const fvMesh& mesh,
    const word& cloudName,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    lagrangianWriter(mesh, cloudName, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::lagrangianWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    // Provide default title

    // XML (inline)

    return vtk::fileWriter::beginFile
    (
        "case='" + mesh_.time().globalCaseName()
      + "' cloud='" + cloudName_
      + "' time='" + mesh_.time().timeName()
      + "' index='" + Foam::name(mesh_.time().timeIndex())
      + "'"
    );
}


bool Foam::vtk::lagrangianWriter::writeGeometry()
{
    enter_Piece();

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }

    // Transcribe cloud into pointField
    pointField cloudPoints(positions());

    numberOfPoints_ = cloudPoints.size();

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
    }


    // beginPiece()
    if (format_)
    {
        if (useVerts_)
        {
            format()
                .tag
                (
                    vtk::fileTag::PIECE,
                    fileAttr::NUMBER_OF_POINTS, numberOfPoints_,
                    fileAttr::NUMBER_OF_VERTS,  numberOfPoints_
                );
        }
        else
        {
            format()
                .tag
                (
                    vtk::fileTag::PIECE,
                    fileAttr::NUMBER_OF_POINTS, numberOfPoints_
                );
        }
    }


    // writePoints()
    {
        if (format_)
        {
            const uint64_t payLoad =
                vtk::sizeofData<float,3>(numberOfPoints_);

            format().tag(vtk::fileTag::POINTS)
                .beginDataArray<float,3>(vtk::dataArrayAttr::POINTS);

            format().writeSize(payLoad);
        }

        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), cloudPoints);
        }
        else
        {
            vtk::writeList(format(), cloudPoints);
        }


        if (format_)
        {
            format().flush();
            format().endDataArray();

            // Non-legacy
            format()
                .endTag(vtk::fileTag::POINTS);
        }
    }

    if (useVerts_) writeVerts();

    return true;
}


bool Foam::vtk::lagrangianWriter::beginParcelData()
{
    if (useVerts_)
    {
        return beginCellData();
    }
    else
    {
        return beginPointData();
    }
}


bool Foam::vtk::lagrangianWriter::endParcelData()
{
    if (useVerts_)
    {
        return endCellData();
    }
    else
    {
        return endPointData();
    }
}


// ************************************************************************* //
