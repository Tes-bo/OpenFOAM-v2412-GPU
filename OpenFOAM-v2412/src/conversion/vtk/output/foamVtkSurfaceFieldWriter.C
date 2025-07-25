/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "foamVtkSurfaceFieldWriter.H"
#include "fvsPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Flatten boundary field values into a contiguous list
template<class Type>
static List<Type> flattenBoundary
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{
    const polyBoundaryMesh& pbm = field.mesh().boundaryMesh();

    List<Type> flat(pbm.nFaces(), Foam::zero{});

    forAll(field.boundaryField(), patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const auto& pfld = field.boundaryField()[patchi];

        // Note: restrict transcribing to actual size of the patch field
        // - handles "empty" patch type etc.
        SubList<Type>(flat, pfld.size(), pp.offset()) = pfld;
    }

    return flat;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::surfaceFieldWriter::surfaceFieldWriter
(
    const fvMesh& mesh,
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts),
    mesh_(mesh),
    numberOfPoints_(0)
{
    opts_.append(false); // No append mode (horrible for streaming)
    opts_.legacy(false); // Disallow legacy (inconvenient)
}


Foam::vtk::surfaceFieldWriter::surfaceFieldWriter
(
    const fvMesh& mesh,
    const fileName& file,
    bool parallel
)
:
    surfaceFieldWriter(mesh)
{
    open(file, parallel);
}


Foam::vtk::surfaceFieldWriter::surfaceFieldWriter
(
    const fvMesh& mesh,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    surfaceFieldWriter(mesh, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::surfaceFieldWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }


    // Provide Default title

    if (legacy())
    {
        return vtk::fileWriter::beginFile("surfaceFields");
    }


    // XML (inline)

    return vtk::fileWriter::beginFile
    (
        "surfaceFields "
        "case='" + mesh_.time().globalCaseName()
      + "' region='" + mesh_.name()
      + "' time='" + mesh_.time().timeName()
      + "' index='" + Foam::name(mesh_.time().timeIndex())
      + "'"
    );
}


bool Foam::vtk::surfaceFieldWriter::writeGeometry()
{
    enter_Piece();

    // Output

    const pointField& centres = mesh_.faceCentres();

    // PointData for each face.
    numberOfPoints_ = centres.size();

    if (parallel_)
    {
        reduce(numberOfPoints_, sumOp<label>());
    }

    // <Piece>
    if (format_)
    {
        format()
            .tag
            (
                vtk::fileTag::PIECE,
                fileAttr::NUMBER_OF_POINTS, numberOfPoints_
            );
    }

    // <Point>
    this->beginPoints(numberOfPoints_);

    if (parallel_)
    {
        // Centres for internal faces
        vtk::writeListParallel
        (
            format_.ref(),
            SubList<point>(centres, mesh_.nInternalFaces())
        );

        // Centres for boundary faces
        vtk::writeListParallel
        (
            format_.ref(),
            SubList<point>(centres, mesh_.boundaryMesh().range())
        );
    }
    else
    {
        // Non-parallel: use a normal write

        vtk::writeList(format(), centres);
    }

    this->endPoints();

    return true;
}


bool Foam::vtk::surfaceFieldWriter::beginCellData(label nFields)
{
    // No legacy, no CellData
    return enter_CellData(0, 0);
}


bool Foam::vtk::surfaceFieldWriter::beginPointData(label nFields)
{
    // No legacy
    return enter_PointData(numberOfPoints_, 0);
}


void Foam::vtk::surfaceFieldWriter::write(const surfaceVectorField& field)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    label nFaces = field.mesh().nFaces();

    if (parallel_)
    {
        reduce(nFaces, sumOp<label>());
    }

    if (nFaces != numberOfPoints_)
    {
        FatalErrorInFunction
            << "Expecting " << numberOfPoints_
            << " faces, but found " << nFaces
            << exit(FatalError);
    }

    this->beginDataArray<vector>(field.name(), nFaces);


    // Internal field
    const SubList<vector> internal(field, mesh_.nInternalFaces());

    // Boundary field (flattened)
    auto boundary(flattenBoundary(field));


    if (parallel_)
    {
        // Internal field
        vtk::writeListParallel(format_.ref(), internal);

        // Boundary field
        vtk::writeListParallel(format_.ref(), boundary);
    }
    else
    {
        // Non-parallel

        // Internal field
        vtk::writeList(format_.ref(), internal);

        // Boundary field
        vtk::writeList(format_.ref(), boundary);
    }


    this->endDataArray();
}


// ************************************************************************* //
