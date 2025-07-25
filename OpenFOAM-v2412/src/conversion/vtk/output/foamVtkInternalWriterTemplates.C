/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "foamVtkInternalWriter.H"
#include "foamVtkOutput.H"
#include "volPointInterpolation.H"
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
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

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();


    this->beginDataArray<Type>(field.name(), numberOfPoints_);

    if (parallel_)
    {
        List<Type> addedValues(addPointCellLabels.size());
        label outi = 0;

        for (const label cellId : addPointCellLabels)
        {
            addedValues[outi++] = interpolatePointToCell(field, cellId);
        }

        vtk::writeListsParallel(format_.ref(), field, addedValues);
    }
    else
    {
        vtk::writeList(format(), field);

        for (const label cellId : addPointCellLabels)
        {
            const Type val = interpolatePointToCell(field, cellId);
            vtk::write(format(), val);
        }
    }

    this->endDataArray();
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const DimensionedField<Type, volMesh>& field
)
{
    writeCellData(field.name(), field.field());
}


template<class Type, template<class> class PatchField>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, PatchField, volMesh>& field
)
{
    writeCellData(field.name(), field.primitiveField());
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const DimensionedField<Type, volMesh>& vfield,
    const volPointInterpolation& pInterp
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << vfield.name() << nl << endl
            << exit(FatalError);
    }

    typedef DimensionedField<Type, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const auto& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();


    this->beginDataArray<Type>(vfield.name(), numberOfPoints_);

    if (parallel_)
    {
        vtk::writeListsParallel
        (
            format_.ref(),
            pfield,
            vfield,
            addPointCellLabels
        );
    }
    else
    {
        vtk::writeLists(format(), pfield, vfield, addPointCellLabels);
    }

    this->endDataArray();
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, fvPatchField, volMesh>& vfield,
    const volPointInterpolation& pInterp
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << vfield.name() << nl << endl
            << exit(FatalError);
    }

    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const auto& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();


    this->beginDataArray<Type>(vfield.name(), numberOfPoints_);

    if (parallel_)
    {
        vtk::writeListsParallel
        (
            format_.ref(),
            pfield,
            vfield,
            addPointCellLabels
        );
    }
    else
    {
        vtk::writeLists(format(), pfield, vfield, addPointCellLabels);
    }

    this->endDataArray();
}


// ************************************************************************* //
