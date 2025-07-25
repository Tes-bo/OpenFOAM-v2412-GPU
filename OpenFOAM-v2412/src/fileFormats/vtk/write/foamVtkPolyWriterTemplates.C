/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtk::polyWriter::writeUniformValue
(
    const label nCellValues,
    const word& fieldName,
    const Type& val
)
{
    label nValues(0);

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
        nValues = nCellValues;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
        nValues = nLocalPoints_;
    }
    else
    {
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )
            << " for uniform field " << fieldName << nl << endl
            << exit(FatalError);

        return;
    }

    vtk::fileWriter::writeUniform<Type>(fieldName, val, nValues);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::polyWriter::write
(
    const word& fieldName,
    const UList<Type>& field
)
{
    // Could check sizes:
    // const label nValues = field.size();
    // CELL_DATA:   nValues == (nLocalPolys | nLocalLines)
    // POINT_DATA:  nValues == nLocalPoints

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
        return;
    }

    vtk::fileWriter::writeBasicField<Type>(fieldName, field);
}


template<class Type>
void Foam::vtk::polyWriter::writeCellData
(
    const word& fieldName,
    const UList<Type>& field
)
{
    // Could check sizes:
    // const label nValues = field.size();
    // CELL_DATA:   nValues == (nLocalPolys | nLocalLines)

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::CELL_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
        return;
    }

    vtk::fileWriter::writeBasicField<Type>(fieldName, field);
}


template<class Type>
void Foam::vtk::polyWriter::writePointData
(
    const word& fieldName,
    const UList<Type>& field
)
{
    // Could check sizes:
    // const label nValues = field.size();
    // POINT_DATA:  nValues == nLocalPoints

    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
        return;
    }

    vtk::fileWriter::writeBasicField<Type>(fieldName, field);
}


// ************************************************************************* //
