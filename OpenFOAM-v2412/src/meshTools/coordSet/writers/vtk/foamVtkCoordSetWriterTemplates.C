/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::coordSetWriter::writePointData
(
    const word& fieldName,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    // Could check sizes:

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

    label nValues = 0;

    for (const Field<Type>& field : fieldPtrs)
    {
        nValues += field.size();
    }

    // if (parallel_)
    // {
    //     reduce(nValues, sumOp<label>());
    // }

    this->beginDataArray<Type>(fieldName, nValues);

    // if (parallel_)
    // {
    //     vtk::writeListParallel(format_.ref(), field);
    // }
    // else
    {
        for (const Field<Type>& field : fieldPtrs)
        {
            vtk::writeList(format(), field);
        }
    }

    this->endDataArray();
}


// ************************************************************************* //
