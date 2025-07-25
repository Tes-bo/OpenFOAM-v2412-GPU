/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include <type_traits>
#include "foamVtkOutput.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtk::fileWriter::beginDataArray
(
    const word& fieldName,
    const label nValues
)
{
    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );

    const direction nCmpt(pTraits<Type>::nComponents);

    if (format_)
    {
        if (std::is_same<label, typename pTraits<Type>::cmptType>::value)
        {
            if (legacy())
            {
                legacy::intField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad =
                    vtk::sizeofData<label, nCmpt>(nValues);

                format().beginDataArray<label, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
        else
        {
            if (legacy())
            {
                legacy::floatField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad =
                    vtk::sizeofData<float, nCmpt>(nValues);

                format().beginDataArray<float, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
    }
}


template<class Type>
void Foam::vtk::fileWriter::writeUniform
(
    const word& fieldName,
    const Type& val,
    const label nLocalValues
)
{
    label nTotal = nLocalValues;

    if (parallel_)
    {
        reduce(nTotal, sumOp<label>());
    }

    this->beginDataArray<Type>(fieldName, nTotal);

    if (parallel_)
    {
        vtk::writeValueParallel(format_.ref(), val, nLocalValues);
    }
    else
    {
        vtk::write(format(), val, nLocalValues);
    }

    this->endDataArray();
}


template<class Type>
void Foam::vtk::fileWriter::writeBasicField
(
    const word& fieldName,
    const UList<Type>& field
)
{
    label nValues = field.size();

    if (parallel_)
    {
        reduce(nValues, sumOp<label>());
    }

    this->beginDataArray<Type>(fieldName, nValues);

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field);
    }
    else
    {
        vtk::writeList(format(), field);
    }

    this->endDataArray();
}


// ************************************************************************* //
