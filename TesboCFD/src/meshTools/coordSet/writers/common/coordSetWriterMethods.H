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

InClass
    Foam::coordSetWriterMethods

Description
    Convenience macros for instantiating coordSetWriter methods.

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriterMethods_H
#define Foam_coordSetWriterMethods_H

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiate templated method for standard types
#define defineCoordSetWriterWriteField(ThisClass, ValueType)                  \
                                                                              \
    Foam::fileName ThisClass::write                                           \
    (                                                                         \
        const word& fieldName,                                                \
        const Field<ValueType>& values                                        \
    )                                                                         \
    {                                                                         \
        return writeTemplate(fieldName, values);                              \
    }                                                                         \
                                                                              \
    Foam::fileName ThisClass::write                                           \
    (                                                                         \
        const word& fieldName,                                                \
        const List<Field<ValueType>>& fieldValues                             \
    )                                                                         \
    {                                                                         \
        return writeTemplate(fieldName, fieldValues);                         \
    }


#define defineCoordSetWriterWriteFields(ThisClass)                            \
    defineCoordSetWriterWriteField(ThisClass, label);                         \
    defineCoordSetWriterWriteField(ThisClass, scalar);                        \
    defineCoordSetWriterWriteField(ThisClass, vector);                        \
    defineCoordSetWriterWriteField(ThisClass, sphericalTensor);               \
    defineCoordSetWriterWriteField(ThisClass, symmTensor);                    \
    defineCoordSetWriterWriteField(ThisClass, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
