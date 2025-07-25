/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::surfaceWriterMethods

Description
    Convenience macros for instantiating surfaceWriter methods.

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriterMethods_H
#define Foam_surfaceWriterMethods_H

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiate templated method for standard types
#define defineSurfaceWriterWriteField(ThisClass, FieldType)                    \
    Foam::fileName ThisClass::write                                            \
    (                                                                          \
        const word& fieldName,                                                 \
        const Field<FieldType>& values                                         \
    )                                                                          \
    {                                                                          \
        return writeTemplate(fieldName, values);                               \
    }


#define defineSurfaceWriterWriteFields(ThisClass)                              \
    defineSurfaceWriterWriteField(ThisClass, label);                           \
    defineSurfaceWriterWriteField(ThisClass, scalar);                          \
    defineSurfaceWriterWriteField(ThisClass, vector);                          \
    defineSurfaceWriterWriteField(ThisClass, sphericalTensor);                 \
    defineSurfaceWriterWriteField(ThisClass, symmTensor);                      \
    defineSurfaceWriterWriteField(ThisClass, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
