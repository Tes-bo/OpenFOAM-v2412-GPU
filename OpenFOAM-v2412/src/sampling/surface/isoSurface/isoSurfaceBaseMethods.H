/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::isoSurfaceBaseMethods

Description
    Convenience macros for instantiating iso-surface interpolate methods.

\*---------------------------------------------------------------------------*/

#ifndef Foam_isoSurfaceBaseMethods_H
#define Foam_isoSurfaceBaseMethods_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiate templated method for standard types
// See note in isoSurfaceBase explaining why these are volume fields

#undef  defineIsoSurfaceInterpolateMethod
#define defineIsoSurfaceInterpolateMethod(ThisClass, Type)                     \
    Foam::tmp<Foam::Field<Type>> ThisClass::interpolate                        \
    (                                                                          \
        const VolumeField<Type>& cellValues,                                   \
        const Field<Type>& pointValues                                         \
    ) const                                                                    \
    {                                                                          \
        return interpolateTemplate(cellValues, pointValues);                   \
    }


#define defineIsoSurfaceInterpolateMethods(ThisClass)                          \
    defineIsoSurfaceInterpolateMethod(ThisClass, Foam::scalar);                \
    defineIsoSurfaceInterpolateMethod(ThisClass, Foam::vector);                \
    defineIsoSurfaceInterpolateMethod(ThisClass, Foam::sphericalTensor);       \
    defineIsoSurfaceInterpolateMethod(ThisClass, Foam::symmTensor);            \
    defineIsoSurfaceInterpolateMethod(ThisClass, Foam::tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
