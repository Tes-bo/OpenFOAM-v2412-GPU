/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "adiosReader.H"
#include "areaFields.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    static inline bool warnUnhandled(const regIOobject& obj, bool carp)
    {
        if (carp)
        {
            Warning
                << "No adios field input for "
                << obj.name() << " type: " << obj.type() << nl;
        }

        return false;
    }


    //- Try dynamic_cast to TargetType (non-cast variant)
    template<class TargetType, class Type>
    static inline TargetType* isA(Type& t)
    {
        return dynamic_cast<TargetType*>(&t);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiosFoam::adiosReader::readFieldObject
(
    regIOobject& obj,
    const adiosFoam::fieldInfo& src,
    const bool carp
)
{
    const word& fieldType = obj.type();
    const word&   srcType = src.type();

    if (fieldType != srcType)
    {
        // probably fatal:
        Info<<"WARNING mismatch on field " << src
            << " (expected: " << fieldType
            << " but had " << srcType << ")\n";
        return false;
    }

    Info<<"Read " << obj.name() << " (type " << fieldType << ") from adios\n";

    const fileName varName(src.fullName());

    return
    (
        // Point fields
        readGeometricField(isA<pointScalarField>(obj), varName)
     || readGeometricField(isA<pointVectorField>(obj), varName)
     || readGeometricField(isA<pointSphericalTensorField>(obj), varName)
     || readGeometricField(isA<pointSymmTensorField>(obj), varName)
     || readGeometricField(isA<pointTensorField>(obj), varName)

        // Surface fields
     || readGeometricField(isA<surfaceScalarField>(obj), varName)
     || readGeometricField(isA<surfaceVectorField>(obj), varName)
     || readGeometricField(isA<surfaceSphericalTensorField>(obj), varName)
     || readGeometricField(isA<surfaceSymmTensorField>(obj), varName)
     || readGeometricField(isA<surfaceTensorField>(obj), varName)

        // Area fields
     || readGeometricField(isA<areaScalarField>(obj), varName)
     || readGeometricField(isA<areaVectorField>(obj), varName)
     || readGeometricField(isA<areaSphericalTensorField>(obj), varName)
     || readGeometricField(isA<areaSymmTensorField>(obj), varName)
     || readGeometricField(isA<areaTensorField>(obj), varName)

        // Volume fields
     || readGeometricField(isA<volScalarField>(obj), varName)
     || readGeometricField(isA<volVectorField>(obj), varName)
     || readGeometricField(isA<volSphericalTensorField>(obj), varName)
     || readGeometricField(isA<volSymmTensorField>(obj), varName)
     || readGeometricField(isA<volTensorField>(obj), varName)

        // Internal volume fields
     || readField(isA<volScalarField::Internal>(obj), varName)
     || readField(isA<volVectorField::Internal>(obj), varName)
     || readField(isA<volSphericalTensorField::Internal>(obj), varName)
     || readField(isA<volSymmTensorField::Internal>(obj), varName)
     || readField(isA<volTensorField::Internal>(obj), varName)

     || warnUnhandled(obj, carp)
    );
}


// ************************************************************************* //
