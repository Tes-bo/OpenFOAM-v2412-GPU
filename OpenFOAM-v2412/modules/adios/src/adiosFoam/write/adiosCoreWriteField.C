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

#include "adiosCoreWrite.H"
#include "areaFields.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    static inline bool warnUnhandled(const regIOobject& obj, bool carp)
    {
        if (carp)
        {
            Warning
                << "No adios field output for "
                << obj.name() << " type: " << obj.type() << nl;
        }

        return false;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::adiosFoam::adiosCoreWrite::writeFieldObject
(
    const regIOobject& obj,
    const bool carp
)
{
    return
    (
        // Point fields
        writeGeometricField(isA<pointScalarField>(obj))
     || writeGeometricField(isA<pointVectorField>(obj))
     || writeGeometricField(isA<pointSphericalTensorField>(obj))
     || writeGeometricField(isA<pointSymmTensorField>(obj))
     || writeGeometricField(isA<pointTensorField>(obj))

        // Surface fields
     || writeGeometricField(isA<surfaceScalarField>(obj))
     || writeGeometricField(isA<surfaceVectorField>(obj))
     || writeGeometricField(isA<surfaceSphericalTensorField>(obj))
     || writeGeometricField(isA<surfaceSymmTensorField>(obj))
     || writeGeometricField(isA<surfaceTensorField>(obj))

        // Area fields
     || writeGeometricField(isA<areaScalarField>(obj))
     || writeGeometricField(isA<areaVectorField>(obj))
     || writeGeometricField(isA<areaSphericalTensorField>(obj))
     || writeGeometricField(isA<areaSymmTensorField>(obj))
     || writeGeometricField(isA<areaTensorField>(obj))

        // Volume fields
     || writeGeometricField(isA<volScalarField>(obj))
     || writeGeometricField(isA<volVectorField>(obj))
     || writeGeometricField(isA<volSphericalTensorField>(obj))
     || writeGeometricField(isA<volSymmTensorField>(obj))
     || writeGeometricField(isA<volTensorField>(obj))

        // Internal volume fields
     || writeInternalField(isA<volScalarField::Internal>(obj))
     || writeInternalField(isA<volVectorField::Internal>(obj))
     || writeInternalField(isA<volSphericalTensorField::Internal>(obj))
     || writeInternalField(isA<volSymmTensorField::Internal>(obj))
     || writeInternalField(isA<volTensorField::Internal>(obj))

     || warnUnhandled(obj, carp)
    );
}


// ************************************************************************* //
