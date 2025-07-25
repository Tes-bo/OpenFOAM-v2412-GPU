/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

\*---------------------------------------------------------------------------*/

#include "ShapeSensitivitiesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void ShapeSensitivitiesBase::constructAndWriteSensitivityField
(
    const autoPtr
    <
        typename GeometricField<Type, fvPatchField, volMesh>::Boundary
    >& sensFieldPtr,
    const word& name
) const
{
    GeometricField<Type, fvPatchField, volMesh> volSensField
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>(dimless, Zero)
    );

    for (const label patchI : sensitivityPatchIDs_)
    {
        volSensField.boundaryFieldRef()[patchI] = sensFieldPtr()[patchI];
    }

    volSensField.write();
}


template<class Type>
void ShapeSensitivitiesBase::constructAndWriteSensitivtyPointField
(
    const autoPtr<List<Field<Type>>>& sensFieldPtr,
    const word& name
) const
{
    GeometricField<Type, pointPatchField, pointMesh> pointSensField
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh_),
        dimensioned<Type>(dimless, Zero)
        //fixedValuePointPatchField<Type>::typeName
    );

    for (const label patchI : sensitivityPatchIDs_)
    {
        // Paraview does not visualise values on fixedValuePointPatchFields
        // Only option is to store in internalField
        //pointSensField.boundaryFieldRef()[patchI] == sensFieldPtr()[patchI];
        pointSensField.boundaryField()[patchI].setInInternalField
        (
            pointSensField.primitiveFieldRef(),
            sensFieldPtr()[patchI]
        );
    }

    pointSensField.write();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ShapeSensitivitiesBase::constructVolSensitivtyField
(
    const autoPtr
    <
        typename GeometricField<Type, fvPatchField, volMesh>::Boundary
    >& sensFieldPtr,
    const word& name
) const
{
        tmp<GeometricField<Type, fvPatchField, volMesh>> tVolSensField
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                pTraits<Type>::zero
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& volSensField =
            tVolSensField.ref();

        typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
            volSensFieldbf = volSensField.boundaryFieldRef();

        forAll(sensitivityPatchIDs_, pI)
        {
            const label patchI = sensitivityPatchIDs_[pI];
            volSensFieldbf[patchI] = sensFieldPtr()[patchI];
        }

        return tVolSensField;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
