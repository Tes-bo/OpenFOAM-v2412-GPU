/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "fvcDiv.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FieldType>
bool Foam::functionObjects::div::calcDiv()
{
    const auto* fieldptr = cfindObject<FieldType>(fieldName_);

    if (!fieldptr)
    {
        return false;
    }

    if (zoneSubSetPtr_)
    {
        const fvMeshSubset& subsetter = zoneSubSetPtr_->subsetter();

        return storeInDb
        (
            resultName_,
            fvc::div(subsetter.interpolate(*fieldptr, false)),
            subsetter.subMesh().thisDb()
        );
    }


    return store
    (
        resultName_,
        fvc::div(*fieldptr)
    );
}


template<class Type>
bool Foam::functionObjects::div::writeField()
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& subMesh = zoneSubSetPtr_->subsetter().subMesh();
    const auto* fieldptr = subMesh.findObject<volFieldType>(resultName_);

    if (fieldptr)
    {
        tmp<volFieldType> tfield = zoneSubSetPtr_->mapToZone<Type>(*fieldptr);
        tfield().write();

        return true;
    }

    return false;
}


// ************************************************************************* //
