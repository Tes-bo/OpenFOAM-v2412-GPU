/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "patchSummaryTemplates.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
Foam::PtrList<GeoField> Foam::readFields
(
    const IOobjectList& objects,
    const typename GeoField::Mesh& mesh
)
{
    const UPtrList<const IOobject> fieldObjects
    (
        objects.csorted<GeoField>()
    );

    PtrList<GeoField> fields(fieldObjects.size());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (!nFields)
        {
            Info<< "    " << GeoField::typeName << " (";
        }
        else
        {
            Info<< ' ';
        }
        Info<< io.name();

        fields.emplace_set(nFields, io, mesh);
        ++nFields;
    }

    if (nFields)
    {
        Info<< ')' << nl;
    }

    return fields;
}


template<class GeoField>
void Foam::outputFieldList
(
    const UPtrList<GeoField>& fieldList,
    const label patchi
)
{
    forAll(fieldList, fieldi)
    {
        if (fieldList.set(fieldi))
        {
            Info<< "    " << pTraits<typename GeoField::value_type>::typeName
                << tab << tab
                << fieldList[fieldi].name() << tab << tab
                << fieldList[fieldi].boundaryField()[patchi].type() << nl;
        }
    }
}


template<class GeoField>
void Foam::collectFieldList
(
    const UPtrList<GeoField>& fieldList,
    const label patchi,
    HashTable<word>& fieldToType
)
{
    forAll(fieldList, fieldi)
    {
        if (fieldList.set(fieldi))
        {
            fieldToType.insert
            (
                fieldList[fieldi].name(),
                fieldList[fieldi].boundaryField()[patchi].type()
            );
        }
    }
}


// ************************************************************************* //
