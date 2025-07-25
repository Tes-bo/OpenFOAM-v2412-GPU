/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "lagrangianFieldDecomposer.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::IOField<Type>>
Foam::lagrangianFieldDecomposer::decomposeField
(
    const word& cloudName,
    const IOField<Type>& field
) const
{
    // Create the field for the processor
    return tmp<IOField<Type>>::New
    (
        IOobject
        (
            field.name(),
            procMesh_.time().timeName(),
            cloud::prefix/cloudName,
            procMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        // Mapping internal field values
        Field<Type>(field, particleIndices_)
    );
}


template<class Type>
Foam::tmp<Foam::CompactIOField<Foam::Field<Type>, Type>>
Foam::lagrangianFieldDecomposer::decomposeFieldField
(
    const word& cloudName,
    const CompactIOField<Field<Type>, Type>& field
) const
{
    // Create the field for the processor
    return tmp<CompactIOField<Field<Type>, Type>>::New
    (
        IOobject
        (
            field.name(),
            procMesh_.time().timeName(),
            cloud::prefix/cloudName,
            procMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        // Mapping internal field values
        Field<Field<Type>>(field, particleIndices_)
    );
}


template<class GeoField>
void Foam::lagrangianFieldDecomposer::decomposeFields
(
    const word& cloudName,
    const PtrList<GeoField>& fields
) const
{
    const bool existsOnProc = (particleIndices_.size() > 0);

    for (const GeoField& fld : fields)
    {
        decomposeField(cloudName, fld)().write(existsOnProc);
    }
}


template<class GeoField>
void Foam::lagrangianFieldDecomposer::decomposeFieldFields
(
    const word& cloudName,
    const PtrList<GeoField>& fields
) const
{
    const bool existsOnProc = (particleIndices_.size() > 0);

    for (const GeoField& fld : fields)
    {
        decomposeFieldField(cloudName, fld)().write(existsOnProc);
    }
}


// ************************************************************************* //
