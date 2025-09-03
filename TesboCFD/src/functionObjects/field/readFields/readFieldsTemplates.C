/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::functionObjects::readFields::loadAndStore(const IOobject& io)
{
    if (io.isHeaderClass<FieldType>())
    {
        // Store field on mesh database
        Log << "    Reading " << io.name()
            << " (" << io.headerClassName() << ')' << endl;

        regIOobject::store(new FieldType(io, mesh_));
        return true;
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::readFields::loadField(const IOobject& io)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    return
    (
        loadAndStore<VolFieldType>(io)
     || loadAndStore<IntVolFieldType>(io)
     || loadAndStore<SurfaceFieldType>(io)
    );
}


// ************************************************************************* //
