/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "fvMeshTools.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshTools::addPatchFields
(
    fvMesh& mesh,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const typename GeoField::value_type& defaultPatchValue
)
{
    for (GeoField& fld : mesh.objectRegistry::sorted<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        const label newPatchi = bfld.size();
        bfld.resize(newPatchi+1);

        const dictionary* dict = patchFieldDict.findDict(fld.name());

        if (dict)
        {
            bfld.set
            (
                newPatchi,
                GeoField::Patch::New
                (
                    mesh.boundary()[newPatchi],
                    fld.internalField(),
                    *dict
                )
            );
        }
        else
        {
            bfld.set
            (
                newPatchi,
                GeoField::Patch::New
                (
                    defaultPatchFieldType,
                    mesh.boundary()[newPatchi],
                    fld.internalField()
                )
            );
            bfld[newPatchi] == defaultPatchValue;
        }
    }
}


template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    for (GeoField& fld : mesh.objectRegistry::sorted<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        const dictionary* dict = patchFieldDict.findDict(fld.name());

        if (dict)
        {
            bfld.set
            (
                patchi,
                GeoField::Patch::New
                (
                    mesh.boundary()[patchi],
                    fld.internalField(),
                    *dict
                )
            );
        }
    }
}


template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    fvMesh& mesh,
    const label patchi,
    const typename GeoField::value_type& value
)
{
    for (GeoField& fld : mesh.objectRegistry::sorted<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        bfld[patchi] == value;
    }
}


// Remove last patch field
template<class GeoField>
void Foam::fvMeshTools::trimPatchFields(fvMesh& mesh, const label nPatches)
{
    for (GeoField& fld : mesh.objectRegistry::sorted<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        bfld.resize(nPatches);
    }
}


// Reorder patch field
template<class GeoField>
void Foam::fvMeshTools::reorderPatchFields
(
    fvMesh& mesh,
    const labelList& oldToNew
)
{
    for (GeoField& fld : mesh.objectRegistry::sorted<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        bfld.reorder(oldToNew);
    }
}


// ************************************************************************* //
