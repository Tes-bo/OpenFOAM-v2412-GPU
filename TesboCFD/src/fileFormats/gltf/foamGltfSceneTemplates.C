/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

template<class Type>
Foam::label Foam::glTF::scene::addField
(
    const Type& fld,
    const word& name,
    const label target
)
{
    const direction nCmpts = pTraits<typename Type::value_type>::nComponents;

    auto& bv = bufferViews_.create(name);
    bv.byteOffset() = bytes_;
    bv.byteLength() = fld.size()*nCmpts*sizeof(float);
    if (target != -1)
    {
        bv.target() = target;
    }
    bytes_ += bv.byteLength();

    auto& acc = accessors_.create(name);
    acc.bufferViewId() = bv.id();
    acc.set(fld);

    auto& obj = objects_.create(name);
    obj.addData(fld);

    return acc.id();
}


template<class Type>
Foam::label Foam::glTF::scene::addMesh(const Type& fld, const word& name)
{
    const label accessorId =
        addField(fld, name, key(targetTypes::ARRAY_BUFFER));

    auto& gmesh = meshes_.create(name);
    gmesh.accessorId() = accessorId;

    return meshes_.size() - 1;
}


template<class Type>
Foam::label Foam::glTF::scene::addFieldToMesh
(
    const Type& fld,
    const word& name,
    const label meshi
)
{
    auto& gmesh = getMesh(meshi);

    const label accessorId = addField(fld, name);

    gmesh.addField(name, accessorId);

    return accessorId;
}


// ************************************************************************* //
