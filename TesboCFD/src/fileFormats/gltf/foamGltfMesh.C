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

#include "foamGltfMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::mesh::mesh()
:
    base(),
    fields_(),
    colours_(),
    accessorId_(-1)
{}


Foam::glTF::mesh::mesh(const word& name)
:
    base(name),
    fields_(),
    colours_(),
    accessorId_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label& Foam::glTF::mesh::accessorId() noexcept
{
    return accessorId_;
}


void Foam::glTF::mesh::addField(const word& name, const label accessorId)
{
    fields_.append(Tuple2<string, label>("_field:" + name, accessorId));
}


void Foam::glTF::mesh::addColour(const label accessorId)
{
    colours_.append
    (
        Tuple2<string, label>
        (
            "COLOR_" + Foam::name(colours_.size()),
            accessorId
        )
    );
}


void Foam::glTF::mesh::write(Ostream& os) const
{
    os  << indent << "\"primitives\" : [{" << nl << incrIndent
        << indent << "\"attributes\" : {" << nl  << incrIndent
        << indent << "\"POSITION\" : " << accessorId_;

    for (const auto& f : fields_)
    {
        os  << "," << nl << indent << f.first() << " : " << f.second();
    }

    for (const auto& c : colours_)
    {
        os  << "," << nl << indent << c.first() <<  " : " << c.second();
    }

    os  << nl << decrIndent << indent << "}," << nl
        << indent << "\"mode\" : " << 0 << nl << decrIndent// 0 = POINTS
        << indent << "}]";

    base::write(os);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const glTF::mesh& mesh)
{
    mesh.write(os);

    return os;
}


// ************************************************************************* //
