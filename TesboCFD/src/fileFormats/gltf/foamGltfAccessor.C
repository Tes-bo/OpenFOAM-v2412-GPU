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

#include "foamGltfAccessor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::accessor::accessor()
:
    base(),
    bufferViewId_(-1),
    byteOffset_(0),
    componentType_(-1),
    count_(-1),
    type_(""),
    max_(""),
    min_(""),
    minMax_(false)
{}


Foam::glTF::accessor::accessor(const word& name)
:
    base("Accessor:"+name),
    bufferViewId_(-1),
    byteOffset_(0),
    componentType_(-1),
    count_(-1),
    type_(""),
    max_(""),
    min_(""),
    minMax_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label& Foam::glTF::accessor::bufferViewId() noexcept
{
    return bufferViewId_;
}


Foam::label& Foam::glTF::accessor::byteOffset() noexcept
{
    return byteOffset_;
}


Foam::label& Foam::glTF::accessor::componentType() noexcept
{
    return componentType_;
}


Foam::label& Foam::glTF::accessor::count() noexcept
{
    return count_;
}


Foam::string& Foam::glTF::accessor::type() noexcept
{
    return type_;
}


void Foam::glTF::accessor::write(Ostream& os) const
{
    os  << indent << "\"bufferView\" : " << bufferViewId_ << ',' << nl
        << indent << "\"byteOffset\" : " << byteOffset_ << ',' << nl
        << indent << "\"componentType\" : " << componentType_ << ',' << nl
        << indent << "\"count\" : " << count_ << ',' << nl
        << indent << "\"type\" : " << type_;

    if (minMax_)
    {
        os  << ',' << nl
            << indent << "\"max\" : " << max_.c_str() << ',' << nl
            << indent << "\"min\" : " << min_.c_str();
    }

    base::write(os);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const glTF::accessor& acc)
{
    acc.write(os);

    return os;
}


// ************************************************************************* //
