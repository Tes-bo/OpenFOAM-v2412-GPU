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

#include "foamGltfBufferView.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::bufferView::bufferView()
:
    base(),
    buffer_(0),
    byteOffset_(-1),
    byteLength_(-1),
    target_(-1)
{}


Foam::glTF::bufferView::bufferView(const word& name)
:
    base("BufferView:"+name),
    buffer_(0),
    byteOffset_(-1),
    byteLength_(-1),
    target_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label& Foam::glTF::bufferView::buffer() noexcept
{
    return buffer_;
}


Foam::label& Foam::glTF::bufferView::byteOffset() noexcept
{
    return byteOffset_;
}


Foam::label& Foam::glTF::bufferView::byteLength() noexcept
{
    return byteLength_;
}


Foam::label& Foam::glTF::bufferView::target() noexcept
{
    return target_;
}


void Foam::glTF::bufferView::write(Ostream& os) const
{
    os  << indent << "\"buffer\" : " << buffer_ << "," << nl
        << indent << "\"byteOffset\" : " << byteOffset_ << "," << nl
        << indent << "\"byteLength\" : " << byteLength_;

    if (target_ != -1)
    {
        os  << "," << nl
            << indent << "\"target\" : " << target_;
    }

    base::write(os);
}


void Foam::glTF::bufferView::operator=(const glTF::bufferView& bv)
{
    base::operator=(bv);

    buffer_ = bv.buffer_;
    byteOffset_ = bv.byteOffset_;
    byteLength_ = bv.byteLength_;
    target_ = bv.target_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const glTF::bufferView& bv)
{
    bv.write(os);

    return os;
}


// ************************************************************************* //
