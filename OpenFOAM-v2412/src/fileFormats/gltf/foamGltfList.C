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

#include "foamGltfList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type& Foam::glTF::List<Type>::create(const word& name)
{
    const label id = data_.size();
    Type& obj = data_.emplace_back(name);
    obj.id() = id;

    return data_.back();
}


template<class Type>
const Foam::DynamicList<Type>& Foam::glTF::List<Type>::data() const noexcept
{
    return data_;
}


template<class Type>
bool Foam::glTF::List<Type>::empty() const noexcept
{
    return data_.empty();
}


template<class Type>
Foam::label Foam::glTF::List<Type>::size() const noexcept
{
    return data_.size();
}


template<class Type>
void Foam::glTF::List<Type>::write
(
    Ostream& os,
    const word& keyword,
    bool firstEntry
)
{
    if (empty())
    {
        return;
    }

    if (!firstEntry)
    {
        os  << "," << nl;
    }

    os  << indent << "\"" << keyword << "\" : [" << nl << incrIndent;

    write(os);

    os  << decrIndent << nl << indent << "]";
}


template<class Type>
void Foam::glTF::List<Type>::write(Ostream& os) const
{
    forAll(data_, i)
    {
        os  << indent << "{"
            << nl << incrIndent
            << data_[i]
            << nl << decrIndent
            << indent << "}";

        if (i != data_.size() - 1) os  << "," << nl;
    }
}


template<class Type>
Type& Foam::glTF::List<Type>::operator[](const label i)
{
    return data_[i];
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const glTF::List<Type>& lst)
{
    lst.write(os);

    return os;
}


// ************************************************************************* //
