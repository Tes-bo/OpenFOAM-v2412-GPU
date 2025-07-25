/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EnumType, int nEnum>
inline Foam::label Foam::NamedEnum<EnumType, nEnum>::size() const
{
    return lookup_.size();
}


template<class EnumType, int nEnum>
inline Foam::wordList Foam::NamedEnum<EnumType, nEnum>::toc() const
{
    return lookup_.toc();
}


template<class EnumType, int nEnum>
inline Foam::wordList Foam::NamedEnum<EnumType, nEnum>::sortedToc() const
{
    return lookup_.sortedToc();
}


template<class EnumType, int nEnum>
inline bool Foam::NamedEnum<EnumType, nEnum>::found
(
    const word& enumName
) const
{
    return lookup_.found(enumName);
}


template<class EnumType, int nEnum>
inline bool Foam::NamedEnum<EnumType, nEnum>::hasEnum
(
    const word& enumName
) const
{
    return lookup_.found(enumName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EnumType, int nEnum>
inline const EnumType Foam::NamedEnum<EnumType, nEnum>::operator[]
(
    const word& name
) const
{
    return EnumType(lookup_[name]);
}


template<class EnumType, int nEnum>
inline const char* Foam::NamedEnum<EnumType, nEnum>::operator[]
(
    const EnumType e
) const
{
    return names[int(e)];
}


// ************************************************************************* //
