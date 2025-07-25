/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "CompactListList.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(Istream& is)
:
    offsets_(is),
    values_(is)
{
    // Optionally: enforceSizeSanity();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::CompactListList<T>::readList(Istream& is)
{
    offsets_.readList(is);
    values_.readList(is);
    // Optionally: enforceSizeSanity();
    return is;
}


template<class T>
Foam::Ostream& Foam::CompactListList<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    offsets_.writeList(os, shortLen);
    values_.writeList(os, shortLen);
    return os;
}


template<class T>
Foam::Ostream& Foam::CompactListList<T>::writeMatrix
(
    Ostream& os,
    const label shortLen
) const
{
    const CompactListList<T>& mat = *this;

    const auto oldFmt = os.format(IOstreamOption::ASCII);

    // Write like multi-line matrix of rows (ASCII)
    // TBD: with/without indenting?

    const label nRows = mat.length();

    os  << nRows << nl << token::BEGIN_LIST << nl;

    for (label i = 0; i < nRows; ++i)
    {
        mat.localList(i).writeList(os, shortLen) << nl;
    }

    os  << token::END_LIST << nl;
    os.format(oldFmt);

    return os;
}


// ************************************************************************* //
