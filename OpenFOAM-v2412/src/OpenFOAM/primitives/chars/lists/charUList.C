/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "charList.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace Foam
{

template<>
void UList<char>::writeEntry(Ostream& os) const
{
    os  << word("List<char>");

    if (this->size())
    {
        // Non-zero size: write as binary, so has leading newline separator.
        os << *this;
    }
    else
    {
        // Zero-sized binary - Write size only
        os << token::SPACE << label(0);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<>
void UList<char>::operator=(const Foam::zero)
{
    UList<char>::operator=(char(0));
}

} // End namespace Foam

// ************************************************************************* //
