/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "Square.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Square<Type>::Square
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Sine<Type>(entryName, dict, obrPtr),
    mark_(dict.getOrDefaultCompat<scalar>("mark", {{"markSpace", 2006}}, 1)),
    space_(dict.getOrDefault<scalar>("space", 1))
{}


template<class Type>
Foam::Function1Types::Square<Type>::Square(const Square<Type>& rhs)
:
    Sine<Type>(rhs),
    mark_(rhs.mark_),
    space_(rhs.space_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Square<Type>::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("mark", 1, mark_);
    os.writeEntryIfDifferent<scalar>("space", 1, space_);
    Sine<Type>::writeEntries(os);
}


template<class Type>
void Foam::Function1Types::Square<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
