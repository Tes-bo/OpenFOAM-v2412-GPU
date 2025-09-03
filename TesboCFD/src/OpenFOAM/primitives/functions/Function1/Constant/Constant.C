/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    const Type& value,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, obrPtr),
    value_(value)
{}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    value_(Zero)
{
    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (eptr && eptr->isStream())
    {
        // Primitive (inline) format. Eg,
        // - key constant 1.2;
        // - key 1.2;

        ITstream& is = eptr->stream();
        if (is.peek().isWord())
        {
            is.skip();  // Discard leading 'constant'
        }
        is >> value_;
        dict.checkITstream(is, entryName);
    }
    else
    {
        // Dictionary format. Eg,
        // key { type constant; value 1.2; }

        dict.readEntry("value", value_);
    }
}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    Istream& is
)
:
    Function1<Type>(entryName),
    value_(pTraits<Type>(is))
{}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant(const Constant<Type>& rhs)
:
    Function1<Type>(rhs),
    value_(rhs.value_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1Types::Constant<Type>::value
(
    const scalarField& x
) const
{
    return tmp<Field<Type>>::New(x.size(), value_);
}


template<class Type>
void Foam::Function1Types::Constant<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << value_;
    os.endEntry();
}


// ************************************************************************* //
