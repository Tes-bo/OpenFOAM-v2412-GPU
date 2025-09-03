/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Lookup<Type>::Lookup
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    lookupBase(dict),
    Function1<Type>(entryName, dict, obrPtr)
{}


template<class Type>
Foam::Function1Types::Lookup<Type>::Lookup(const Lookup<Type>& rhs)
:
    lookupBase(rhs),
    Function1<Type>(rhs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Lookup<Type>::value(const scalar t) const
{
    const objectRegistry& db = function1Base::obr();
    const auto& obj =
        db.lookupObject<UniformDimensionedField<Type>>(name_, true);

    return obj.value();
}


template<class Type>
Type Foam::Function1Types::Lookup<Type>::integrate
(
    const scalar t1,
    const scalar t2
) const
{
    return (t2-t1)*value(0.5*(t1+t2));
}


template<class Type>
void Foam::Function1Types::Lookup<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
