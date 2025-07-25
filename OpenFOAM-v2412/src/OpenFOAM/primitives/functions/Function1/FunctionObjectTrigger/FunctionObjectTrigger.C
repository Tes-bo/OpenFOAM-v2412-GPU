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

#include "FunctionObjectTrigger.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::FunctionObjectTrigger<Type>::read
(
    const dictionary& coeffs
)
{
    triggers_ = coeffs.get<labelList>("triggers", keyType::LITERAL);
    defaultValue_ =
        coeffs.getOrDefault("defaultValue", false, keyType::LITERAL);
}


template<class Type>
Foam::Function1Types::FunctionObjectTrigger<Type>::FunctionObjectTrigger
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    triggers_(),
    defaultValue_(false)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::FunctionObjectTrigger<Type>::FunctionObjectTrigger
(
    const FunctionObjectTrigger<Type>& rhs
)
:
    Function1<Type>(rhs),
    triggers_(rhs.triggers_),
    defaultValue_(rhs.defaultValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::FunctionObjectTrigger<Type>::writeEntries
(
    Ostream& os
) const
{
    os.writeKeyword("triggers");
    os << flatOutput(triggers_);
    os.endEntry();

    if (defaultValue_)
    {
        os.writeEntry("default", "true");
    }
}


template<class Type>
void Foam::Function1Types::FunctionObjectTrigger<Type>::writeData
(
    Ostream& os
) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
