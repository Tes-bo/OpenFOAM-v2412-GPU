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

#include "FunctionObjectValue.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::FunctionObjectValue<Type>::read
(
    const dictionary& coeffs
)
{
    foName_ = coeffs.get<word>("functionObject");
    foResultName_ = coeffs.get<word>("functionObjectResult");
    haveDefaultValue_ = coeffs.readIfPresent("defaultValue", defaultValue_);
}


template<class Type>
Foam::Function1Types::FunctionObjectValue<Type>::FunctionObjectValue
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    foName_(),
    foResultName_(),
    defaultValue_(Zero),
    haveDefaultValue_(false)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::FunctionObjectValue<Type>::FunctionObjectValue
(
    const FunctionObjectValue<Type>& rhs
)
:
    Function1<Type>(rhs),
    foName_(rhs.foName_),
    foResultName_(rhs.foResultName_),
    defaultValue_(rhs.defaultValue_),
    haveDefaultValue_(rhs.haveDefaultValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::FunctionObjectValue<Type>::writeEntries
(
    Ostream& os
) const
{
    os.writeEntry("functionObject", foName_);
    os.writeEntry("functionObjectResult", foResultName_);

    if (haveDefaultValue_)
    {
        os.writeEntry("defaultValue", defaultValue_);
    }
}


template<class Type>
void Foam::Function1Types::FunctionObjectValue<Type>::writeData
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
