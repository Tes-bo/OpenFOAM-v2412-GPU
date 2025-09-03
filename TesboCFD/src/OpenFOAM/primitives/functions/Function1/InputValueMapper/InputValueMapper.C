/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "InputValueMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum
<
    typename Foam::Function1Types::InputValueMapper<Type>::mappingMode
>
Foam::Function1Types::InputValueMapper<Type>::mappingModeNames_
({
    { mappingMode::NONE, "none" },
    { mappingMode::FUNCTION1, "function" },
    { mappingMode::MINMAX, "minMax" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::InputValueMapper<Type>::read
(
    const dictionary& coeffs
)
{
    mappingMode_ = mappingModeNames_.get("mode", coeffs);

    switch (mappingMode_)
    {
        case mappingMode::NONE:
        {
            break;
        }
        case mappingMode::FUNCTION1:
        {
            mappingValuePtr_.reset
            (
                Function1<scalar>::New("function", coeffs, this->obrPtr_)
            );

            break;
        }
        case mappingMode::MINMAX:
        {
            min_ = coeffs.get<scalar>("min");
            max_ = coeffs.get<scalar>("max");

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << mappingModeNames_[mappingMode_]
                << ".  Available options are: " << mappingModeNames_.sortedToc()
                << abort(FatalError);
        }
    }

    value_ = Function1<Type>::New("value", coeffs, this->obrPtr_);
}


template<class Type>
Foam::Function1Types::InputValueMapper<Type>::InputValueMapper
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, obrPtr),
    mappingMode_(mappingMode::NONE),
    mappingValuePtr_(nullptr),
    min_(0),
    max_(0),
    value_(nullptr)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::InputValueMapper<Type>::InputValueMapper
(
    const InputValueMapper<Type>& rhs
)
:
    Function1<Type>(rhs),
    mappingMode_(rhs.mappingMode_),
    mappingValuePtr_(rhs.mappingValuePtr_.clone()),
    min_(rhs.min_),
    max_(rhs.max_),
    value_(rhs.value_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::InputValueMapper<Type>::writeEntries
(
    Ostream& os
) const
{
    os.writeEntry("mode", mappingModeNames_[mappingMode_]);

    switch (mappingMode_)
    {
        case mappingMode::NONE:
        {
            break;
        }
        case mappingMode::FUNCTION1:
        {
            mappingValuePtr_->writeData(os);

            break;
        }
        case mappingMode::MINMAX:
        {
            os.writeEntry("min", min_);
            os.writeEntry("max", max_);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << mappingModeNames_[mappingMode_]
                << ".  Available options are: " << mappingModeNames_.sortedToc()
                << abort(FatalError);
        }
    }

    value_->writeData(os);
}


template<class Type>
void Foam::Function1Types::InputValueMapper<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
