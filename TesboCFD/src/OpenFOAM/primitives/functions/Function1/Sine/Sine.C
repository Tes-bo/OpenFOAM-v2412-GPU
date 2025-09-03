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

#include "Sine.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Sine<Type>::Sine
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    t0_(dict.getOrDefault<scalar>("t0", 0)),
    amplitude_
    (
        Function1<scalar>::NewIfPresent("amplitude", dict, word::null, obrPtr)
    ),
    period_
    (
        Function1<scalar>::NewIfPresent("period", dict, word::null, obrPtr)
    ),
    frequency_(nullptr),
    scale_(Function1<Type>::New("scale", dict, obrPtr)),
    level_(Function1<Type>::New("level", dict, obrPtr))
{
    if (!period_)
    {
        frequency_ = Function1<scalar>::New("frequency", dict, obrPtr);
    }
}


template<class Type>
Foam::Function1Types::Sine<Type>::Sine(const Sine<Type>& rhs)
:
    Function1<Type>(rhs),
    t0_(rhs.t0_),
    amplitude_(rhs.amplitude_.clone()),
    period_(rhs.period_.clone()),
    frequency_(rhs.frequency_.clone()),
    scale_(rhs.scale_.clone()),
    level_(rhs.level_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sine<Type>::userTimeToTime(const Time& t)
{
    t0_ = t.userTimeToTime(t0_);
}


template<class Type>
void Foam::Function1Types::Sine<Type>::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("t0", 0, t0_);
    if (amplitude_)
    {
        amplitude_->writeData(os);
    }
    if (period_)
    {
        period_->writeData(os);
    }
    if (frequency_)
    {
        frequency_->writeData(os);
    }
    scale_->writeData(os);
    level_->writeData(os);
}


template<class Type>
void Foam::Function1Types::Sine<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
