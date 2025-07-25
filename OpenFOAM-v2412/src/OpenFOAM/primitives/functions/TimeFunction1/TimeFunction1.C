/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const Time& runTime,
    const word& entryName,
    const dictionary& dict
)
:
    time_(runTime),
    name_(entryName),
    entry_(Function1<Type>::New(entryName, dict, &runTime))
{
    // Time conversion now handled by Function1 directly
    // entry_->userTimeToTime(runTime);
}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const Time& runTime,
    const word& entryName
)
:
    time_(runTime),
    name_(entryName),
    entry_(nullptr)
{}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const TimeFunction1<Type>& rhs
)
:
    time_(rhs.time_),
    name_(rhs.name_),
    entry_(nullptr) // steal/reuse (missing clone!)
{
    if (rhs.entry_)
    {
        entry_.reset(rhs.entry_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TimeFunction1<Type>::reset(const dictionary& dict)
{
    entry_ = Function1<Type>::New(name_, dict, &time_);
    entry_->userTimeToTime(time_);
}


template<class Type>
const Foam::word& Foam::TimeFunction1<Type>::name() const
{
    return entry_->name();
}


template<class Type>
Type Foam::TimeFunction1<Type>::value(const scalar x) const
{
    return entry_->value(x);
}


template<class Type>
Type Foam::TimeFunction1<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return entry_->integrate(x1, x2);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TimeFunction1<Type>& rhs
)
{
    if (rhs.entry_)
    {
        rhs.entry_->operator<<(os, rhs);
    }
    return os;
}


template<class Type>
void Foam::TimeFunction1<Type>::writeData(Ostream& os) const
{
    entry_->writeData(os);
}


// ************************************************************************* //
