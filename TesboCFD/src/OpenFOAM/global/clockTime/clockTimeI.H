/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::clockTime::clockTime()
:
    start_(true),  // == clockValue::now()
    last_(start_)
{}


inline Foam::clockTime::clockTime(const clockValue& clockval)
:
    start_(clockval),
    last_(start_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::clockTime::resetTime()
{
    last_.update();
    start_ = last_;
}


inline void Foam::clockTime::resetTimeIncrement() const
{
    last_.update();
}


inline double Foam::clockTime::elapsedTime() const
{
    last_.update();
    return static_cast<double>(last_ - start_);
}


inline double Foam::clockTime::timeIncrement() const
{
    const auto prev(last_);
    last_.update();
    return static_cast<double>(last_ - prev);
}


// ************************************************************************* //
