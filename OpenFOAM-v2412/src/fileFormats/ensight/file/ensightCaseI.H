/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::ensightCase::options& Foam::ensightCase::option() const
{
    return *options_;
}


inline Foam::IOstreamOption::streamFormat Foam::ensightCase::format() const
{
    return options_->format();
}


inline const Foam::fileName& Foam::ensightCase::path() const noexcept
{
    return ensightDir_;
}


inline const Foam::word& Foam::ensightCase::mask() const
{
    return options_->mask();
}


inline Foam::word Foam::ensightCase::padded(const label i) const
{
    return options_->padded(i);
}


inline bool Foam::ensightCase::nodeValues() const
{
    return options_->nodeValues();
}


inline bool Foam::ensightCase::separateCloud() const
{
    return options_->separateCloud();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::ensightCase::operator()() const
{
    return (os_ ? *os_ : Foam::Snull);
}


// ************************************************************************* //
