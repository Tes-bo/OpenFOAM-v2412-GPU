/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::cellSetOption::timeStart() const noexcept
{
    return timeStart_;
}


inline Foam::scalar Foam::fv::cellSetOption::duration() const noexcept
{
    return duration_;
}


inline bool Foam::fv::cellSetOption::inTimeLimits(const scalar timeValue) const
{
    return
    (
        (timeStart_ < 0)
     ||
        (
            (mesh_.time().value() >= timeStart_)
         && (mesh_.time().value() <= (timeStart_ + duration_))
        )
    );
}


inline bool Foam::fv::cellSetOption::useSubMesh() const noexcept
{
    return selectionMode_ != selectionModeType::smAll;
}


inline Foam::fv::cellSetOption::selectionModeType
Foam::fv::cellSetOption::selectionMode() const noexcept
{
    return selectionMode_;
}


inline const Foam::wordRes&
Foam::fv::cellSetOption::selectionNames() const noexcept
{
    return selectionNames_;
}


inline const Foam::wordRe& Foam::fv::cellSetOption::zoneName() const
{
    return (selectionNames_.empty() ? wordRe::null : selectionNames_.first());
}


inline Foam::scalar Foam::fv::cellSetOption::V() const noexcept
{
    return V_;
}


inline const Foam::labelList& Foam::fv::cellSetOption::cells() const noexcept
{
    return cells_;
}


inline Foam::scalar Foam::fv::cellSetOption::timeStart(scalar val) noexcept
{
    scalar old(timeStart_);
    timeStart_ = val;
    return old;
}


inline Foam::scalar Foam::fv::cellSetOption::duration(scalar val) noexcept
{
    scalar old(duration_);
    duration_ = val;
    return old;
}


// ************************************************************************* //
