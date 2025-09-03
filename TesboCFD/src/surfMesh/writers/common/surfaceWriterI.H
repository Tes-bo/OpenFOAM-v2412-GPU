/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

inline bool Foam::surfaceWriter::is_open() const noexcept
{
    return !outputPath_.empty();
}


inline Foam::label Foam::surfaceWriter::nFields() const noexcept
{
    return nFields_;
}


inline Foam::label Foam::surfaceWriter::nFields(const label n) noexcept
{
    label old(nFields_);
    nFields_ = n;
    return old;
}


inline bool Foam::surfaceWriter::isPointData() const noexcept
{
    return isPointData_;
}


inline bool Foam::surfaceWriter::isPointData(const bool on) noexcept
{
    bool old(isPointData_);
    isPointData_ = on;
    return old;
}


inline bool Foam::surfaceWriter::useTimeDir() const noexcept
{
    return useTimeDir_;
}


inline bool Foam::surfaceWriter::useTimeDir(const bool on) noexcept
{
    bool old(useTimeDir_);
    useTimeDir_ = on;
    return old;
}


inline bool Foam::surfaceWriter::verbose() const noexcept
{
    return verbose_;
}


inline bool Foam::surfaceWriter::verbose(const bool on) noexcept
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


inline Foam::scalar Foam::surfaceWriter::mergeDim() const noexcept
{
    return mergeDim_;
}


inline Foam::scalar Foam::surfaceWriter::mergeDim(const scalar dist) noexcept
{
    scalar old(mergeDim_);
    mergeDim_ = dist;
    return old;
}


inline Foam::scalar Foam::surfaceWriter::scale() const noexcept
{
    return geometryScale_;
}


inline Foam::scalar Foam::surfaceWriter::scale(const scalar factor) noexcept
{
    // This is probably not yet needed -> adjustedSurf_.clear();
    scalar old(geometryScale_);
    geometryScale_ = factor;
    return old;
}


inline const Foam::coordSystem::cartesian&
Foam::surfaceWriter::transform() const noexcept
{
    return geometryTransform_;
}


inline bool Foam::surfaceWriter::hasTime() const
{
    return currTime_.name().size();
}


inline const Foam::word& Foam::surfaceWriter::timeName() const
{
    return currTime_.name();
}


inline Foam::scalar Foam::surfaceWriter::timeValue() const
{
    return currTime_.name().empty() ? 0 : currTime_.value();
}


// ************************************************************************* //
