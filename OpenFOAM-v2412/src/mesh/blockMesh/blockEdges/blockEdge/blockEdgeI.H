/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

inline bool Foam::blockEdge::good() const noexcept
{
    return (start_ != end_ && start_ >= 0 && end_ >= 0);
}


inline Foam::label Foam::blockEdge::start() const noexcept
{
    return start_;
}


inline Foam::label Foam::blockEdge::end() const noexcept
{
    return end_;
}


inline const Foam::point& Foam::blockEdge::firstPoint() const
{
    return points_[start_];
}


inline const Foam::point& Foam::blockEdge::lastPoint() const
{
    return points_[end_];
}


inline int Foam::blockEdge::compare(const label start, const label end) const
{
    if (start_ == start && end_ == end)
    {
        return 1;
    }
    else if (start_ == end && end_ == start)
    {
        return -1;
    }

    return 0;
}


inline int Foam::blockEdge::compare(const blockEdge& e) const
{
    return Foam::blockEdge::compare(e.start(), e.end());
}


inline int Foam::blockEdge::compare(const edge& e) const
{
    return Foam::blockEdge::compare(e.start(), e.end());
}


inline Foam::point Foam::blockEdge::linearPosition(const scalar lambda) const
{
    #ifdef FULLDEBUG
    if (lambda < -SMALL || lambda > 1 + SMALL)
    {
        InfoInFunction
            << "Limit parameter to [0-1] range: " << lambda << nl;
    }
    #endif

    if (lambda < SMALL)
    {
        return firstPoint();
    }
    else if (lambda >= 1 - SMALL)
    {
        return lastPoint();
    }

    return firstPoint() + lambda * (lastPoint() - firstPoint());
}


// ************************************************************************* //
