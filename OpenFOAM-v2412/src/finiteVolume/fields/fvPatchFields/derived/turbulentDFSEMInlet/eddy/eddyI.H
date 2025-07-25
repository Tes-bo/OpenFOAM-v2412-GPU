/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::eddy::epsi(Random& rndGen) const
{
    // Random number with zero mean and unit variance
    return rndGen.sample01<scalar>() > 0.5 ? 1 : -1;
}


inline Foam::label Foam::eddy::patchFaceI() const noexcept
{
    return patchFaceI_;
}


inline const Foam::point& Foam::eddy::position0() const noexcept
{
    return position0_;
}


inline Foam::scalar Foam::eddy::x() const noexcept
{
    return x_;
}


inline const Foam::vector& Foam::eddy::sigma() const noexcept
{
    return sigma_;
}


inline const Foam::vector& Foam::eddy::alpha() const noexcept
{
    return alpha_;
}


inline const Foam::tensor& Foam::eddy::Rpg() const noexcept
{
    return Rpg_;
}


inline Foam::scalar Foam::eddy::c1() const noexcept
{
    return c1_;
}


inline Foam::point Foam::eddy::position(const vector& n) const
{
    return position0_ + n*x_;
}


Foam::vector Foam::eddy::epsilon(Random& rndGen) const
{
    return vector(epsi(rndGen), epsi(rndGen), epsi(rndGen));
}


inline Foam::scalar Foam::eddy::volume() const
{
    return mathematical::pi*4.0/3.0*cmptProduct(sigma_);
}


void Foam::eddy::move(const scalar dx)
{
    x_ += dx;
}


Foam::boundBox Foam::eddy::bounds(const bool global) const
{
    boundBox bb;

    if (global)
    {
        bb.min() = Rpg_ & -sigma_;
        bb.max() = Rpg_ & sigma_;
    }
    else
    {
        bb.min() = -sigma_;
        bb.max() = sigma_;
    }

    return bb;
}


// ************************************************************************* //
