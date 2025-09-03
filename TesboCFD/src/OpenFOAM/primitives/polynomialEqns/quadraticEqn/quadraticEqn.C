/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "linearEqn.H"
#include "quadraticEqn.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Roots<2> Foam::quadraticEqn::roots() const
{
    const scalar a = this->a();
    const scalar b = this->b();
    const scalar c = this->c();

    // Check the leading term in the quadratic eqn exists
    if (mag(a) < VSMALL)
    {
        return Roots<2>(linearEqn(b, c).roots(), roots::nan, 0);
    }

    // (JLM:p. 2246) [discriminant = b*b/4 - a*c]
    const scalar w = a*c;
    const scalar numDiscr = fma(-a, c, w) + fma(b, b/4, -w);
    const scalar discr = (mag(numDiscr) > VSMALL) ? numDiscr : 0;

    // Find how many roots of what types are available
    const bool twoReal = discr > 0;
    const bool twoComplex = discr < 0;
    //const bool oneReal = discr == 0;

    if (twoReal)
    {
        // (F:Exp. 8.9)
        const scalar x = -b/2 - sign(b)*sqrt(discr);
        return Roots<2>(linearEqn(-a, x).roots(), linearEqn(-x, c).roots());
    }
    else if (twoComplex)
    {
        const Roots<1> xRe(roots::type::complex, -b/2/a);
        const Roots<1> xIm(roots::type::complex, sign(b)*sqrt(mag(discr))/a);
        return Roots<2>(xRe, xIm);
    }
    else // (oneReal)
    {
        const Roots<1> r(linearEqn(a, b/2).roots());
        return Roots<2>(r, r);
    }
}

// ************************************************************************* //
