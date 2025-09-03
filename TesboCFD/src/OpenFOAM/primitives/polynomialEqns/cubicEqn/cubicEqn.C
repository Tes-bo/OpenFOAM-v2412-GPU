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
#include "cubicEqn.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Roots<3> Foam::cubicEqn::roots() const
{
    const scalar a = this->a();
    const scalar b = this->b();
    const scalar c = this->c();
    const scalar d = this->d();

    #ifdef FULLDEBUG
    Info<< "#DEBUG#" << nl
        << "Coefficients of the characteristic cubic polynomial:" << nl
        << "a = " << a << nl
        << "b = " << b << nl
        << "c = " << c << nl
        << "d = " << d << nl
        << "#######" << endl;
    #endif

    // Check the leading term in the cubic eqn exists
    if (mag(a) < VSMALL)
    {
        return Roots<3>(quadraticEqn(b, c, d).roots(), roots::nan, 0);
    }

    // (JLM:p. 2246) [p = a*c - b*b/3]
    const scalar w = a*c;
    const scalar p = -(fma(-a, c, w) + fma(b, b/3.0, -w));
    const scalar q = b*b*b*scalar(2)/27 - b*c*a/3 + d*a*a;
    const scalar numDiscr = p*p*p/27 + q*q/4;
    const scalar discr = (mag(numDiscr) > VSMALL) ? numDiscr : 0;

    // Determine the number and types of the roots
    const bool threeReal = discr < 0;
    const bool oneRealTwoComplex = discr > 0;
    const bool twoReal = //p != 0; && discr == 0;
        (mag(p) > sqrt(SMALL)) && !(threeReal || oneRealTwoComplex);
    // const bool oneReal = p == 0 && discr == 0;

    #ifdef FULLDEBUG
    Info<< "#DEBUG#" << nl
        << "Numerical discriminant:" << tab << numDiscr << nl
        << "Adjusted discriminant:" << tab << discr << nl
        << "Number and types of the roots:" << nl
        << "threeReal = " << threeReal << nl
        << "oneRealTwoComplex = " << oneRealTwoComplex << nl
        << "twoReal = " << twoReal << nl
        << "oneReal = " << !(threeReal || oneRealTwoComplex || twoReal) << nl
        << "#######" << endl;
    #endif

    static const scalar sqrt3 = sqrt(3.0);

    scalar x = 0;

    if (threeReal)
    {
        const scalar wCbRe = -q/2;
        const scalar wCbIm = sqrt(-discr);
        const scalar wAbs = cbrt(hypot(wCbRe, wCbIm));
        const scalar wArg = atan2(wCbIm, wCbRe)/3;
        const scalar wRe = wAbs*cos(wArg);
        const scalar wIm = wAbs*sin(wArg);

        if (b > 0)
        {
            x = -wRe - mag(wIm)*sqrt3 - b/3;
        }
        else
        {
            x = 2*wRe - b/3;
        }
    }
    else if (oneRealTwoComplex)
    {
        const scalar wCb = -q/2 - sign(q)*sqrt(discr);
        const scalar w = cbrt(wCb);
        const scalar t = w - p/(3*w);

        if (p + t*b < 0)
        {
            x = t - b/3;
        }
        else
        {
            const scalar xRe = -t/2 - b/3;
            const scalar xIm = sqrt3/2*(w + p/3/w);

            return
                Roots<3>
                (
                    Roots<1>(roots::real, -a*d/(xRe*xRe + xIm*xIm)),
                    Roots<2>
                    (
                        Roots<1>(roots::complex, xRe),
                        Roots<1>(roots::complex, xIm)
                    )
                );
        }
    }
    else if (twoReal)
    {
        if (q*b > 0)
        {
            x = -2*cbrt(q/2) - b/3;
        }
        else
        {
            x = cbrt(q/2) - b/3;
            const Roots<1> r(linearEqn(-a, x).roots());
            return Roots<3>(Roots<2>(r, r), linearEqn(x*x, a*d).roots());
        }
    }
    else // (oneReal)
    {
        const Roots<1> r(linearEqn(a, b/3).roots());
        return Roots<3>(r.type(0), r[0]);
    }

    #if FULLDEBUG
    Info<< "#DEBUG#" << nl
        << "x = " << x << nl
        << "#######" << endl;
    #endif

    return
        Roots<3>
        (
            linearEqn(-a, x).roots(),
            quadraticEqn(-x*x, c*x + a*d, d*x).roots()
        );
}


// ************************************************************************* //
