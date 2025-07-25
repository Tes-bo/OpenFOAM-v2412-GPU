/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

Application
    Test-Polynomial

Description
    Test application for the templated Polynomial class

\*---------------------------------------------------------------------------*/

#include "Polynomial.H"
#include "FixedList.H"
#include "polynomialFunction.H"
#include "ITstream.H"
#include "OTstream.H"
#include "Random.H"
#include "cpuTime.H"

using namespace Foam;

std::initializer_list<scalar> coeffs
{
    0.11, 0.45, -0.94, 1.58, -2.58, 0.08, 3.15, -4.78
};

const FixedList<scalar, 8> coeff(coeffs);


scalar polyValue(const scalar x)
{
    // Hard-coded polynomial 8 coeff (7th order)
    return
        coeff[0]
      + coeff[1]*x
      + coeff[2]*sqr(x)
      + coeff[3]*pow3(x)
      + coeff[4]*pow4(x)
      + coeff[5]*pow5(x)
      + coeff[6]*pow6(x)
      + coeff[7]*x*pow6(x);
}


scalar intPolyValue(const scalar x)
{
    // Hard-coded integral form of above polynomial
    return
        coeff[0]*x
      + coeff[1]/2.0*sqr(x)
      + coeff[2]/3.0*pow3(x)
      + coeff[3]/4.0*pow4(x)
      + coeff[4]/5.0*pow5(x)
      + coeff[5]/6.0*pow6(x)
      + coeff[6]/7.0*x*pow6(x)
      + coeff[7]/8.0*x*x*pow6(x);
}


scalar polyValue1(const scalar x)
{
    // Naive evaluation using pow()
    scalar value = coeff[0];

    for (label i=1; i < coeff.size(); ++i)
    {
        value += coeff[i]*pow(x, i);
    }

    return value;
}


scalar polyValue2(const scalar x)
{
    // Calculation avoiding pow()
    scalar value = coeff[0];

    scalar powX = x;
    for (label i=1; i < coeff.size(); ++i)
    {
        value += coeff[i] * powX;
        powX *= x;
    }

    return value;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    constexpr label n = 10000;
    constexpr label nIters = 1000;
    scalar sum = 0;

    Info<< "null poly = " << (Polynomial<8>{}) << nl
        << "null poly = " << (polynomialFunction{8}) << nl
        << endl;

    Polynomial<8> poly{coeffs};
    Polynomial<9> intPoly{poly.integral(0)};

    polynomialFunction polyfunc;

    // Could profit from a bi-directional stream
    {
        OTstream os;
        os << poly;

        ITstream is(std::move(os.tokens()));
        is >> polyfunc;
    }

    Info<< "poly = " << poly << nl
        << "intPoly = " << intPoly << nl
        << endl;

    Info<< "2*poly = " << 2*poly << nl
        << "poly+poly = " << poly + poly << nl
        << "3*poly = " << 3*poly << nl
        << "poly+poly+poly = " << poly + poly + poly << nl
        << "3*poly - 2*poly = " << 3*poly - 2*poly << nl
        << endl;

    Info<< nl << "--- as polynomialFunction" << nl << endl;
    Info<< "polyf = " << polyfunc << nl
        << "intPoly = " << poly.integral(0.0) << nl
        << endl;

    Info<< "2*polyf = " << 2*polyfunc << nl
        << "polyf+polyf = " << polyfunc + polyfunc << nl
        << "3*polyf = " << 3*polyfunc << nl
        << "polyf+polyf+polyf = " << polyfunc + polyfunc + polyfunc << nl
        << "3*polyf - 2*polyf = " << 3*polyfunc - 2*polyfunc << nl
        << endl;


    Polynomial<8> polyCopy = poly;
    Info<< "poly, polyCopy = " << poly << ", " << polyCopy << nl << endl;
    polyCopy = 2.5*poly;
    Info<< "2.5*poly = " << polyCopy << nl << endl;

    Random rnd(123456);
    for (label i=0; i<10; ++i)
    {
        scalar x = rnd.sample01<scalar>()*100;

        scalar px = polyValue(x);
        scalar ipx = intPolyValue(x);

        scalar pxTest = poly.value(x);
        scalar ipxTest = intPoly.value(x);

        Info<<"\nx = " << x << nl
            << "    px, pxTest = " << px << ", " << pxTest << nl
            << "    ipx, ipxTest = " << ipx << ", " << ipxTest << nl;

        if (mag(px - pxTest) > SMALL)
        {
            Info<< "    *** WARNING: px != pxTest: " << px - pxTest << nl;
        }

        if (mag(ipx - ipxTest) > SMALL)
        {
            Info<< "    *** WARNING: ipx != ipxTest: " << ipx - ipxTest << nl;
        }

        Info<< endl;
    }


    //
    // test speed of Polynomial:
    //
    Info<< "start timing loops" << nl
        << "~~~~~~~~~~~~~~~~~~" << endl;

    cpuTime timer;

    for (label loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += poly.value(loop+iter);
        }
    }
    Info<< "value:        " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";

    for (label loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyfunc.value(loop+iter);
        }
    }
    Info<< "via function: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    for (label loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue(loop+iter);
        }
    }
    Info<< "hard-coded 0: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    for (label loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue1(loop+iter);
        }
    }
    Info<< "hard-coded 1: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";

    for (label loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue2(loop+iter);
        }
    }
    Info<< "hard-coded 2: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    Info<< nl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
