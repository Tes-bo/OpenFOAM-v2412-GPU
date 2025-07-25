/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016 OpenCFD Ltd.
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

Namespace
    Foam::Elliptic

Description

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Elliptic
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void ellipticIntegralsKE(const scalar m, scalar& K, scalar& E)
{
    if (m == 0)
    {
        K = 0.5*mathematical::pi;
        E = 0.5*mathematical::pi;
        return;
    }

    scalar a = 1;
    scalar g = Foam::sqrt(1 - m);
    scalar ga = g*a;
    scalar aux = 1;
    scalar sum = 2 - m;

    while (true)
    {
        scalar aOld = a;
        scalar gOld = g;
        ga = gOld*aOld;
        a = 0.5*(gOld + aOld);
        aux += aux;
        sum -= aux*(a*a - ga);

        if (mag(aOld - gOld) < SMALL)
        {
            break;
        }

        g = sqrt(ga);
    }

    K = 0.5*mathematical::pi/a;
    E = 0.25*mathematical::pi/a*sum;
}


Foam::scalar JacobiAmp(const scalar u, const scalar mIn)
{
    static const label ITER = 25;
    FixedList<scalar, ITER + 1> a, g, c;
    scalar aux, amp;
    label n;

    const scalar m = mag(mIn);

    if (m == 0)
    {
        return u;
    }

    if (m == 1)
    {
        return 2*atan(exp(u)) - mathematical::pi/2;
    }

    a[0] = 1.0;
    g[0] = Foam::sqrt(1.0 - m);
    c[0] = Foam::sqrt(m);

    aux = 1.0;

    for (n = 0; n < ITER; n++)
    {
        if (mag(a[n] - g[n]) < SMALL)
        {
            break;
        }

        aux += aux;
        a[n+1] = 0.5*(a[n] + g[n]);
        g[n+1] = Foam::sqrt(a[n]*g[n]);
        c[n+1] = 0.5*(a[n] - g[n]);
    }

    amp = aux*a[n]*u;

    for (; n > 0; n--)
    {
        amp = 0.5*(amp + asin(c[n]*sin(amp)/a[n]));
    }

    return scalar(amp);
}


void JacobiSnCnDn
(
    const scalar u,
    const scalar m,
    scalar& Sn,
    scalar& Cn,
    scalar& Dn
)
{
    const scalar amp = Elliptic::JacobiAmp(u, m);

    Sn = sin(amp);
    Cn = cos(amp);
    Dn = sqrt(1.0 - m*sin(amp)*sin(amp));
    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Elliptic
} // End namespace Foam

// ************************************************************************* //
