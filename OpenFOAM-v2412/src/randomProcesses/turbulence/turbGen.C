/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "fft.H"
#include "turbGen.H"
#include "Kmesh.H"
#include "primitiveFields.H"
#include "Ek.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbGen::turbGen(const Kmesh& k, const scalar EA, const scalar K0)
:
    K(k),
    Ea(EA),
    k0(K0),
    RanGen(label(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::turbGen::U()
{
    vectorField s(K.size());
    scalarField rndPhases(K.size());

    forAll(K, i)
    {
        s[i] = RanGen.sample01<vector>();
        rndPhases[i] = RanGen.sample01<scalar>();
    }

    s = K ^ s;
    s = s/(mag(s) + 1.0e-20);

    s = Ek(Ea, k0, mag(K))*s;

    label ntot = 1;
    forAll(K.nn(), idim)
    {
        ntot *= K.nn()[idim];
    }
    const scalar recRootN = 1.0/sqrt(scalar(ntot));

    complexVectorField up
    (
        fft::reverseTransform
        (
            ComplexField
            (
                cos(constant::mathematical::twoPi*rndPhases)*s,
                sin(constant::mathematical::twoPi*rndPhases)*s
            ),
            K.nn()
        )*recRootN
    );

    return ReImSum(up);
}


// ************************************************************************* //
