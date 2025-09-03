/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "error.H"

#include "UOprocess.H"
#include "Kmesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

complexVector UOprocess::WeinerProcess()
{
    return RootDeltaT*complexVector
    (
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>()),
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>()),
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

UOprocess::UOprocess
(
    const Kmesh& kmesh,
    const scalar deltaT,
    const dictionary& UOdict
)
:
    GaussGen(),
    Mesh(kmesh),
    DeltaT(deltaT),
    RootDeltaT(sqrt(DeltaT)),
    UOfield(Mesh.size()),

    Alpha(UOdict.get<scalar>("UOalpha")),
    Sigma(UOdict.get<scalar>("UOsigma")),
    Kupper(UOdict.get<scalar>("UOKupper")),
    Klower(UOdict.get<scalar>("UOKlower")),
    Scale((Kupper - Klower)*pow(scalar(Mesh.size()), 1.0/vector::dim))
{
    const vectorField& K = Mesh;

    scalar sqrKupper = sqr(Kupper);
    scalar sqrKlower = sqr(Klower) + SMALL;
    scalar sqrK;

    forAll(UOfield, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrKupper && sqrK > sqrKlower)
        {
            UOfield[i] = Scale*Sigma*WeinerProcess();
        }
        else
        {
            UOfield[i] = complexVector::zero;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const complexVectorField& UOprocess::newField()
{
    const vectorField& K = Mesh;

    label count = 0;
    scalar sqrKupper = sqr(Kupper);
    scalar sqrKlower = sqr(Klower) + SMALL;
    scalar sqrK;

    forAll(UOfield, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrKupper && sqrK > sqrKlower)
        {
            count++;
            UOfield[i] =
                (1.0 - Alpha*DeltaT)*UOfield[i]
              + Scale*Sigma*WeinerProcess();
        }
    }

    Info<< "    Number of forced K = " << count << nl;

    return UOfield;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
