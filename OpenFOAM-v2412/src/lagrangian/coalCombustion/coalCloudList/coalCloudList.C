/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
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

#include "coalCloudList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalCloudList::coalCloudList
(
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& slgThermo
)
:
    PtrList<coalCloud>(),
    mesh_(rho.mesh())
{
    IOdictionary props
    (
        IOobject
        (
            "coalCloudList",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ
        )
    );

    const wordHashSet cloudNames(props.get<wordList>("clouds"));

    setSize(cloudNames.size());

    label i = 0;
    for (const word& name : cloudNames)
    {
        Info<< "creating cloud: " << name << endl;

        set
        (
            i++,
            new coalCloud
            (
                name,
                rho,
                U,
                g,
                slgThermo
            )
        );

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalCloudList::evolve()
{
    forAll(*this, i)
    {
        operator[](i).evolve();
    }
}


// ************************************************************************* //
