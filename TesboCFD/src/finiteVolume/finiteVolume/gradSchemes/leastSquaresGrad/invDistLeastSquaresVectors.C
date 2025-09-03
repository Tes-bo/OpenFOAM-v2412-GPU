/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "leastSquaresVectors.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, Zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, Zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectors::calcLeastSquaresVectors()
{
    DebugInFunction << "Calculating least square gradient vectors" << nl;

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), Zero);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);
        const scalar magSqrDist = d.magSqr();

        if (magSqrDist > ROOTVSMALL)
        {
            const symmTensor wdd(sqr(d)/magSqrDist);
            dd[own] += wdd;
            dd[nei] += wdd;
        }
    }

    auto& blsP = pVectors_.boundaryField();

    forAll(blsP, patchi)
    {
        const fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvPatch& p = patchLsP.patch();
        const labelUList& faceCells = p.patch().faceCells();

        // Build the d-vectors
        const vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];
            const scalar magSqrDist = d.magSqr();

            if (magSqrDist > ROOTVSMALL)
            {
                dd[faceCells[patchFacei]] += sqr(d)/magSqrDist;
            }
        }
    }


    // Invert the dd tensors - including failsafe checks
    const symmTensorField invDd(inv(dd));


    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);
        const scalar magSqrDist = d.magSqr();

        if (magSqrDist > ROOTVSMALL)
        {
            pVectors_[facei] = (invDd[own] & d)/magSqrDist;
            nVectors_[facei] = -(invDd[nei] & d)/magSqrDist;
        }
        else
        {
            pVectors_[facei] = Zero;
            nVectors_[facei] = Zero;
        }
    }

    forAll(blsP, patchi)
    {
        fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvPatch& p = patchLsP.patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        const vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];
            const scalar magSqrDist = d.magSqr();

            if (magSqrDist > ROOTVSMALL)
            {
                patchLsP[patchFacei] =
                    (invDd[faceCells[patchFacei]] & d)/magSqrDist;
            }
            else
            {
                patchLsP[patchFacei] = Zero;
            }
        }
    }

    DebugInfo << "Finished calculating least square gradient vectors" << endl;
}


bool Foam::leastSquaresVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
