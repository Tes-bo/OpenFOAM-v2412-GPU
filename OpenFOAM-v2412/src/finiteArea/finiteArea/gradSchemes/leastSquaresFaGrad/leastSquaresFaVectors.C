/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "leastSquaresFaVectors.H"
#include "edgeFields.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresFaVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresFaVectors::leastSquaresFaVectors(const faMesh& mesh)
:
    MeshObject_type(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresFaVectors::makeLeastSquaresVectors() const
{
    DebugInFunction
        << "Constructing finite area (invDist) least square gradient vectors"
        << nl;

    pVectorsPtr_ = std::make_unique<edgeVectorField>
    (
        IOobject
        (
            "LeastSquaresP",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedVector(dimless/dimLength, Zero)
    );
    auto& lsP = *pVectorsPtr_;

    nVectorsPtr_ = std::make_unique<edgeVectorField>
    (
        IOobject
        (
            "LeastSquaresN",
            mesh().pointsInstance(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh(),
        dimensionedVector(dimless/dimLength, Zero)
    );
    auto& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const areaVectorField& C = mesh().areaCentres();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nFaces(), Zero);

    // No contribution when mag(val) < SMALL
    const scalar minLenSqr(SMALL*SMALL);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);

        // No contribution when mag(val) < SMALL
        const scalar magSqrDist = d.magSqr();
        if (magSqrDist >= minLenSqr)
        {
            const symmTensor wdd(sqr(d)/magSqrDist);
            dd[own] += wdd;
            dd[nei] += wdd;
        }
    }


    for (const auto& patchLsP : lsP.boundaryField())
    {
        const faPatch& p = patchLsP.patch();
        const labelUList& edgeFaces = p.edgeFaces();

        // Build the d-vectors
        const vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];

            // No contribution when mag(val) < SMALL
            const scalar magSqrDist = d.magSqr();
            if (magSqrDist >= minLenSqr)
            {
                dd[edgeFaces[patchFacei]] += sqr(d)/magSqrDist;
            }
        }
    }


    // Invert the dd tensors - including failsafe checks
    const symmTensorField invDd(inv(dd));


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector d(C[nei] - C[own]);

        // No contribution when mag(val) < SMALL
        const scalar magSqrDist = d.magSqr();
        if (magSqrDist >= minLenSqr)
        {
            lsP[facei] = (invDd[own] & d)/magSqrDist;
            lsN[facei] = -(invDd[nei] & d)/magSqrDist;
        }
        // ... already zero from initialisation
        // else
        // {
        //     lsP[facei] = Zero;
        //     lsN[facei] = Zero;
        // }
    }

    for (auto& patchLsP : lsP.boundaryFieldRef())
    {
        const faPatch& p = patchLsP.patch();
        const labelUList& edgeFaces = p.edgeFaces();

        // Build the d-vectors
        const vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];

            // No contribution when mag(val) < SMALL
            const scalar magSqrDist = d.magSqr();
            if (magSqrDist >= minLenSqr)
            {
                patchLsP[patchFacei] =
                    (invDd[edgeFaces[patchFacei]] & d)/magSqrDist;
            }
            // ... already zero from initialisation
            // else
            // {
            //     patchLsP[patchFacei] = Zero;
            // }
        }
    }

    DebugInfo
        << "Done constructing finite area least square gradient vectors" << nl;
}


const Foam::edgeVectorField& Foam::leastSquaresFaVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::edgeVectorField& Foam::leastSquaresFaVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresFaVectors::movePoints()
{
    DebugInFunction
        << "Clearing least square data" << nl;

    pVectorsPtr_.reset(nullptr);
    nVectorsPtr_.reset(nullptr);

    return true;
}


// ************************************************************************* //
