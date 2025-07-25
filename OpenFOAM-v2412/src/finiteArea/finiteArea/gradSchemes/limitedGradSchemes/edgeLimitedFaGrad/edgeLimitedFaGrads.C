/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "edgeLimitedFaGrad.H"
#include "gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFaGradScheme(edgeLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void edgeLimitedGrad<Type>::limitEdge
(
    scalar& limiter,
    const scalar maxDelta,
    const scalar minDelta,
    const scalar extrapolate
) const
{
    if (extrapolate > maxDelta + VSMALL)
    {
        limiter = min(limiter, maxDelta/extrapolate);
    }
    else if (extrapolate < minDelta - VSMALL)
    {
        limiter = min(limiter, minDelta/extrapolate);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<areaVectorField> edgeLimitedGrad<scalar>::calcGrad
(
    const areaScalarField& vsf,
    const word& name
) const
{
    const faMesh& mesh = vsf.mesh();

    tmp<areaVectorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaVectorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // Create limiter field
    scalarField limiter(vsf.internalField().size(), 1.0);

    const scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        const label own = owner[edgei];
        const label nei = neighbour[edgei];

        const scalar vsfOwn = vsf[own];
        const scalar vsfNei = vsf[nei];

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        const scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        // owner side
        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn,
            minEdge - vsfOwn,
            (Cf[edgei] - C[own]) & g[own]
        );

        // neighbour side
        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei,
            minEdge - vsfNei,
            (Cf[edgei] - C[nei]) & g[nei]
        );
    }

    // Lambda expression to update limiter for boundary edges
    auto updateLimiter = [&](const label patchi, const scalarField& fld) -> void
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, edgei)
        {
            const label own = pOwner[edgei];

            const scalar vsfOwn = vsf[own];
            const scalar vsfNei = fld[edgei];

            scalar maxEdge = max(vsfOwn, vsfNei);
            scalar minEdge = min(vsfOwn, vsfNei);
            const scalar maxMinEdge = rk*(maxEdge - minEdge);
            maxEdge += maxMinEdge;
            minEdge -= maxMinEdge;

            limitEdge
            (
                limiter[own],
                maxEdge - vsfOwn,
                minEdge - vsfOwn,
                (pCf[edgei] - C[own]) & g[own]
            );
        }
    };

    const areaScalarField::Boundary& bsf = vsf.boundaryField();
    forAll(bsf, patchi)
    {
        const faPatchScalarField& psf = bsf[patchi];

        if (psf.coupled())
        {
            updateLimiter(patchi, psf.patchNeighbourField());
        }
        else if (psf.fixesValue())
        {
            updateLimiter(patchi, psf);
        }
    }

    if (fa::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.primitiveFieldRef() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


template<>
tmp<areaTensorField> edgeLimitedGrad<vector>::calcGrad
(
    const areaVectorField& vvf,
    const word& name
) const
{
    const faMesh& mesh = vvf.mesh();

    tmp<areaTensorField> tGrad = basicGradScheme_().grad(vvf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    // Create limiter
    scalarField limiter(vvf.internalField().size(), 1.0);

    const scalar rk = (1.0/k_ - 1.0);

    forAll(owner, edgei)
    {
        const label own = owner[edgei];
        const label nei = neighbour[edgei];

        // owner side
        vector gradf((Cf[edgei] - C[own]) & g[own]);

        scalar vsfOwn = gradf & vvf[own];
        scalar vsfNei = gradf & vvf[nei];

        scalar maxEdge = max(vsfOwn, vsfNei);
        scalar minEdge = min(vsfOwn, vsfNei);
        const scalar maxMinEdge = rk*(maxEdge - minEdge);
        maxEdge += maxMinEdge;
        minEdge -= maxMinEdge;

        limitEdge
        (
            limiter[own],
            maxEdge - vsfOwn,
            minEdge - vsfOwn,
            magSqr(gradf)
        );


        // neighbour side
        gradf = (Cf[edgei] - C[nei]) & g[nei];

        vsfOwn = gradf & vvf[own];
        vsfNei = gradf & vvf[nei];

        maxEdge = max(vsfOwn, vsfNei);
        minEdge = min(vsfOwn, vsfNei);

        limitEdge
        (
            limiter[nei],
            maxEdge - vsfNei,
            minEdge - vsfNei,
            magSqr(gradf)
        );
    }


    // Lambda expression to update limiter for boundary edges
    auto updateLimiter = [&](const label patchi, const vectorField& fld) -> void
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, edgei)
        {
            const label own = pOwner[edgei];

            const vector gradf((pCf[edgei] - C[own]) & g[own]);

            const scalar vsfOwn = gradf & vvf[own];
            const scalar vsfNei = gradf & fld[edgei];

            scalar maxEdge = max(vsfOwn, vsfNei);
            scalar minEdge = min(vsfOwn, vsfNei);
            const scalar maxMinEdge = rk*(maxEdge - minEdge);
            maxEdge += maxMinEdge;
            minEdge -= maxMinEdge;

            limitEdge
            (
                limiter[own],
                maxEdge - vsfOwn,
                minEdge - vsfOwn,
                magSqr(gradf)
            );
        }
    };


    const areaVectorField::Boundary& bvf = vvf.boundaryField();
    forAll(bvf, patchi)
    {
        const faPatchVectorField& psf = bvf[patchi];

        if (psf.coupled())
        {
            updateLimiter(patchi, psf.patchNeighbourField());
        }
        else if (psf.fixesValue())
        {
            updateLimiter(patchi, psf);
        }
    }

    if (fa::debug)
    {
        Info<< "gradient limiter for: " << vvf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.primitiveFieldRef() *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vvf, g);

    return tGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
