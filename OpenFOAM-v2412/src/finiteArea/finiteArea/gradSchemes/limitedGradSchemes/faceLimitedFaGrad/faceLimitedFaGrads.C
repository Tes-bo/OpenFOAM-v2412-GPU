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

#include "faceLimitedFaGrad.H"
#include "gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFaGradScheme(faceLimitedGrad)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline void faceLimitedGrad<scalar>::limitEdge
(
    scalar& limiter,
    const scalar& maxDelta,
    const scalar& minDelta,
    const scalar& extrapolate
)
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


template<class Type>
inline void faceLimitedGrad<Type>::limitEdge
(
    Type& limiter,
    const Type& maxDelta,
    const Type& minDelta,
    const Type& extrapolate
)
{
    for (direction cmpt = 0; cmpt < Type::nComponents; ++cmpt)
    {
        faceLimitedGrad<scalar>::limitEdge
        (
            limiter.component(cmpt),
            maxDelta.component(cmpt),
            minDelta.component(cmpt),
            extrapolate.component(cmpt)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<areaVectorField> faceLimitedGrad<scalar>::calcGrad
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

    scalarField maxVsf(vsf.internalField());
    scalarField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const scalar vsfOwn = vsf[own];
        const scalar vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    // Lambda expression to update limiter for boundary edges
    auto updateLimiter = [&](const label patchi, const scalarField& fld) -> void
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();

        forAll(pOwner, facei)
        {
            const label own = pOwner[facei];
            const scalar vsf = fld[facei];

            maxVsf[own] = max(maxVsf[own], vsf);
            minVsf[own] = min(minVsf[own], vsf);
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
        else
        {
            updateLimiter(patchi, psf);
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const scalarField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;
    }


    // Create limiter
    scalarField limiter(vsf.internalField().size(), 1.0);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitEdge
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitEdge
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitEdge
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                (pCf[pFacei] - C[own]) & g[own]
            );
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
tmp<areaTensorField> faceLimitedGrad<vector>::calcGrad
(
    const areaVectorField& vsf,
    const word& name
) const
{
    const faMesh& mesh = vsf.mesh();

    tmp<areaTensorField> tGrad = basicGradScheme_().grad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    areaTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const areaVectorField& C = mesh.areaCentres();
    const edgeVectorField& Cf = mesh.edgeCentres();

    vectorField maxVsf(vsf.internalField());
    vectorField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector& vsfOwn = vsf[own];
        const vector& vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    // Lambda expression to update limiter for boundary edges
    auto updateLimiter = [&](const label patchi, const vectorField& fld) -> void
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();

        forAll(pOwner, facei)
        {
            const label own = pOwner[facei];
            const vector& vsf = fld[facei];

            maxVsf[own] = max(maxVsf[own], vsf);
            minVsf[own] = min(minVsf[own], vsf);
        }
    };

    const areaVectorField::Boundary& bsf = vsf.boundaryField();
    forAll(bsf, patchi)
    {
        const faPatchVectorField& psf = bsf[patchi];

        if (psf.coupled())
        {
            updateLimiter(patchi, psf.patchNeighbourField());
        }
        else
        {
            updateLimiter(patchi, psf);
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const vectorField maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;

        //maxVsf *= 1.0/k_;
        //minVsf *= 1.0/k_;
    }


    // Create limiter
    vectorField limiter(vsf.internalField().size(), vector::one);

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitEdge
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitEdge
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].edgeFaces();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitEdge
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
    }

    if (fa::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    tensorField& gIf = g.primitiveFieldRef();

    forAll(gIf, celli)
    {
        gIf[celli] = tensor
        (
            cmptMultiply(limiter[celli], gIf[celli].x()),
            cmptMultiply(limiter[celli], gIf[celli].y()),
            cmptMultiply(limiter[celli], gIf[celli].z())
        );
    }

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
