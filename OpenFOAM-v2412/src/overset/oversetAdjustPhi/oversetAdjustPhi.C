/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "oversetAdjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cellCellStencilObject.H"
#include "syncTools.H"
#include "fv.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::oversetAdjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const label zoneId
)
{
    const fvMesh& mesh = U.mesh();
    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelList& zoneID = overlap.zoneID();

    // Pass1: accumulate all fluxes, calculate correction factor

    scalar massIn = 0;
    scalar massOut = 0;

    surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

    // Check all faces on the outside of interpolated cells
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    forAll(own, facei)
    {
        const label zonei = zoneID[own[facei]];
        const label ownType = cellTypes[own[facei]];
        const label neiType = cellTypes[nei[facei]];

        const bool ownCalc =
            (ownType == cellCellStencil::CALCULATED)
         && (neiType == cellCellStencil::INTERPOLATED);

        const bool neiCalc =
            (ownType == cellCellStencil::INTERPOLATED)
         && (neiType == cellCellStencil::CALCULATED);

        const bool ownOrCalc = (ownCalc || neiCalc);

        if
        (
            (ownOrCalc && (zonei == zoneId))
         || (ownOrCalc && (zoneId == -1))
        )
        {
            // Calculate flux w.r.t. calculated cell
            scalar flux = phi[facei];

            if (ownCalc)
            {
                flux = -flux;
            }

            if (flux < 0.0)
            {
                massIn -= flux;
            }
            else
            {
                massOut += flux;
            }
        }
    }


    // Check all coupled faces on the outside of interpolated cells
    labelList neiCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, neiCellTypes);

    forAll(bphi, patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phip = bphi[patchi];
        const labelUList& fc = Up.patch().faceCells();

        const label start = Up.patch().start();

        forAll(fc, i)
        {
            const label facei = start + i;
            const label celli = fc[i];
            const label ownType = cellTypes[celli];
            const label neiType = neiCellTypes[facei - mesh.nInternalFaces()];

            const label zonei = zoneID[celli];

            const bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            if (ownCalc && (zonei == zoneId))
            {
                // Calculate flux w.r.t. calculated cell
                scalar flux = phip[i];

                if (flux < 0.0)
                {
                    massIn -= flux;
                }
                else
                {
                    massOut += flux;
                }
            }
        }
    }
    reduce(massIn, sumOp<scalar>());
    reduce(massOut, sumOp<scalar>());

    const scalar massCorr = massIn/(massOut + SMALL);

    if (fv::debug)
    {
        Info<< "Zone                    : " << zoneId << nl
            << "mass outflow            : " << massOut << nl
            << "mass inflow             : " << massIn << nl
            << "correction factor       : " << massCorr << endl;
    }


    // Pass2: adjust fluxes
    forAll(own, facei)
    {
        const label zonei = zoneID[own[facei]];  // own and nei in same zone

        const label ownType = cellTypes[own[facei]];
        const label neiType = cellTypes[nei[facei]];

        const bool ownCalc =
            (ownType == cellCellStencil::CALCULATED)
         && (neiType == cellCellStencil::INTERPOLATED);

        const bool neiCalc =
            (ownType == cellCellStencil::INTERPOLATED)
         && (neiType == cellCellStencil::CALCULATED);

        const bool ownOrCalc = (ownCalc || neiCalc);

        if
        (
            (ownOrCalc && (zonei == zoneId)) || (ownOrCalc && (zoneId == -1))
        )
        {
            scalar flux = phi[facei];

            if (ownCalc)
            {
                flux = -flux;
            }

            if (flux < 0.0)
            {
                phi[facei] /= Foam::sqrt(massCorr);
            }
            else
            {
                phi[facei] *= Foam::sqrt(massCorr);
            }
        }
    }

    forAll(bphi, patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        fvsPatchScalarField& phip = bphi[patchi];
        const labelUList& fc = Up.patch().faceCells();

        const label start = Up.patch().start();

        forAll(fc, i)
        {
            const label facei = start + i;
            const label celli = fc[i];
            const label zonei = zoneID[celli];   // note:own and nei in same zone
            const label ownType = cellTypes[celli];
            const label neiType = neiCellTypes[facei - mesh.nInternalFaces()];

            const bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            const bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);

            const bool ownOrCalc = (ownCalc || neiCalc);

            if (ownOrCalc && (zonei == zoneId))
            {
                // Calculate flux w.r.t. calculated cell
                scalar flux = phip[i];

                if (neiCalc)
                {
                    flux = -flux;
                }

                if (flux < 0.0)
                {
                    // phip[i] /= Foam::sqrt(massCorr[zonei]);
                }
                else
                {
                    phip[i] *= massCorr;
                }
            }
        }
    }


    // Check correction
    #ifdef FULLDEBUG
    if (fv::debug)
    {
        scalar massOutCheck = 0;
        scalar massInCheck = 0;

        forAll(own, facei)
        {
            const label zonei = zoneID[own[facei]];
            const label ownType = cellTypes[own[facei]];
            const label neiType = cellTypes[nei[facei]];

            const bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            const bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);

            const bool ownOrCalc = (ownCalc || neiCalc);

            if
            (
                (ownOrCalc && (zonei == zoneId))||(ownOrCalc && (zoneId == -1))
            )
            {
                scalar flux = phi[facei];

                if (ownCalc)
                {
                    flux = -flux;
                }

                if (flux < 0.0)
                {
                    massInCheck -= flux;
                }
                else
                {
                    massOutCheck += flux;
                }
            }
        }

        Info<< "mass in:" << massInCheck << " out:" << massOutCheck << nl;
    }
    #endif  /* FULLDEBUG */

    return true;
}


// ************************************************************************* //
