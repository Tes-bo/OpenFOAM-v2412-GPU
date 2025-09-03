/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "viewFactor2AI.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(viewFactor2AI, 0);
    addToRunTimeSelectionTable(viewFactorModel, viewFactor2AI, mesh);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::VF::viewFactor2AI::calculateFij
(
    const point& xi,
    const point& xj,
    const vector& dAi,
    const vector& dAj
)
{
    const vector r(xj - xi);
    const scalar rMag = mag(r);
    const scalar dAiMag(mag(dAi));
    const scalar dAjMag(mag(dAj));

    if (rMag > SMALL && dAiMag > SMALL && dAjMag > SMALL)
    {
        const vector nr(r/rMag);
        const vector ni(dAi/dAiMag);
        const vector nj(dAj/dAjMag);

        const scalar Fij =
            -(nr & ni)*(nr & nj)/sqr(rMag)*dAjMag*dAiMag/mathematical::pi;

        if (Fij > 0) return Fij;
    }

    return 0;
}


Foam::scalarListList Foam::VF::viewFactor2AI::calculate
(
    const labelListList& visibleFaceFaces,
    const pointField& compactCf,
    const vectorField& compactSf,
    const List<List<vector>>& compactFineSf,
    const List<List<point>>& compactFineCf,
    const DynamicList<List<point>>& compactPoints,
    const DynamicList<label>& compactPatchId
) const
{
    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarListList Fij(visibleFaceFaces.size());

    forAll(visibleFaceFaces, facei)
    {
        if (debug > 1)
        {
            Pout<< "facei:" << facei << "/" << visibleFaceFaces.size()
                << endl;
        }

        const labelList& visibleFaces = visibleFaceFaces[facei];

        Fij[facei].resize(visibleFaces.size(), Zero);

        const point& dCi = compactCf[facei];
        const vector& Ai = compactSf[facei];
        const scalar magAi = mag(Ai);

        if (magAi < SMALL) continue;

        forAll(visibleFaces, visibleFacei)
        {
            const label sloti = visibleFaces[visibleFacei];
            const point& dCj = compactCf[sloti];
            const vector& Aj = compactSf[sloti];

            const scalar Fij2AI = calculateFij(dCi, dCj, Ai, Aj);

            Fij[facei][visibleFacei] = Fij2AI/magAi;
        }
    }


    return Fij;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::viewFactor2AI::viewFactor2AI
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    viewFactorModel(mesh, dict)
{}


// ************************************************************************* //
