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

#include "viewFactor2LI.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(viewFactor2LI, 0);
    addToRunTimeSelectionTable(viewFactorModel, viewFactor2LI, mesh);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::VF::viewFactor2LI::calculateFij
(
    const List<point>& lPoints,
    const List<point>& rPoints,
    const scalar alpha
)
{
    scalar Fij = 0;

    forAll(lPoints, i)
    {
        // Edge vector and centroid of edge i
        const vector si(lPoints[i] - lPoints.rcValue(i));
        const point ci(0.5*(lPoints[i] + lPoints.rcValue(i)));

        forAll(rPoints, j)
        {
            // Edge vector and centroid of edge j
            const vector sj(rPoints[j] - rPoints.rcValue(j));
            const point cj(0.5*(rPoints[j] + rPoints.rcValue(j)));

            vector r(ci - cj);
            if (mag(r) < SMALL)
            {
                r = alpha*si;
            }

            Fij += (si & sj)*Foam::log(r & r);
        }
    }

    return max(0, 0.25*Fij/mathematical::pi);
}


Foam::scalarListList Foam::VF::viewFactor2LI::calculate
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

        const vector& Ai = compactSf[facei];
        const scalar magAi = mag(Ai);

        forAll(visibleFaces, visibleFacei)
        {
            const label sloti = visibleFaces[visibleFacei];
            const List<point>& lPoints = compactPoints[facei];
            const List<point>& rPoints = compactPoints[sloti];

            const scalar Fij2LI = calculateFij(lPoints, rPoints, alpha_);

            Fij[facei][visibleFacei] = Fij2LI/magAi;
        }
    }

    return Fij;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::viewFactor2LI::viewFactor2LI
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    viewFactorModel(mesh, dict),
    alpha_(dict.getOrDefault("alpha", 0.21))
{}


// ************************************************************************* //
