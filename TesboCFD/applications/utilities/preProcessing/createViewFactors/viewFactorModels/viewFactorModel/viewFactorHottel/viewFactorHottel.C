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

#include "viewFactorHottel.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"
#include "meshTools.H"
//#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(viewFactorHottel, 0);
//    addToRunTimeSelectionTable(viewFactorModel, viewFactorHottel, mesh);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::VF::viewFactorHottel::calculateFij
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    return 0.5*(mag(p2-p1) + mag(p3-p0) - mag(p2-p0) - mag(p3-p1));
}


Foam::scalarListList Foam::VF::viewFactorHottel::calculate
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

        Fij[facei].resize_nocopy(visibleFaces.size());

        const point& dCi = compactCf[facei];
        const vector& Ai = compactSf[facei];
        const scalar magAi = mag(Ai);

        const vector d1((Ai/magAi) ^ emptyDir_);
        const vector l1(0.5*magAi/w_*d1);
        const point p0(dCi + l1);
        const point p1(dCi - l1);

        forAll(visibleFaces, visibleFacei)
        {
            const label sloti = visibleFaces[visibleFacei];

            const point& dCj = compactCf[sloti];
            const vector& Aj = compactSf[sloti];
            const scalar magAj = mag(Aj);

            const vector d2((Aj/magAj) ^ emptyDir_);
            const vector l2(0.5*magAj/w_*d2);
            const point p2(dCj - l2);
            const point p3(dCj + l2);

            const scalar FijH = calculateFij(p0, p1, p2, p3);

            Fij[facei][visibleFacei] = FijH/(magAi/w_);
        }
    }


    return Fij;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::viewFactorHottel::viewFactorHottel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    viewFactorModel(mesh, dict),
    emptyDir_(vector::one),
    w_(0)
{
    if (mesh.nSolutionD() != 2)
    {
        FatalErrorInFunction
            << "Hottel crossed strings method only applicable to 2D cases"
            << exit(FatalError);
    }

    meshTools::constrainDirection(mesh, mesh.solutionD(), emptyDir_);
    emptyDir_ = vector::one - emptyDir_;
    emptyDir_.normalise();

    // 2D width - assume slab
    // TODO: warn wedge/axisymmetric?
    w_ = mesh.bounds().span() & emptyDir_;

    Info<< "\nEmpty direction: " << emptyDir_
        << "\nWidth: " << w_ << endl;
}


// ************************************************************************* //
