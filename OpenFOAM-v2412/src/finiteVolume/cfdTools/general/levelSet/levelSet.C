/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "levelSet.H"
#include "cut.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::levelSetFraction
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const bool above
)
{
    auto tResult = DimensionedField<scalar, volMesh>::New
    (
        "levelSetFraction",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    auto& result = tResult.ref();

    forAll(result, cI)
    {
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cI);

        scalar v = 0, r = 0;

        forAll(cellTetIs, cellTetI)
        {
            const triFace triIs = cellTetIs[cellTetI].faceTriIs(mesh);

            const FixedList<point, 4>
                tet =
                {
                    mesh.cellCentres()[cI],
                    mesh.points()[triIs[0]],
                    mesh.points()[triIs[1]],
                    mesh.points()[triIs[2]]
                };
            const FixedList<scalar, 4>
                level =
                {
                    levelC[cI],
                    levelP[triIs[0]],
                    levelP[triIs[1]],
                    levelP[triIs[2]]
                };

            v += cut::volumeOp()(tet);

            if (above)
            {
                r += tetCut(tet, level, cut::volumeOp(), cut::noOp());
            }
            else
            {
                r += tetCut(tet, level, cut::noOp(), cut::volumeOp());
            }
        }

        result[cI] = r/v;
    }

    return tResult;
}


Foam::tmp<Foam::scalarField> Foam::levelSetFraction
(
    const fvPatch& patch,
    const scalarField& levelF,
    const scalarField& levelP,
    const bool above
)
{
    auto tResult = tmp<scalarField>::New(patch.size(), Zero);
    auto& result = tResult.ref();

    forAll(result, fI)
    {
        const face& f = patch.patch().localFaces()[fI];

        vector a(Zero);
        vector r(Zero);

        for (label edgei = 0; edgei < f.nEdges(); ++edgei)
        {
            const edge e = f.edge(edgei);

            const FixedList<point, 3>
                tri =
                {
                    patch.patch().faceCentres()[fI],
                    patch.patch().localPoints()[e[0]],
                    patch.patch().localPoints()[e[1]]
                };
            const FixedList<scalar, 3>
                level =
                {
                    levelF[fI],
                    levelP[e[0]],
                    levelP[e[1]]
                };

            a += cut::areaOp()(tri);

            if (above)
            {
                r += triCut(tri, level, cut::areaOp(), cut::noOp());
            }
            else
            {
                r += triCut(tri, level, cut::noOp(), cut::areaOp());
            }
        }

        result[fI] = a/magSqr(a) & r;
    }

    return tResult;
}

// ************************************************************************* //
