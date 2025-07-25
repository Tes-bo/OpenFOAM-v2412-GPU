/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "equipotentialPointSmoother.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointSmoothers
{
    defineTypeNameAndDebug(equipotentialPointSmoother, 0);
    addToRunTimeSelectionTable
    (
        pointSmoother,
        equipotentialPointSmoother,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSmoothers::equipotentialPointSmoother::equipotentialPointSmoother
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointSmoother(mesh, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointSmoothers::equipotentialPointSmoother::calculate
(
    const labelList& facesToMove,
    const pointField& oldPoints,
    const pointField& currentPoints,
    const pointField& faceCentres,
    const vectorField& faceAreas,
    const pointField& cellCentres,
    const scalarField& cellVolumes,
    vectorField& pointDisplacement
) const
{
    // Number of points used in each average
    scalarField weights(mesh().nPoints(), 0);

    // Reset the displacements which are about to be calculated
    reset(facesToMove, weights, pointDisplacement);

    // Sum the non-internal face displacements
    forAll(facesToMove, faceToMoveI)
    {
        const label faceI(facesToMove[faceToMoveI]);

        if (!isInternalOrProcessorFace(faceI))
        {
            const face& fPoints(mesh().faces()[faceI]);

            const scalar area(mag(mesh().faceAreas()[faceI]));

            forAll(fPoints, fPointI)
            {
                const label pointI(fPoints[fPointI]);

                pointDisplacement[pointI] +=
                    area
                   *(
                        faceCentres[faceI]
                      - oldPoints[pointI]
                    );

                weights[pointI] += area;
            }
        }
    }

    // Sum the internal face displacements
    forAll(facesToMove, faceToMoveI)
    {
        const label faceI(facesToMove[faceToMoveI]);

        if (isInternalOrProcessorFace(faceI))
        {
            const face& fPoints(mesh().faces()[faceI]);

            forAll(fPoints, fPointI)
            {
                const label pointI(fPoints[fPointI]);

                if (weights[pointI] < SMALL)
                {
                    const labelList& pCells(mesh().pointCells()[pointI]);

                    forAll(pCells, pCellI)
                    {
                        const label cellI(pCells[pCellI]);

                        const scalar volume(mesh().cellVolumes()[cellI]);

                        pointDisplacement[pointI] +=
                            volume
                           *(
                                cellCentres[cellI]
                              - oldPoints[pointI]
                            );

                        weights[pointI] += volume;
                    }
                }
            }
        }
    }

    // Average
    average(facesToMove, weights, pointDisplacement);
}


// ************************************************************************* //
