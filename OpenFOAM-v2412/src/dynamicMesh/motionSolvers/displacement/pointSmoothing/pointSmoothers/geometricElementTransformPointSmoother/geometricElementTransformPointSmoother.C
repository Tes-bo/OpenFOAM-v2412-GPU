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

#include "geometricElementTransformPointSmoother.H"
#include "cellPointConnectivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointSmoothers
{
    defineTypeNameAndDebug(geometricElementTransformPointSmoother, 0);
    addToRunTimeSelectionTable
    (
        pointSmoother,
        geometricElementTransformPointSmoother,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::pointSmoothers::geometricElementTransformPointSmoother::cellQuality
(
    const polyMesh& mesh,
    const pointField& currentPoints,
    const cellPointConnectivity& connectivity,
    const label celli
)
{
    const cell& cFaces = mesh.cells()[celli];
    const labelList cPoints(cFaces.labels(mesh.faces()));

    // Calculate a transformed point for each cell point

    scalar cellQ = 0.0;
    label nTets = 0;

    forAll(cPoints, cPointi)
    {
        const label pointi(cPoints[cPointi]);
        const point& pt = currentPoints[pointi];
        const labelList& pPoints
        (
            connectivity.cellPointPoints()[celli][cPointi]
        );
        if (pPoints.size() != 3)
        {
            WarningInFunction<< "Cell:" << celli
                << " point:" << pointi
                << " connected points:" << pPoints << endl;
        }
        else
        {
            const Tensor<scalar> mA
            (
                currentPoints[pPoints[0]] - pt,
                currentPoints[pPoints[1]] - pt,
                currentPoints[pPoints[2]] - pt
            );
            const scalar sigma(det(mA));

            if (sigma < ROOTVSMALL)
            {
                return 0;
            }

            // 3 * pow(sigma, 2.0/3.0)/magSqr(mA)
            const scalar tetQ =
                scalar(3) * Foam::cbrt(Foam::sqr(sigma)) / magSqr(mA);
            cellQ += tetQ;
            nTets++;
        }
    }

    return cellQ/nTets;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSmoothers::geometricElementTransformPointSmoother::
geometricElementTransformPointSmoother
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointSmoother(mesh, dict),
    transformationParameter_
    (
        readScalar(dict.lookup("transformationParameter"))
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointSmoothers::geometricElementTransformPointSmoother::calculate
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
    // Lookup or generate the cell-point connectivity/
    const cellPointConnectivity& connectivity =
        MeshObject<polyMesh, MoveableMeshObject, cellPointConnectivity>::New
        (
            mesh()
        );

    // Number of points used in each average
    labelField counts(mesh().nPoints(), -1);

    // Reset the displacements which are about to be calculated
    reset(facesToMove, counts, pointDisplacement);

    // Identify the cells which are to be moved
    labelHashSet cellsToMove(facesToMove.size()*2/3);
    forAll(facesToMove, faceToMoveI)
    {
        const label faceI(facesToMove[faceToMoveI]);

        cellsToMove.insert(mesh().faceOwner()[faceI]);

        if (mesh().isInternalFace(faceI))
        {
            cellsToMove.insert(mesh().faceNeighbour()[faceI]);
        }
    }

    // Transformed point field
    pointField transformedPoints(currentPoints);

    // Calculate the internal transformations
    forAllConstIter(labelHashSet, cellsToMove, iter)
    {
        const label cellI(iter.key());

        const cell& cFaces
        (
            mesh().cells()[cellI]
        );
        const labelList cPoints
        (
            cFaces.labels(mesh().faces())
        );
        const edgeList cEdges
        (
            cFaces.edges(mesh().faces())
        );

        // Calculate a transformed point for each cell point
        forAll(cPoints, cPointI)
        {
            const label pointI(cPoints[cPointI]);

            if (counts[pointI] == -1) continue;

            const labelList& pPoints
            (
                connectivity.cellPointPoints()[cellI][cPointI]
            );
            const labelList& pFaces
            (
                connectivity.cellPointFaces()[cellI][cPointI]
            );
            const label nPPoints(pPoints.size());

            // Initial guess of the dual face centre
            vector dualAverage(vector::zero);
            forAll(pPoints, pPointI)
            {
                dualAverage +=
                    currentPoints[pPoints[pPointI]]
                  + faceCentres[pFaces[pPointI]];
            }
            dualAverage /= 2*nPPoints;

            // Calculate the dual face centre and normal
            vector dualNormal(vector::zero);
            forAll(pPoints, pPointI)
            {
                const label nextPPointI((pPointI + 1) % nPPoints);

                point edgeCentre
                (
                    0.5
                   *(
                        currentPoints[pPoints[pPointI]]
                      + currentPoints[pointI]
                    )
                );
                point nextFaceCentre
                (
                    faceCentres[pFaces[nextPPointI]]
                );
                point nextEdgeCentre
                (
                    0.5
                   *(
                        currentPoints[pPoints[nextPPointI]]
                      + currentPoints[pointI]
                    )
                );

                dualNormal +=
                    (nextFaceCentre - edgeCentre)
                   ^(edgeCentre - dualAverage);
                dualNormal +=
                    (nextEdgeCentre - nextFaceCentre)
                   ^(nextFaceCentre - dualAverage);
            }
            vector dualNormalHat(dualNormal/max(mag(dualNormal), ROOTVSMALL));

            scalar sumA(0);
            vector sumAc(vector::zero);
            forAll(pPoints, pPointI)
            {
                const label nextPPointI((pPointI + 1) % nPPoints);

                point edgeCentre
                (
                    0.5
                   *(
                        currentPoints[pPoints[pPointI]]
                      + currentPoints[pointI]
                    )
                );
                point nextFaceCentre
                (
                    faceCentres[pFaces[nextPPointI]]
                );
                point nextEdgeCentre
                (
                    0.5
                   *(
                        currentPoints[pPoints[nextPPointI]]
                      + currentPoints[pointI]
                    )
                );

                vector c1 = edgeCentre + nextFaceCentre + dualAverage;
                vector c2 = nextFaceCentre + nextEdgeCentre + dualAverage;

                vector n1 =
                    (nextFaceCentre - edgeCentre)
                   ^(edgeCentre - dualAverage);
                vector n2 =
                    (nextEdgeCentre - nextFaceCentre)
                   ^(nextFaceCentre - dualAverage);

                scalar a1 = n1 & dualNormalHat;
                scalar a2 = n2 & dualNormalHat;

                sumA += a1 + a2;
                sumAc += a1*c1 + a2*c2;
            }

            const vector dualCentre(sumAc/max(sumA, ROOTVSMALL)/3);

            // Calculate the transformed point
            transformedPoints[pointI] =
                dualCentre
              + transformationParameter_
               *dualNormal/sqrt(max(mag(dualNormal), ROOTVSMALL));
        }

        // Length scale
        scalar lengthScale(0), transformedLengthScale(0);
        forAll(cEdges, cEdgeI)
        {
            lengthScale +=
                cEdges[cEdgeI].mag(currentPoints);
            transformedLengthScale +=
                cEdges[cEdgeI].mag(transformedPoints);
        }
        lengthScale /= cEdges.size();
        transformedLengthScale /= cEdges.size();

        const scalar lengthScaleRatio =
        (
             (transformedLengthScale > SMALL)
           ? lengthScale/transformedLengthScale
           : scalar(0)
        );

        // Add the displacement to the average
        forAll(cPoints, cPointI)
        {
            const label pointI(cPoints[cPointI]);

            if (counts[pointI] == -1) continue;

            const vector newPoint
            (
                cellCentres[cellI]
              + lengthScaleRatio
               *(
                    transformedPoints[pointI]
                  - cellCentres[cellI]
                )
            );

            ++ counts[pointI];

            pointDisplacement[pointI] += newPoint - oldPoints[pointI];
        }
    }

    // Reset all the boundary faces
    reset(facesToMove, counts, pointDisplacement, false);

    // Calculate the boundary transformations
    forAll(facesToMove, faceToMoveI)
    {
        const label faceI(facesToMove[faceToMoveI]);

        if (!isInternalOrProcessorFace(faceI))
        {
            const labelList& fPoints(mesh().faces()[faceI]);

            // Face normal
            vector faceNormalHat(faceAreas[faceI]);
            faceNormalHat /= max(SMALL, mag(faceNormalHat));

            // Calculate a transformed point for each face point
            forAll(fPoints, fPointI)
            {
                const label pointI(fPoints[fPointI]);

                const label fPointIPrev
                (
                    (fPointI - 1 + fPoints.size()) % fPoints.size()
                );
                const label fPointINext
                (
                    (fPointI + 1) % fPoints.size()
                );

                const vector dualCentre
                (
                    currentPoints[pointI]/2
                  + (
                        currentPoints[fPoints[fPointINext]]
                      + currentPoints[fPoints[fPointIPrev]]
                    )/4
                );
                const vector dualNormal
                (
                    (
                        currentPoints[fPoints[fPointINext]]
                      - currentPoints[fPoints[fPointIPrev]]
                    )/2
                   ^faceNormalHat
                );

                transformedPoints[pointI] =
                    dualCentre
                  + transformationParameter_
                   *dualNormal/sqrt(max(mag(dualNormal), ROOTVSMALL));
            }

            // Length scale
            scalar lengthScale(0), transformedLengthScale(0);
            forAll(fPoints, fPointI)
            {
                const label fPointINext((fPointI + 1)%fPoints.size());

                const label pointI(fPoints[fPointI]);
                const label pointINext(fPoints[fPointINext]);

                lengthScale += mag
                (
                    currentPoints[pointINext] - currentPoints[pointI]
                );
                transformedLengthScale += mag
                (
                    transformedPoints[pointINext] - transformedPoints[pointI]
                );
            }
            lengthScale /= fPoints.size();
            transformedLengthScale /= fPoints.size();

            const scalar lengthScaleRatio
            (
                 (transformedLengthScale > SMALL)
               ? lengthScale/transformedLengthScale
               : scalar(0)
            );

            // Add the displacement to the average
            forAll(fPoints, fPointI)
            {
                const label pointI(fPoints[fPointI]);

                const vector newPoint
                (
                    faceCentres[faceI]
                  + lengthScaleRatio
                   *(
                        transformedPoints[pointI]
                      - faceCentres[faceI]
                    )
                );

                ++ counts[pointI];

                pointDisplacement[pointI] += newPoint - oldPoints[pointI];
            }
        }
    }

    // Average
    average(facesToMove, counts, pointDisplacement);
}


Foam::tmp<Foam::scalarField>
Foam::pointSmoothers::geometricElementTransformPointSmoother::cellQuality
(
    const polyMesh& mesh,
    const pointField& currentPoints
)
{
    const cellPointConnectivity& connectivity =
        MeshObject<polyMesh, MoveableMeshObject, cellPointConnectivity>::New
        (
            mesh
        );

    tmp<scalarField> tfld(tmp<scalarField>::New(mesh.nCells()));
    scalarField& fld = tfld.ref();

    forAll(fld, celli)
    {
        fld[celli] = cellQuality(mesh, currentPoints, connectivity, celli);
    }

    return tfld;
}


// ************************************************************************* //
