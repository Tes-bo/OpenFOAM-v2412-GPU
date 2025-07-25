/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class Type>
void Foam::fieldSmoother::minSmoothField
(
    const label nIter,
    const bitSet& isPatchMasterPoint,
    const bitSet& isPatchMasterEdge,
    const indirectPrimitivePatch& adaptPatch,
    const scalarField& fieldMinMag,
    Field<Type>& field
) const
{
    const edgeList& edges = adaptPatch.edges();
    const labelList& meshPoints = adaptPatch.meshPoints();

    scalarField edgeWeights(edges.size());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh_,
        isPatchMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< typeName << " : Smoothing field ..." << endl;

    for (label iter = 0; iter < nIter; iter++)
    {
        Field<Type> average(adaptPatch.nPoints());
        meshRefinement::weightedSum
        (
            mesh_,
            isPatchMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            field,
            average
        );
        average *= invSumWeight;

        // Transfer to field
        forAll(field, pointI)
        {
            //full smoothing neighbours + point value
            average[pointI] = 0.5*(field[pointI]+average[pointI]);

            // perform monotonic smoothing
            if
            (
                mag(average[pointI]) < mag(field[pointI])
             && mag(average[pointI]) >= mag(fieldMinMag[pointI])
            )
            {
                field[pointI] = average[pointI];
            }
        }

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isPatchMasterPoint,
                mag(field-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}


// ************************************************************************* //
