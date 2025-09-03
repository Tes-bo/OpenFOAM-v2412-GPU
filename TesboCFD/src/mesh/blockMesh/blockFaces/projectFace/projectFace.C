/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "projectFace.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "blockDescriptor.H"
#include "OBJstream.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{
    defineTypeNameAndDebug(projectFace, 0);
    addToRunTimeSelectionTable(blockFace, projectFace, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::searchableSurface& Foam::blockFaces::projectFace::lookupSurface
(
    const searchableSurfaces& geometry,
    Istream& is
) const
{
    const word name(is);

    for (const searchableSurface& geom : geometry)
    {
        if (geom.name() == name)
        {
            return geom;
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot find surface " << name << " in geometry"
        << exit(FatalIOError);

    return geometry[0];
}


Foam::label Foam::blockFaces::projectFace::index
(
    const labelPair& n,
    const labelPair& coord
)
{
    return coord.first() + coord.second()*n.first();
}


void Foam::blockFaces::projectFace::calcLambdas
(
    const labelPair& n,
    const pointField& points,
    scalarField& lambdaI,
    scalarField& lambdaJ
) const
{
    lambdaI.setSize(points.size());
    lambdaI = 0.0;
    lambdaJ.setSize(points.size());
    lambdaJ = 0.0;

    for (label i = 1; i < n.first(); i++)
    {
        for (label j = 1; j < n.second(); j++)
        {
            label ij = index(n, labelPair(i, j));
            label iMin1j = index(n, labelPair(i-1, j));
            lambdaI[ij] = lambdaI[iMin1j] + mag(points[ij]-points[iMin1j]);

            label ijMin1 = index(n, labelPair(i, j-1));
            lambdaJ[ij] = lambdaJ[ijMin1] + mag(points[ij]-points[ijMin1]);
        }
    }

    for (label i = 1; i < n.first(); i++)
    {
        label ijLast = index(n, labelPair(i, n.second()-1));
        for (label j = 1; j < n.second(); j++)
        {
            label ij = index(n, labelPair(i, j));
            lambdaJ[ij] /= lambdaJ[ijLast];
        }
    }
    for (label j = 1; j < n.second(); j++)
    {
        label iLastj = index(n, labelPair(n.first()-1, j));
        for (label i = 1; i < n.first(); i++)
        {
            label ij = index(n, labelPair(i, j));
            lambdaI[ij] /= lambdaI[iLastj];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFaces::projectFace::projectFace
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
:
    blockFace(dict, index, is),
    surface_(lookupSurface(geometry, is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFaces::projectFace::project
(
    const blockDescriptor& desc,
    const label blockFacei,
    pointField& points
) const
{
    // For debugging to tag the output
    static label fIter = 0;

    autoPtr<OBJstream> debugStr;
    if (debug)
    {
        debugStr.reset
        (
            new OBJstream("projectFace_" + Foam::name(fIter++) + ".obj")
        );
        Info<< "Face:" << blockFacei << " on block:" << desc.blockShape()
            << " with verts:" << desc.vertices()
            << " writing lines from start points"
            << " to projected points to " << debugStr().name() << endl;
    }


    // Find out range of vertices in face
    labelPair n(-1, -1);
    switch (blockFacei)
    {
        case 0:
        case 1:
        {
            n.first() = desc.density().y() + 1;
            n.second() = desc.density().z() + 1;
        }
        break;

        case 2:
        case 3:
        {
            n.first() = desc.density().x() + 1;
            n.second() = desc.density().z() + 1;
        }
        break;

        case 4:
        case 5:
        {
            n.first() = desc.density().x() + 1;
            n.second() = desc.density().y() + 1;
        }
        break;
    }


    // Calculate initial normalised edge lengths (= u,v coordinates)
    scalarField lambdaI(points.size(), Zero);
    scalarField lambdaJ(points.size(), Zero);
    calcLambdas(n, points, lambdaI, lambdaJ);


    // Upper limit for number of iterations
    const label maxIter = 10;
    // Residual tolerance
    const scalar relTol = 0.1;

    scalar initialResidual = 0.0;
    scalar iResidual = 0.0;
    scalar jResidual = 0.0;

    for (label iter = 0; iter < maxIter;  iter++)
    {
        // Do projection
        {
            List<pointIndexHit> hits;
            scalarField nearestDistSqr
            (
                points.size(),
                magSqr(points[0] - points[points.size()-1])
            );
            surface_.findNearest(points, nearestDistSqr, hits);

            forAll(hits, i)
            {
                if (hits[i].hit())
                {
                    if (debugStr)
                    {
                        debugStr().writeLine(points[i], hits[i].point());
                    }
                    points[i] = hits[i].point();
                }
            }
        }

        if (debug)
        {
            Pout<< "Iter:" << iter << " initialResidual:" << initialResidual
                << " iResidual+jResidual:" << iResidual+jResidual << endl;
        }


        if
        (
            iter > 0
         && (
                initialResidual < ROOTVSMALL
             || ((iResidual+jResidual)/initialResidual < relTol)
            )
        )
        {
            break;
        }


        // Predict along i
        vectorField residual(points.size(), Zero);

        // Work arrays for interpolation
        labelList indices;
        scalarField weights;
        for (label j = 1; j < n.second()-1; j++)
        {
            // Calculate actual lamdba along constant j line
            scalarField projLambdas(n.first());
            projLambdas[0] = 0.0;
            for (label i = 1; i < n.first(); i++)
            {
                label ij = index(n, labelPair(i, j));
                label iMin1j = index(n, labelPair(i-1, j));
                projLambdas[i] =
                    projLambdas[i-1]
                   +mag(points[ij]-points[iMin1j]);
            }
            projLambdas /= projLambdas.last();

            linearInterpolationWeights interpolator(projLambdas);

            for (label i = 1; i < n.first()-1; i++)
            {
                label ij = index(n, labelPair(i, j));

                interpolator.valueWeights(lambdaI[ij], indices, weights);

                point predicted(Zero);
                forAll(indices, indexi)
                {
                    label ptIndex = index(n, labelPair(indices[indexi], j));
                    predicted += weights[indexi]*points[ptIndex];
                }
                residual[ij] = predicted-points[ij];
            }
        }

        if (debugStr)
        {
            forAll(points, i)
            {
                const point predicted(points[i] + residual[i]);
                debugStr().writeLine(points[i], predicted);
            }
        }

        iResidual = sum(mag(residual));

        // Update points before doing j. Note: is this needed? Complicates
        // residual checking.
        points += residual;


        // Predict along j
        residual = vector::zero;
        for (label i = 1; i < n.first()-1; i++)
        {
            // Calculate actual lamdba along constant i line
            scalarField projLambdas(n.second());
            projLambdas[0] = 0.0;
            for (label j = 1; j < n.second(); j++)
            {
                label ij = index(n, labelPair(i, j));
                label ijMin1 = index(n, labelPair(i, j-1));
                projLambdas[j] =
                    projLambdas[j-1]
                   +mag(points[ij]-points[ijMin1]);
            }

            projLambdas /= projLambdas.last();

            linearInterpolationWeights interpolator(projLambdas);

            for (label j = 1; j < n.second()-1; j++)
            {
                label ij = index(n, labelPair(i, j));

                interpolator.valueWeights(lambdaJ[ij], indices, weights);

                point predicted(Zero);
                forAll(indices, indexi)
                {
                    label ptIndex = index(n, labelPair(i, indices[indexi]));
                    predicted += weights[indexi]*points[ptIndex];
                }
                residual[ij] = predicted-points[ij];
            }
        }

        if (debugStr)
        {
            forAll(points, i)
            {
                const point predicted(points[i] + residual[i]);
                debugStr().writeLine(points[i], predicted);
            }
        }

        jResidual = sum(mag(residual));

        if (iter == 0)
        {
            initialResidual = iResidual + jResidual;
        }

        points += residual;
    }
}


// ************************************************************************* //
