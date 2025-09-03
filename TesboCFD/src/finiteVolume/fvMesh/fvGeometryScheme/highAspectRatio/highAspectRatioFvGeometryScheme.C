/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "highAspectRatioFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "triangle.H"
#include "syncTools.H"
#include "cellAspectRatio.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(highAspectRatioFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        highAspectRatioFvGeometryScheme,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//void Foam::highAspectRatioFvGeometryScheme::cellClosedness
//(
//    const vectorField& areas,
//    const scalarField& vols,
//    const tensorField& cellCoords,
//
//    scalarField& aratio
//) const
//{
//    // From primitiveMeshTools::cellClosedness:
//    //  calculate aspect ratio in given direction
//    const labelList& own = mesh_.faceOwner();
//    const labelList& nei = mesh_.faceNeighbour();
//
//    // Loop through cell faces and sum up the face area vectors for each cell.
//    // This should be zero in all vector components
//
//    vectorField sumMagClosed(mesh_.nCells(), Zero);
//
//    forAll(own, facei)
//    {
//        // Add to owner
//        vector& v = sumMagClosed[own[facei]];
//        v += cmptMag(cellCoords[own[facei]] & areas[facei]);
//    }
//
//    forAll(nei, facei)
//    {
//        // Subtract from neighbour
//        vector& v = sumMagClosed[nei[facei]];
//        v += cmptMag(cellCoords[nei[facei]] & areas[facei]);
//    }
//
//    // Check the sums
//    aratio.setSize(mesh_.nCells());
//
//    forAll(sumMagClosed, celli)
//    {
//        // Calculate the aspect ration as the maximum of Cartesian component
//        // aspect ratio to the total area hydraulic area aspect ratio
//        scalar minCmpt = VGREAT;
//        scalar maxCmpt = -VGREAT;
//        for (direction dir = 0; dir < vector::nComponents; dir++)
//        {
//            minCmpt = min(minCmpt, sumMagClosed[celli][dir]);
//            maxCmpt = max(maxCmpt, sumMagClosed[celli][dir]);
//        }
//
//        scalar aspectRatio = maxCmpt/(minCmpt + ROOTVSMALL);
//        const scalar v = max(ROOTVSMALL, vols[celli]);
//
//        aspectRatio = max
//        (
//            aspectRatio,
//            1.0/6.0*cmptSum(sumMagClosed[celli])/Foam::pow(v, 2.0/3.0)
//        );
//
//        aratio[celli] = aspectRatio;
//    }
//}
//
//
//void Foam::highAspectRatioFvGeometryScheme::cellDirections
//(
//    tensorField& T,
//    vectorField& lambda
//) const
//{
//    // Calculate principal directions in increasing order
//
//    T.setSize(mesh_.nCells());
//    lambda.setSize(mesh_.nCells());
//
//    forAll(T, celli)
//    {
//        tensor J = Zero;
//        {
//            const List<tetIndices> cellTets
//            (
//                polyMeshTetDecomposition::cellTetIndices
//                (
//                    mesh_,
//                    celli
//                )
//            );
//            triFaceList faces(cellTets.size());
//            forAll(cellTets, cTI)
//            {
//                faces[cTI] = cellTets[cTI].faceTriIs(mesh_);
//            }
//
//            scalar m = 0.0;
//            vector cM = Zero;
//            J = Zero;
//            momentOfInertia::massPropertiesShell
//            (
//                mesh_.points(),
//                faces,
//                1.0,
//                m,
//                cM,
//                J
//            );
//        }
//
//        lambda[celli] = Foam::eigenValues(J);
//        T[celli] = Foam::eigenVectors(J, lambda[celli]);
//    }
//}

void Foam::highAspectRatioFvGeometryScheme::calcAspectRatioWeights
(
    scalarField& cellWeight,
    scalarField& faceWeight
) const
{
    //scalarField aratio;
    //{
    //    tensorField principalDirections;
    //    vectorField lambdas;
    //    cellDirections(principalDirections, lambdas);
    //
    //    cellClosedness
    //    (
    //        mesh_.faceAreas(),
    //        mesh_.cellVolumes(),
    //        principalDirections,
    //        aratio
    //    );
    //}
    const cellAspectRatio aratio(mesh_);

    // Weighting for correction
    // - 0 if aratio < minAspect_
    // - 1 if aratio >= maxAspect_

    scalar delta(maxAspect_-minAspect_);
    if (delta < ROOTVSMALL)
    {
        delta = SMALL;
    }

    cellWeight = clamp((aratio-minAspect_)/delta, zero_one{});

    faceWeight.setSize(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const label own = mesh_.faceOwner()[facei];
        const label nei = mesh_.faceNeighbour()[facei];
        faceWeight[facei] = max(cellWeight[own], cellWeight[nei]);
    }
    scalarField nbrCellWeight;
    syncTools::swapBoundaryCellList
    (
        mesh_,
        cellWeight,
        nbrCellWeight
    );
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        const label own = mesh_.faceOwner()[facei];
        const label bFacei = facei-mesh_.nInternalFaces();
        faceWeight[facei] = max(cellWeight[own], nbrCellWeight[bFacei]);
    }
}


void Foam::highAspectRatioFvGeometryScheme::makeAverageCentres
(
    const polyMesh& mesh,
    const pointField& p,
    const pointField& faceAreas,
    const scalarField& magFaceAreas,
    pointField& faceCentres,
    pointField& cellCentres
)
{
    if (debug)
    {
        Pout<< "highAspectRatioFvGeometryScheme::makeAverageCentres() : "
            << "calculating weighted average face/cell centre" << endl;
    }

    const faceList& fs = mesh.faces();

    // Start off from primitiveMesh faceCentres (preserved for triangles)
    faceCentres.setSize(mesh.nFaces());

    forAll(fs, facei)
    {
        const labelList& f = fs[facei];
        const label nPoints = f.size();

        if (nPoints == 3)
        {
            faceCentres[facei] = triPointRef::centre(p[f[0]], p[f[1]], p[f[2]]);
        }
        else
        {
            solveScalar sumA = 0.0;
            solveVector sumAc = Zero;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const label nextPi(pi == nPoints-1 ? 0 : pi+1);
                const solveVector nextPoint(p[f[nextPi]]);
                const solveVector thisPoint(p[f[pi]]);

                const solveVector eMid = 0.5*(thisPoint+nextPoint);
                const solveScalar a = mag(nextPoint-thisPoint);

                sumAc += a*eMid;
                sumA += a;
            }
            // This is to deal with zero-area faces. Mark very small faces
            // to be detected in e.g. processorPolyPatch.
            if (sumA >= ROOTVSMALL)
            {
                faceCentres[facei] = sumAc/sumA;
            }
            else
            {
                // Unweighted average of points
                sumAc = Zero;
                for (label pi = 0; pi < nPoints; pi++)
                {
                    sumAc += static_cast<solveVector>(p[f[pi]]);
                }
                faceCentres[facei] = sumAc/nPoints;
            }
        }
    }


    cellCentres.setSize(mesh.nCells());
    cellCentres = Zero;
    {
        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();

        Field<solveScalar> cellWeights(mesh.nCells(), Zero);
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            const solveScalar magfA(magFaceAreas[facei]);
            const vector weightedFc(magfA*faceCentres[facei]);

            // Accumulate area-weighted face-centre
            cellCentres[own[facei]] += weightedFc;
            cellCentres[nei[facei]] += weightedFc;

            // Accumulate weights
            cellWeights[own[facei]] += magfA;
            cellWeights[nei[facei]] += magfA;
        }

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        for (const polyPatch& pp : pbm)
        {
            if (!isA<emptyPolyPatch>(pp) && !isA<wedgePolyPatch>(pp))
            {
                for
                (
                    label facei = pp.start();
                    facei < pp.start()+pp.size();
                    facei++
                )
                {
                    const solveScalar magfA(magFaceAreas[facei]);
                    const vector weightedFc(magfA*faceCentres[facei]);

                    cellCentres[own[facei]] += weightedFc;
                    cellWeights[own[facei]] += magfA;
                }
            }
        }

        forAll(cellCentres, celli)
        {
            if (mag(cellWeights[celli]) > VSMALL)
            {
                cellCentres[celli] /= cellWeights[celli];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highAspectRatioFvGeometryScheme::highAspectRatioFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicFvGeometryScheme(mesh, dict),
    minAspect_(dict.get<scalar>("minAspect")),
    maxAspect_(dict.get<scalar>("maxAspect"))
{
    if (maxAspect_ < minAspect_)
    {
        FatalIOErrorInFunction(dict)
            << "minAspect " << minAspect_
            << " has to be less than maxAspect " << maxAspect_
            << exit(FatalIOError);
    }
    if (minAspect_ < 0 || maxAspect_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Illegal aspect ratio : minAspect:" << minAspect_
            << " maxAspect:" << maxAspect_
            << exit(FatalIOError);
    }

    // Force local calculation
    movePoints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::highAspectRatioFvGeometryScheme::movePoints()
{
    //basicFvGeometryScheme::movePoints();
    fvGeometryScheme::movePoints();

    if (debug)
    {
        Pout<< "highAspectRatioFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    if
    (
       !mesh_.hasCellCentres()
    && !mesh_.hasFaceCentres()
    && !mesh_.hasCellVolumes()
    && !mesh_.hasFaceAreas()
    )
    {
        // Use lower level to calculate the geometry
        const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();

        pointField avgFaceCentres;
        pointField avgCellCentres;
        makeAverageCentres
        (
            mesh_,
            mesh_.points(),
            mesh_.faceAreas(),
            mag(mesh_.faceAreas()),
            avgFaceCentres,
            avgCellCentres
        );


        // Calculate aspectratio weights
        // - 0 if aratio < minAspect_
        // - 1 if aratio >= maxAspect_
        scalarField cellWeight, faceWeight;
        calcAspectRatioWeights(cellWeight, faceWeight);


        // Weight with average ones
        vectorField faceCentres
        (
            (1.0-faceWeight)*mesh_.faceCentres()
          + faceWeight*avgFaceCentres
        );
        vectorField cellCentres
        (
            (1.0-cellWeight)*mesh_.cellCentres()
          + cellWeight*avgCellCentres
        );


        if (debug)
        {
            Pout<< "highAspectRatioFvGeometryScheme::movePoints() :"
                << " highAspectRatio weight"
                << " max:" << gMax(cellWeight) << " min:" << gMin(cellWeight)
                << " average:" << gAverage(cellWeight) << endl;
        }

        vectorField faceAreas(mesh_.faceAreas());
        scalarField cellVolumes(mesh_.cellVolumes());

        // Store on primitiveMesh
        //const_cast<fvMesh&>(mesh_).clearGeom();
        const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
        (
            std::move(faceCentres),
            std::move(faceAreas),
            std::move(cellCentres),
            std::move(cellVolumes)
        );
    }
}


bool Foam::highAspectRatioFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    // Basic
    basicFvGeometryScheme::updateGeom
    (
        points,
        oldPoints,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes
    );

    // Average opposite faces
    pointField avgFaceCentres;
    pointField avgCellCentres;
    makeAverageCentres
    (
        mesh_,
        points,
        faceAreas,
        mag(faceAreas),
        avgFaceCentres,
        avgCellCentres
    );

    faceCentres = std::move(avgFaceCentres);
    cellCentres = std::move(avgCellCentres);

    // Assume something has changed.
    return true;
}


// ************************************************************************* //
