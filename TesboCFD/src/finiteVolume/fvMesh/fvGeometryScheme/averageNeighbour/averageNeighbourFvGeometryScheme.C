/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "averageNeighbourFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "cellAspectRatio.H"
#include "syncTools.H"
#include "polyMeshTools.H"
#include "unitConversion.H"
#include "OBJstream.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(averageNeighbourFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        averageNeighbourFvGeometryScheme,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::averageNeighbourFvGeometryScheme::clipFaceTet
(
    const scalar minRatio,
    const vectorField& faceCentres,
    const vectorField& faceNormals,

    vectorField& faceCorrection
) const
{
    // Clip correction vector if any triangle becomes too small. Return number
    // of correction vectors clipped

    const pointField& p = mesh_.points();

    label nClipped = 0;
    for (label facei = 0; facei < mesh_.nFaces(); facei++)
    {
        #ifdef WM_SPDP
        const solveVector fcCorr(faceCorrection[facei]);
        #else
        const vector& fcCorr = faceCorrection[facei];
        #endif
        if (fcCorr != solveVector::zero)
        {
            #ifdef WM_SPDP
            const solveVector fn(faceNormals[facei]);
            const solveVector fc(faceCentres[facei]);
            #else
            const vector& fn = faceNormals[facei];
            const point& fc = faceCentres[facei];
            #endif
            const face& f = mesh_.faces()[facei];

            forAll(f, fp)
            {
                const solveVector thisPt(p[f[fp]]);
                const solveVector nextPt(p[f.fcValue(fp)]);
                const solveVector d(nextPt-thisPt);

                // Calculate triangle area with correction
                const solveVector nCorr(d^(fc+fcCorr - thisPt));

                if ((nCorr & fn) < 0)
                {
                    // Triangle points wrong way
                    faceCorrection[facei] = vector::zero;
                    nClipped++;
                    break;
                }
                else
                {
                    // Calculate triangle area without correction
                    const solveVector n(d^(fc - thisPt));
                    if ((n & fn) < 0)
                    {
                        // Original triangle points the wrong way, new one is ok
                    }
                    else
                    {
                        // Both point correctly. Make sure triangle doesn't get
                        // too small
                        if (mag(nCorr) < minRatio*mag(n))
                        {
                            faceCorrection[facei] = vector::zero;
                            nClipped++;
                            break;
                        }
                    }
                }
            }
        }
    }
    return returnReduce(nClipped, sumOp<label>());
}


void Foam::averageNeighbourFvGeometryScheme::makePyrHeights
(
    const pointField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceNormals,

    scalarField& ownHeight,
    scalarField& neiHeight
) const
{
    ownHeight.setSize(mesh_.nFaces());
    neiHeight.setSize(mesh_.nInternalFaces());

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const solveVector n = faceNormals[facei];
        const solveVector fc = faceCentres[facei];
        const solveVector ownCc = cellCentres[own[facei]];
        const solveVector neiCc = cellCentres[nei[facei]];

        ownHeight[facei] = ((fc-ownCc)&n);
        neiHeight[facei] = ((neiCc-fc)&n);
    }

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const solveVector n = faceNormals[facei];
        const solveVector fc = faceCentres[facei];
        const solveVector ownCc = cellCentres[own[facei]];

        ownHeight[facei] = ((fc-ownCc)&n);
    }
}


Foam::label Foam::averageNeighbourFvGeometryScheme::clipPyramids
(
    const pointField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceNormals,

    const scalarField& minOwnHeight,
    const scalarField& minNeiHeight,

    vectorField& correction
) const
{
    // Clip correction vector if any pyramid becomes too small. Return number of
    // cells clipped

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    label nClipped = 0;
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        #ifdef WM_SPDP
        const solveVector n(faceNormals[facei]);
        const solveVector fc(faceCentres[facei]);
        #else
        const vector& n = faceNormals[facei];
        const point& fc = faceCentres[facei];
        #endif

        const label ownCelli = own[facei];
        if (correction[ownCelli] != vector::zero)
        {
            const solveVector ownCc(cellCentres[ownCelli]+correction[ownCelli]);
            const scalar ownHeight = ((fc-ownCc)&n);
            if (ownHeight < minOwnHeight[facei])
            {
                //Pout<< "    internalface:" << fc
                //    << " own:" << ownCc
                //    << " pyrHeight:" << ownHeight
                //    << " minHeight:" << minOwnHeight[facei]
                //    << endl;
                correction[ownCelli] = vector::zero;
                nClipped++;
            }
        }

        const label neiCelli = nei[facei];
        if (correction[neiCelli] != vector::zero)
        {
            const solveVector neiCc(cellCentres[neiCelli]+correction[neiCelli]);
            const scalar neiHeight = ((neiCc-fc)&n);
            if (neiHeight < minNeiHeight[facei])
            {
                //Pout<< "    internalface:" << fc
                //    << " nei:" << neiCc
                //    << " pyrHeight:" << neiHeight
                //    << " minHeight:" << minNeiHeight[facei]
                //    << endl;
                correction[neiCelli] = vector::zero;
                nClipped++;
            }
        }
    }

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        #ifdef WM_SPDP
        const solveVector n(faceNormals[facei]);
        const solveVector fc(faceCentres[facei]);
        #else
        const vector& n = faceNormals[facei];
        const point& fc = faceCentres[facei];
        #endif

        const label ownCelli = own[facei];
        if (correction[ownCelli] != vector::zero)
        {
            const solveVector ownCc(cellCentres[ownCelli]+correction[ownCelli]);
            const scalar ownHeight = ((fc-ownCc)&n);
            if (ownHeight < minOwnHeight[facei])
            {
                //Pout<< "    boundaryface:" << fc
                //    << " own:" << ownCc
                //    << " pyrHeight:" << ownHeight
                //    << " minHeight:" << minOwnHeight[facei]
                //    << endl;
                correction[ownCelli] = vector::zero;
                nClipped++;
            }
        }
    }
    return returnReduce(nClipped, sumOp<label>());
}


Foam::tmp<Foam::pointField>
Foam::averageNeighbourFvGeometryScheme::averageNeighbourCentres
(
    const pointField& cellCentres,
    const vectorField& faceNormals,
    const scalarField& faceWeights
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    auto tcc = tmp<pointField>::New(mesh_.nCells(), Zero);
    auto& cc = tcc.ref();

    Field<solveScalar> cellWeights(mesh_.nCells(), Zero);

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        #ifdef WM_SPDP
        const solveVector n(faceNormals[facei]);
        #else
        const vector& n = faceNormals[facei];
        #endif
        const point& ownCc = cellCentres[own[facei]];
        const point& neiCc = cellCentres[nei[facei]];

        solveVector d(neiCc-ownCc);

        // 1. Normalise contribution. This increases actual non-ortho
        // since it does not 'see' the tangential offset of neighbours
        //neiCc = ownCc + (d&n)*n;

        // 2. Remove normal contribution, i.e. get tangential vector
        //    (= non-ortho correction vector?)
        d -= (d&n)*n;

        // Apply half to both sides (as a correction)
        // Note: should this be linear weights instead of 0.5?
        const scalar w = 0.5*faceWeights[facei];
        cc[own[facei]] += point(w*d);
        cellWeights[own[facei]] += w;

        cc[nei[facei]] -= point(w*d);
        cellWeights[nei[facei]] += w;
    }


    // Boundary faces. Bypass stored cell centres
    pointField neiCellCentres;
    syncTools::swapBoundaryCellPositions(mesh_, cellCentres, neiCellCentres);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled())
        {
            const labelUList& fc = pp.faceCells();

            forAll(fc, i)
            {
                const label meshFacei = pp.start()+i;
                const label bFacei = meshFacei-mesh_.nInternalFaces();

                #ifdef WM_SPDP
                const solveVector n(faceNormals[meshFacei]);
                #else
                const vector& n = faceNormals[meshFacei];
                #endif

                const point& ownCc = cellCentres[fc[i]];
                const point& neiCc = neiCellCentres[bFacei];

                solveVector d(neiCc-ownCc);

                // 1. Normalise contribution. This increases actual non-ortho
                // since it does not 'see' the tangential offset of neighbours
                //neiCc = ownCc + (d&n)*n;

                // 2. Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                // Apply half to both sides (as a correction)
                const scalar w = 0.5*faceWeights[meshFacei];
                cc[fc[i]] += point(w*d);
                cellWeights[fc[i]] += w;
            }
        }
    }

    // Now cc is still the correction vector. Add to cell original centres.
    forAll(cc, celli)
    {
        if (cellWeights[celli] > VSMALL)
        {
            cc[celli] = cellCentres[celli] + cc[celli]/cellWeights[celli];
        }
        else
        {
            cc[celli] = cellCentres[celli];
        }
    }

    return tcc;
}


Foam::tmp<Foam::pointField>
Foam::averageNeighbourFvGeometryScheme::averageCentres
(
//    const scalar ratio,             // Amount of change in face-triangles area
    const pointField& cellCentres,
    const pointField& faceCentres,
    const vectorField& faceNormals
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    auto tnewFc = tmp<pointField>::New(faceCentres);
    auto& newFc = tnewFc.ref();

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        #ifdef WM_SPDP
        const solveVector n(faceNormals[facei]);
        const solveVector oldFc(faceCentres[facei]);
        #else
        const vector& n = faceNormals[facei];
        const point& oldFc = faceCentres[facei];
        #endif

        const solveVector ownCc(cellCentres[own[facei]]);
        const solveVector neiCc(cellCentres[nei[facei]]);

        solveVector deltaCc(neiCc-ownCc);
        solveVector deltaFc(oldFc-ownCc);

        //solveVector d(neiCc-ownCc);
        //// 1. Normalise contribution. This increases actual non-ortho
        //// since it does not 'see' the tangential offset of neighbours
        ////neiCc = ownCc + s*n;
        //
        //// 2. Remove normal contribution, i.e. get tangential vector
        ////    (= non-ortho correction vector?)
        //d -= s*n;
        //newFc[facei] = faceCentres[facei]+d;

        // Get linear weight (normal distance to face)
        const solveScalar f = (deltaFc&n)/(deltaCc&n);
        const solveVector avgCc((1.0-f)*ownCc + f*neiCc);

        solveVector d(avgCc-oldFc);
        // Remove normal contribution, i.e. get tangential vector
        //    (= non-ortho correction vector?)
        d -= (d&n)*n;

//        // Clip to limit change in
//        d *= ratio;


        newFc[facei] = oldFc + d;
    }


    // Boundary faces. Bypass stored cell centres
    pointField neiCellCentres;
    syncTools::swapBoundaryCellPositions(mesh_, cellCentres, neiCellCentres);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    for (const polyPatch& pp : pbm)
    {
        const labelUList& fc = pp.faceCells();

        if (pp.coupled())
        {
            forAll(fc, i)
            {
                // Same as internal faces
                const label facei = pp.start()+i;
                const label bFacei = facei-mesh_.nInternalFaces();

                #ifdef WM_SPDP
                const solveVector n(faceNormals[facei]);
                const solveVector oldFc(faceCentres[facei]);
                #else
                const vector& n = faceNormals[facei];
                const point& oldFc = faceCentres[facei];
                #endif

                const solveVector ownCc(cellCentres[fc[i]]);
                const solveVector neiCc(neiCellCentres[bFacei]);

                solveVector deltaCc(neiCc-ownCc);
                solveVector deltaFc(oldFc-ownCc);

                // Get linear weight (normal distance to face)
                const solveScalar f = (deltaFc&n)/(deltaCc&n);
                const solveVector avgCc((1.0-f)*ownCc + f*neiCc);

                solveVector d(avgCc-oldFc);
                // Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                newFc[facei] = oldFc + d;
            }
        }
        else
        {
            // Zero-grad?
            forAll(fc, i)
            {
                const label facei = pp.start()+i;

                #ifdef WM_SPDP
                const solveVector n(faceNormals[facei]);
                const solveVector oldFc(faceCentres[facei]);
                #else
                const vector& n = faceNormals[facei];
                const point& oldFc = faceCentres[facei];
                #endif

                const solveVector ownCc(cellCentres[fc[i]]);

                solveVector d(ownCc-oldFc);
                // Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                newFc[facei] = oldFc+d;
            }
        }
    }

    return tnewFc;
}


void Foam::averageNeighbourFvGeometryScheme::makeNonOrthoWeights
(
    const pointField& cellCentres,
    const vectorField& faceNormals,

    scalarField& cosAngles,
    scalarField& faceWeights
) const
{
    cosAngles = clamp
    (
        polyMeshTools::faceOrthogonality(mesh_, faceNormals, cellCentres),
        zero_one{}
    );

    // Make weight: 0 for ortho faces, 1 for 90degrees non-ortho
    //const scalarField faceWeights(scalar(1)-cosAngles);
    faceWeights.setSize(cosAngles.size());
    {
        const scalar minCos = Foam::cos(degToRad(80));
        const scalar maxCos = Foam::cos(degToRad(10));

        forAll(cosAngles, facei)
        {
            const scalar cosAngle = cosAngles[facei];
            if (cosAngle < minCos)
            {
                faceWeights[facei] = 1.0;
            }
            else if (cosAngle > maxCos)
            {
                faceWeights[facei] = 0.0;
            }
            else
            {
                faceWeights[facei] =
                    1.0-(cosAngle-minCos)/(maxCos-minCos);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::averageNeighbourFvGeometryScheme::averageNeighbourFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    highAspectRatioFvGeometryScheme(mesh, dict),
    nIters_
    (
        dict.getCheckOrDefault<label>
        (
            "nIters",
            1,
            [](label val) { return val >= 0; }
        )
    ),
    relax_
    (
        dict.getCheck<scalar>
        (
            "relax",
            [](scalar val) { return val > 0 && val <= 1; }
        )
    ),
    minRatio_
    (
        dict.getCheckOrDefault<scalar>
        (
            "minRatio",
            0.5,
            [](scalar val) { return val >= 0 && val <= 1; }
        )
    )
{
    if (debug)
    {
        Pout<< "averageNeighbourFvGeometryScheme :"
            << " nIters:" << nIters_
            << " relax:" << relax_
            << " minRatio:" << minRatio_ << endl;
    }

    // Force local calculation
    movePoints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::averageNeighbourFvGeometryScheme::movePoints()
{
    if (debug)
    {
        Pout<< "averageNeighbourFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    //if
    //(
    //   !mesh_.hasCellCentres()
    //&& !mesh_.hasFaceCentres()
    //&& !mesh_.hasCellVolumes()
    //&& !mesh_.hasFaceAreas()
    //)
    {
        highAspectRatioFvGeometryScheme::movePoints();

        // Note: at this point the highAspectRatioFvGeometryScheme constructor
        //       will have already reset the primitive geometry!

        vectorField faceAreas(mesh_.faceAreas());
        const scalarField magFaceAreas(mag(faceAreas));
        const vectorField faceNormals(faceAreas/magFaceAreas);


        // Calculate aspectratio weights
        // - 0 if aratio < minAspect_
        // - 1 if aratio >= maxAspect_
        scalarField cellWeight, faceWeight;
        calcAspectRatioWeights(cellWeight, faceWeight);

        // Relaxation
        cellWeight *= relax_;
        //faceWeight *= relax_;

        // Calculate current pyramid heights
        scalarField minOwnHeight;
        scalarField minNeiHeight;
        makePyrHeights
        (
            mesh_.cellCentres(),
            mesh_.faceCentres(),
            faceNormals,

            minOwnHeight,
            minNeiHeight
        );

        // How much is the cell centre to vary inside the cell.
        minOwnHeight *= minRatio_;
        minNeiHeight *= minRatio_;



        autoPtr<OBJstream> osPtr;
        autoPtr<surfaceWriter> writerPtr;
        if (debug)
        {
            osPtr.reset
            (
                new OBJstream
                (
                    mesh_.time().timePath()
                  / "cellCentre_correction.obj"
                )
            );
            Pout<< "averageNeighbourFvGeometryScheme::movePoints() :"
                << " writing cell centre path to " << osPtr().name() << endl;


            // Write current non-ortho
            fileName outputDir
            (
                mesh_.time().globalPath()
              / functionObject::outputPrefix
              / mesh_.pointsInstance()
            );
            outputDir.clean();  // Remove unneeded ".."
            writerPtr = surfaceWriter::New
            (
                "ensight" //"vtk"
                // options
            );

            // Use outputDir/TIME/surface-name
            writerPtr->useTimeDir(true);

            writerPtr->beginTime(mesh_.time());

            writerPtr->open
            (
                mesh_.points(),
                mesh_.faces(),
                (outputDir / "cosAngle"),
                true  // merge parallel bits
            );

            writerPtr->endTime();
        }


        // Current cellCentres. These get adjusted to lower the
        // non-orthogonality
        pointField cellCentres(mesh_.cellCentres());

        // Modify cell centres to be more in-line with the face normals
        for (label iter = 0; iter < nIters_; iter++)
        {
            // Get neighbour average (weighted by face area). This gives
            // preference to the dominant faces. However if the non-ortho
            // is not caused by the dominant faces this moves to the wrong
            // direction.
            //tmp<pointField> tcc
            //(
            //    averageNeighbourCentres
            //    (
            //        cellCentres,
            //        faceNormals,
            //        magFaceAreas
            //    )
            //);

            // Get neighbour average weighted by non-ortho. Question: how to
            // weight boundary faces?
            tmp<pointField> tcc;
            {
                scalarField cosAngles;
                scalarField faceWeights;
                makeNonOrthoWeights
                (
                    cellCentres,
                    faceNormals,

                    cosAngles,
                    faceWeights
                );

                if (writerPtr)
                {
                    writerPtr->beginTime(instant(scalar(iter)));
                    writerPtr->write("cosAngles", cosAngles);
                    writerPtr->endTime();
                }

                if (debug)
                {
                    forAll(cosAngles, facei)
                    {
                        if (cosAngles[facei] < Foam::cos(degToRad(85.0)))
                        {
                            Pout<< "    face:" << facei
                                << " at:" << mesh_.faceCentres()[facei]
                                << " cosAngle:" << cosAngles[facei]
                                << " faceWeight:" << faceWeights[facei]
                                << endl;
                        }
                    }
                }

                tcc = averageNeighbourCentres
                (
                    cellCentres,
                    faceNormals,
                    faceWeights
                );
            }


            // Calculate correction for cell centres. Leave low-aspect
            // ratio cells unaffected (i.e. correction = 0)
            vectorField correction(cellWeight*(tcc-cellCentres));

            // Clip correction vector if pyramid becomes too small
            const label nClipped = clipPyramids
            (
                cellCentres,
                mesh_.faceCentres(),
                faceNormals,

                minOwnHeight, // minimum owner pyramid height. Usually fraction
                minNeiHeight, // of starting mesh

                correction
            );

            cellCentres += correction;

            if (debug)
            {
                forAll(cellCentres, celli)
                {
                    const point& cc = cellCentres[celli];
                    osPtr().writeLine(cc-correction[celli], cc);
                }

                const scalarField magCorrection(mag(correction));
                const scalarField faceOrthogonality
                (
                    min
                    (
                        scalar(1),
                        polyMeshTools::faceOrthogonality
                        (
                            mesh_,
                            faceAreas,
                            cellCentres
                        )
                    )
                );
                const scalarField nonOrthoAngle
                (
                    radToDeg
                    (
                        Foam::acos(faceOrthogonality)
                    )
                );
                Pout<< "    iter:" << iter
                    << " nClipped:" << nClipped
                    << " average displacement:" << gAverage(magCorrection)
                    << " non-ortho angle : average:" << gAverage(nonOrthoAngle)
                    << " max:" << gMax(nonOrthoAngle) << endl;
            }
        }

        tmp<pointField> tfc
        (
            averageCentres
            (
                cellCentres,
                mesh_.faceCentres(),
                faceNormals
            )
        );
        vectorField faceCorrection(faceWeight*(tfc-mesh_.faceCentres()));
        // Clip correction vector to not have min triangle shrink
        // by more than minRatio
        clipFaceTet
        (
            minRatio_,
            mesh_.faceCentres(),
            faceNormals,
            faceCorrection
        );
        vectorField faceCentres(mesh_.faceCentres() + faceCorrection);

        if (debug)
        {
            Pout<< "averageNeighbourFvGeometryScheme::movePoints() :"
                << " averageNeighbour weight"
                << " max:" << gMax(cellWeight) << " min:" << gMin(cellWeight)
                << " average:" << gAverage(cellWeight) << endl;

            // Dump lines from old to new location
            const fileName tp(mesh_.time().timePath());
            mkDir(tp);
            OBJstream str(tp/"averageNeighbourCellCentres.obj");
            Pout<< "Writing lines from old to new cell centre to " << str.name()
                << endl;
            forAll(mesh_.cellCentres(), celli)
            {
                const point& oldCc = mesh_.cellCentres()[celli];
                const point& newCc = cellCentres[celli];
                str.writeLine(oldCc, newCc);
            }
        }
        if (debug)
        {
            // Dump lines from old to new location
            const fileName tp(mesh_.time().timePath());
            OBJstream str(tp/"averageFaceCentres.obj");
            Pout<< "Writing lines from old to new face centre to " << str.name()
                << endl;
            forAll(mesh_.faceCentres(), facei)
            {
                const point& oldFc = mesh_.faceCentres()[facei];
                const point& newFc = faceCentres[facei];
                str.writeLine(oldFc, newFc);
            }
        }

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


// ************************************************************************* //
