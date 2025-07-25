/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "polyMeshGeometry.H"
#include "polyMeshTetDecomposition.H"
#include "pyramid.H"
#include "tetrahedron.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshGeometry, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMeshGeometry::updateFaceCentresAndAreas
(
    const pointField& p,
    const labelList& changedFaces
)
{
    const faceList& fcs = mesh_.faces();

    for (const label facei : changedFaces)
    {
        const face& f = fcs[facei];
        const label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            faceCentres_[facei] =
                triPointRef::centre(p[f[0]], p[f[1]], p[f[2]]);
            faceAreas_[facei] =
                triPointRef::areaNormal(p[f[0]], p[f[1]], p[f[2]]);
        }
        else
        {
            solveVector sumN = Zero;
            solveScalar sumA = Zero;
            solveVector sumAc = Zero;

            solveVector fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; ++pi)
            {
                fCentre += solveVector(p[f[pi]]);
            }
            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; ++pi)
            {
                const solveVector thisPoint(p[f.thisLabel(pi)]);
                const solveVector nextPoint(p[f.nextLabel(pi)]);

                solveVector c = thisPoint + nextPoint + fCentre;
                solveVector n = (nextPoint - thisPoint)^(fCentre - thisPoint);
                solveScalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            if (sumA < ROOTVSMALL)
            {
                faceCentres_[facei] = fCentre;
                faceAreas_[facei] = Zero;
            }
            else
            {
                faceCentres_[facei] = (1.0/3.0)*sumAc/sumA;
                faceAreas_[facei] = 0.5*sumN;
            }
        }
    }
}


void Foam::polyMeshGeometry::updateCellCentresAndVols
(
    const labelList& changedCells,
    const labelList& changedFaces
)
{
    const labelList& own = mesh().faceOwner();
    const cellList& cells = mesh().cells();

    // Clear the fields for accumulation
    UIndirectList<vector>(cellCentres_, changedCells) = Zero;
    UIndirectList<scalar>(cellVolumes_, changedCells) = Zero;


    // Re-calculate the changed cell centres and volumes
    for (const label celli : changedCells)
    {
        const labelList& cFaces = cells[celli];

        // Estimate the cell centre and bounding box using the face centres
        vector cEst(Zero);
        boundBox bb;

        for (const label facei : cFaces)
        {
            const point& fc = faceCentres_[facei];
            cEst += fc;
            bb.add(fc);
        }
        cEst /= cFaces.size();


        // Sum up the face-pyramid contributions
        for (const label facei : cFaces)
        {
            // Calculate 3* the face-pyramid volume
            scalar pyr3Vol = faceAreas_[facei] & (faceCentres_[facei] - cEst);

            if (own[facei] != celli)
            {
                pyr3Vol = -pyr3Vol;
            }

            // Accumulate face-pyramid volume
            cellVolumes_[celli] += pyr3Vol;

            // Calculate the face-pyramid centre
            const vector pCtr = (3.0/4.0)*faceCentres_[facei] + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            cellCentres_[celli] += pyr3Vol*pCtr;
        }

        // Average the accumulated quantities

        if (mag(cellVolumes_[celli]) > VSMALL)
        {
            point cc = cellCentres_[celli] / cellVolumes_[celli];

            // Do additional check for collapsed cells since some volumes
            // (e.g. 1e-33) do not trigger above but do return completely
            // wrong cell centre
            if (bb.contains(cc))
            {
                cellCentres_[celli] = cc;
            }
            else
            {
                cellCentres_[celli] = cEst;
            }
        }
        else
        {
            cellCentres_[celli] = cEst;
        }

        cellVolumes_[celli] *= (1.0/3.0);
    }
}


Foam::labelList Foam::polyMeshGeometry::affectedCells
(
    const polyMesh& mesh,
    const labelList& changedFaces
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    labelHashSet affectedCells(2*changedFaces.size());

    for (const label facei : changedFaces)
    {
        affectedCells.insert(own[facei]);

        if (mesh.isInternalFace(facei))
        {
            affectedCells.insert(nei[facei]);
        }
    }
    return affectedCells.toc();
}


Foam::scalar Foam::polyMeshGeometry::checkNonOrtho
(
    const polyMesh& mesh,
    const bool report,
    const scalar severeNonorthogonalityThreshold,
    const label facei,
    const vector& s,    // face area vector
    const vector& d,    // cc-cc vector

    label& severeNonOrth,
    label& errorNonOrth,
    labelHashSet* setPtr
)
{
    scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

    if (dDotS < severeNonorthogonalityThreshold)
    {
        label nei = -1;

        if (mesh.isInternalFace(facei))
        {
            nei = mesh.faceNeighbour()[facei];
        }

        if (dDotS > SMALL)
        {
            if (report)
            {
                // Severe non-orthogonality but mesh still OK
                Pout<< "Severe non-orthogonality for face " << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            ++severeNonOrth;
        }
        else
        {
            // Non-orthogonality greater than 90 deg
            if (report)
            {
                WarningInFunction
                    << "Severe non-orthogonality detected for face "
                    << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            ++errorNonOrth;
        }

        if (setPtr)
        {
            setPtr->insert(facei);
        }
    }
    return dDotS;
}


// Create the neighbour pyramid - it will have positive volume
bool Foam::polyMeshGeometry::checkFaceTet
(
    const polyMesh& mesh,
    const bool report,
    const scalar minTetQuality,
    const pointField& p,
    const label facei,
    const point& fc,    // face centre
    const point& cc,    // cell centre

    labelHashSet* setPtr
)
{
    const face& f = mesh.faces()[facei];

    forAll(f, fp)
    {
        scalar tetQual = tetPointRef
        (
            p[f[fp]],
            p[f.nextLabel(fp)],
            fc,
            cc
        ).quality();

        if (tetQual < minTetQuality)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFaceTets("
                    << "const bool, const scalar, const pointField&"
                    << ", const pointField&"
                    << ", const labelList&, labelHashSet*) : "
                    << "face " << facei
                    << " has a triangle that points the wrong way." << nl
                    << "Tet quality: " << tetQual
                    << " Face " << facei
                    << endl;
            }
            if (setPtr)
            {
                setPtr->insert(facei);
            }
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMeshGeometry::polyMeshGeometry(const polyMesh& mesh)
:
    mesh_(mesh)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMeshGeometry::correct()
{
    faceAreas_ = mesh_.faceAreas();
    faceCentres_ = mesh_.faceCentres();
    cellCentres_ = mesh_.cellCentres();
    cellVolumes_ = mesh_.cellVolumes();
}


void Foam::polyMeshGeometry::correct
(
    const pointField& p,
    const labelList& changedFaces
)
{
    // Update face quantities
    updateFaceCentresAndAreas(p, changedFaces);
    // Update cell quantities from face quantities
    updateCellCentresAndVols(affectedCells(mesh_, changedFaces), changedFaces);
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool report,
    const scalar orthWarn,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // for all internal and coupled faces check that the d dot S product
    // is positive

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold = ::cos(degToRad(orthWarn));

    // Calculate coupled cell centre
    pointField neiCc(mesh.nBoundaryFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;
    label nDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    for (const label facei : checkFaces)
    {
        const point& ownCc = cellCentres[own[facei]];

        if (mesh.isInternalFace(facei))
        {
            scalar dDotS = checkNonOrtho
            (
                mesh,
                report,
                severeNonorthogonalityThreshold,
                facei,
                faceAreas[facei],
                cellCentres[nei[facei]] - ownCc,

                severeNonOrth,
                errorNonOrth,
                setPtr
            );

            if (dDotS < minDDotS)
            {
                minDDotS = dDotS;
            }

            sumDDotS += dDotS;
            ++nDDotS;
        }
        else
        {
            const label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dDotS = checkNonOrtho
                (
                    mesh,
                    report,
                    severeNonorthogonalityThreshold,
                    facei,
                    faceAreas[facei],
                    neiCc[facei-mesh.nInternalFaces()] - ownCc,

                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                );

                if (dDotS < minDDotS)
                {
                    minDDotS = dDotS;
                }

                sumDDotS += dDotS;
                ++nDDotS;
            }
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        const point& ownCc = cellCentres[own[face0]];

        scalar dDotS = checkNonOrtho
        (
            mesh,
            report,
            severeNonorthogonalityThreshold,
            face0,
            faceAreas[face0],
            cellCentres[own[face1]] - ownCc,

            severeNonOrth,
            errorNonOrth,
            setPtr
         );

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }

        sumDDotS += dDotS;
        ++nDDotS;
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(nDDotS, sumOp<label>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    // Only report if there are some internal faces
    if (nDDotS > 0)
    {
        if (report && minDDotS < severeNonorthogonalityThreshold)
        {
            Info<< "Number of non-orthogonality errors: " << errorNonOrth
                << ". Number of severely non-orthogonal faces: "
                << severeNonOrth  << "." << endl;
        }
    }

    if (report)
    {
        if (nDDotS > 0)
        {
            Info<< "Mesh non-orthogonality Max: "
                << radToDeg(::acos(minDDotS))
                << " average: " << radToDeg(::acos(sumDDotS/nDDotS))
                << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in non-orthogonality detected" << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Non-orthogonality check OK.\n" << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // check whether face area vector points to the cell with higher label
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const faceList& f = mesh.faces();

    label nErrorPyrs = 0;

    for (const label facei : checkFaces)
    {
        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef
        (
            f[facei],
            cellCentres[own[facei]]
        ).mag(p);

        if (pyrVol > -minPyrVol)
        {
            ++nErrorPyrs;

            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << facei << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << f[facei] << " area: " << f[facei].mag(p)
                    << " Owner cell: " << own[facei] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[facei]].labels(f)
                    << endl;
            }
            if (setPtr)
            {
                setPtr->insert(facei);
            }
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[facei], cellCentres[nei[facei]]).mag(p);

            if (pyrVol < minPyrVol)
            {
                ++nErrorPyrs;

                if (report)
                {
                    Pout<< "bool polyMeshGeometry::checkFacePyramids("
                        << "const bool, const scalar, const pointField&"
                        << ", const labelList&, labelHashSet*): "
                        << "face " << facei << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << f[facei] << " area: " << f[facei].mag(p)
                        << " Neighbour cell: " << nei[facei] << endl
                        << "Neighbour cell vertex labels: "
                        << mesh.cells()[nei[facei]].labels(f)
                        << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        const point& ownCc = cellCentres[own[face0]];

        // Create the owner pyramid - it will have negative volume
        scalar pyrVolOwn = pyramidPointFaceRef
        (
            f[face0],
            ownCc
        ).mag(p);

        if (pyrVolOwn > -minPyrVol)
        {
            ++nErrorPyrs;

            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolOwn
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Owner cell: " << own[face0] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[face0]].labels(f)
                    << endl;
            }
            if (setPtr)
            {
                setPtr->insert(face0);
            }
        }

        // Create the neighbour pyramid - it will have positive volume
        scalar pyrVolNbr =
            pyramidPointFaceRef(f[face0], cellCentres[own[face1]]).mag(p);

        if (pyrVolNbr < minPyrVol)
        {
            ++nErrorPyrs;

            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolNbr
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Neighbour cell: " << own[face1] << endl
                    << "Neighbour cell vertex labels: "
                    << mesh.cells()[own[face1]].labels(f)
                    << endl;
            }
            if (setPtr)
            {
                setPtr->insert(face0);
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in face pyramids: faces pointing the wrong way."
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Face pyramids OK.\n" << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool report,
    const scalar minTetQuality,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // check whether decomposing each cell into tets results in
    // positive volume, non-flat tets
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nBoundaryFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        neiCc[facei - mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    label nErrorTets = 0;

    for (const label facei : checkFaces)
    {
        // Create the owner pyramid - note: exchange cell and face centre
        // to get positive volume.
        bool tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            facei,
            cellCentres[own[facei]],    // face centre
            faceCentres[facei],         // cell centre
            setPtr
        );

        if (tetError)
        {
            ++nErrorTets;
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour tets - they will have positive volume
            bool tetError = checkFaceTet
            (
                mesh,
                report,
                minTetQuality,
                p,
                facei,
                faceCentres[facei],         // face centre
                cellCentres[nei[facei]],    // cell centre
                setPtr
            );

            if (tetError)
            {
                ++nErrorTets;
            }

            if
            (
                polyMeshTetDecomposition::findSharedBasePoint
                (
                    mesh,
                    facei,
                    minTetQuality,
                    report
                ) == -1
            )
            {
                ++nErrorTets;
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                if
                (
                    polyMeshTetDecomposition::findSharedBasePoint
                    (
                        mesh,
                        facei,
                        neiCc[facei - mesh.nInternalFaces()],
                        minTetQuality,
                        report
                    ) == -1
                )
                {
                    ++nErrorTets;
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }
                }
            }
            else
            {
                if
                (
                    polyMeshTetDecomposition::findBasePoint
                    (
                        mesh,
                        facei,
                        minTetQuality,
                        report
                    ) == -1
                )
                {
                    ++nErrorTets;
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }
                }
            }
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        bool tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            face0,
            cellCentres[own[face0]],    // face centre
            faceCentres[face0],         // cell centre
            setPtr
        );

        if (tetError)
        {
            ++nErrorTets;
        }

        // Create the neighbour tets - they will have positive volume
        tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            face0,
            faceCentres[face0],         // face centre
            cellCentres[own[face1]],    // cell centre
            setPtr
        );

        if (tetError)
        {
            ++nErrorTets;
        }

        if
        (
            polyMeshTetDecomposition::findSharedBasePoint
            (
                mesh,
                face0,
                cellCentres[own[face1]],
                minTetQuality,
                report
            ) == -1
        )
        {
            ++nErrorTets;
            if (setPtr)
            {
                setPtr->insert(face0);
            }
        }
    }

    reduce(nErrorTets, sumOp<label>());

    if (nErrorTets > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in face decomposition: negative tets."
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Face tets OK.\n" << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceSkewness
(
    const bool report,
    const scalar internalSkew,
    const scalar boundarySkew,
    const polyMesh& mesh,
    const pointField& points,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the skew correction vector is more than skew times
    // larger than the face area vector

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCentres, neiCc);

    scalar maxSkew = 0;

    label nWarnSkew = 0;

    for (const label facei : checkFaces)
    {
        if (mesh.isInternalFace(facei))
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                faces,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]],
                cellCentres[nei[facei]]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                ++nWarnSkew;

                if (report)
                {
                    Pout<< "Severe skewness for face " << facei
                        << " skewness = " << skewness << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }

            maxSkew = max(maxSkew, skewness);
        }
        else if (patches[patches.whichPatch(facei)].coupled())
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                faces,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]],
                neiCc[facei-mesh.nInternalFaces()]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                ++nWarnSkew;

                if (report)
                {
                    Pout<< "Severe skewness for coupled face " << facei
                        << " skewness = " << skewness << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }

            maxSkew = max(maxSkew, skewness);
        }
        else
        {
            scalar skewness = primitiveMeshTools::boundaryFaceSkewness
            (
                faces,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]]
            );


            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > boundarySkew)
            {
                ++nWarnSkew;

                if (report)
                {
                    Pout<< "Severe skewness for boundary face " << facei
                        << " skewness = " << skewness << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }

            maxSkew = max(maxSkew, skewness);
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        const point& ownCc = cellCentres[own[face0]];
        const point& neiCc = cellCentres[own[face1]];

        scalar skewness = primitiveMeshTools::faceSkewness
        (
            faces,
            points,
            faceCentres,
            faceAreas,

            face0,
            ownCc,
            neiCc
        );

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor
        // mesh.
        if (skewness > internalSkew)
        {
            ++nWarnSkew;

            if (report)
            {
                Pout<< "Severe skewness for face " << face0
                    << " skewness = " << skewness << endl;
            }
            if (setPtr)
            {
                setPtr->insert(face0);
            }
        }

        maxSkew = max(maxSkew, skewness);
    }


    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        if (report)
        {
            WarningInFunction
                << 100*maxSkew
                << " percent.\nThis may impair the quality of the result." << nl
                << nWarnSkew << " highly skew faces detected."
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Max skewness = " << 100*maxSkew
            << " percent.  Face skewness OK.\n" << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool report,
    const scalar warnWeight,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the delta factor (0..1) is too large.

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nBoundaryFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);


    scalar minWeight = GREAT;

    label nWarnWeight = 0;

    for (const label facei : checkFaces)
    {
        const point& fc = faceCentres[facei];
        const vector& fa = faceAreas[facei];

        scalar dOwn = mag(fa & (fc-cellCentres[own[facei]]));

        if (mesh.isInternalFace(facei))
        {
            scalar dNei = mag(fa & (cellCentres[nei[facei]]-fc));
            scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

            if (weight < warnWeight)
            {
                ++nWarnWeight;

                if (report)
                {
                    Pout<< "Small weighting factor for face " << facei
                        << " weight = " << weight << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }

            minWeight = min(minWeight, weight);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dNei = mag(fa & (neiCc[facei-mesh.nInternalFaces()]-fc));
                scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

                if (weight < warnWeight)
                {
                    ++nWarnWeight;

                    if (report)
                    {
                        Pout<< "Small weighting factor for face " << facei
                            << " weight = " << weight << endl;
                    }
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }
                }

                minWeight = min(minWeight, weight);
            }
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        const point& ownCc = cellCentres[own[face0]];
        const point& fc = faceCentres[face0];
        const vector& fa = faceAreas[face0];

        scalar dOwn = mag(fa & (fc-ownCc));
        scalar dNei = mag(fa & (cellCentres[own[face1]]-fc));
        scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

        if (weight < warnWeight)
        {
            ++nWarnWeight;

            if (report)
            {
                Pout<< "Small weighting factor for face " << face0
                    << " weight = " << weight << endl;
            }
            if (setPtr)
            {
                setPtr->insert(face0);
            }
        }

        minWeight = min(minWeight, weight);
    }

    reduce(minWeight, minOp<scalar>());
    reduce(nWarnWeight, sumOp<label>());

    if (minWeight < warnWeight)
    {
        if (report)
        {
            WarningInFunction
                << minWeight << '.' << nl
                << nWarnWeight << " faces with small weights detected."
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Min weight = " << minWeight
            << ".  Weights OK.\n" << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool report,
    const scalar warnRatio,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the volume ratio between neighbouring cells is too large

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell vol
    scalarField neiVols(mesh.nBoundaryFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        neiVols[facei-mesh.nInternalFaces()] = cellVolumes[own[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiVols);


    scalar minRatio = GREAT;

    label nWarnRatio = 0;

    for (const label facei : checkFaces)
    {
        scalar ownVol = mag(cellVolumes[own[facei]]);

        scalar neiVol = -GREAT;

        if (mesh.isInternalFace(facei))
        {
            neiVol = mag(cellVolumes[nei[facei]]);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                neiVol = mag(neiVols[facei-mesh.nInternalFaces()]);
            }
        }

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

            if (ratio < warnRatio)
            {
                ++nWarnRatio;

                if (report)
                {
                    Pout<< "Small ratio for face " << facei
                        << " ratio = " << ratio << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }

            minRatio = min(minRatio, ratio);
        }
    }

    for (const labelPair& baffle : baffles)
    {
        const label face0 = baffle.first();
        const label face1 = baffle.second();

        scalar ownVol = mag(cellVolumes[own[face0]]);

        scalar neiVol = mag(cellVolumes[own[face1]]);

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

            if (ratio < warnRatio)
            {
                ++nWarnRatio;

                if (report)
                {
                    Pout<< "Small ratio for face " << face0
                        << " ratio = " << ratio << endl;
                }
                if (setPtr)
                {
                    setPtr->insert(face0);
                }
            }

            minRatio = min(minRatio, ratio);
        }
    }

    reduce(minRatio, minOp<scalar>());
    reduce(nWarnRatio, sumOp<label>());

    if (minRatio < warnRatio)
    {
        if (report)
        {
            WarningInFunction
                << minRatio << '.' << nl
                << nWarnRatio << " faces with small ratios detected."
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "Min ratio = " << minRatio
            << ".  Ratios OK.\n" << endl;
    }

    return false;
}


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorInFunction
            << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(degToRad(maxDeg));

    const faceList& fcs = mesh.faces();

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorFacei = -1;

    for (const label facei : checkFaces)
    {
        const face& f = fcs[facei];

        const vector faceNormal = normalised(faceAreas[facei]);

        // Normalized vector from f[size-1] to f[0];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        forAll(f, fp0)
        {
            // Normalized vector between two consecutive points
            vector e10(p[f.nextLabel(fp0)] - p[f.thisLabel(fp0)]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormal) < SMALL)
                    {
                        if (facei != errorFacei)
                        {
                            // Count only one error per face.
                            errorFacei = facei;
                            ++nConcave;
                        }

                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (report)
    {
        if (maxEdgeSin > SMALL)
        {
            scalar maxConcaveDegr =
                radToDeg(Foam::asin(Foam::min(1.0, maxEdgeSin)));

            Info<< "There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees.\n" << endl;
        }
        else
        {
            Info<< "All angles in faces are convex or less than "  << maxDeg
                << " degrees concave.\n" << endl;
        }
    }

    if (nConcave > 0)
    {
        if (report)
        {
            WarningInFunction
                << nConcave  << " face points with severe concave angle (> "
                << maxDeg << " deg) found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


// Check twist of faces. Is calculated as the difference between normals of
// individual triangles and the cell-cell centre edge.
bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }


    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nBoundaryFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    for (const label facei : checkFaces)
    {
        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            vector nf(Zero);

            if (mesh.isInternalFace(facei))
            {
                nf =
                    normalised
                    (
                        cellCentres[nei[facei]]
                      - cellCentres[own[facei]]
                    );
            }
            else if (patches[patches.whichPatch(facei)].coupled())
            {
                nf =
                    normalised
                    (
                        neiCc[facei-mesh.nInternalFaces()]
                      - cellCentres[own[facei]]
                    );
            }
            else
            {
                nf =
                    normalised
                    (
                        faceCentres[facei]
                      - cellCentres[own[facei]]
                    );
            }

            if (nf != vector::zero)
            {
                const point& fc = faceCentres[facei];

                forAll(f, fpI)
                {
                    vector triArea
                    (
                        triPointRef::areaNormal
                        (
                            p[f[fpI]],
                            p[f.nextLabel(fpI)],
                            fc
                        )
                    );

                    scalar magTri = mag(triArea);

                    if (magTri > VSMALL && ((nf & triArea/magTri) < minTwist))
                    {
                        ++nWarped;

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        break;
                    }
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between triangle normal and "
                << "face normal < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


// Like checkFaceTwist but compares normals of consecutive triangles.
bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    for (const label facei : checkFaces)
    {
        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Find starting triangle (at startFp) with non-zero area
            label startFp = -1;
            vector prevN;

            forAll(f, fp)
            {
                prevN = triPointRef::areaNormal
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                );

                scalar magTri = mag(prevN);

                if (magTri > VSMALL)
                {
                    startFp = fp;
                    prevN /= magTri;
                    break;
                }
            }

            if (startFp != -1)
            {
                label fp = startFp;

                do
                {
                    fp = f.fcIndex(fp);

                    vector triN
                    (
                        triPointRef::areaNormal
                        (
                            p[f[fp]],
                            p[f.nextLabel(fp)],
                            fc
                        )
                    );
                    scalar magTri = mag(triN);

                    if (magTri > VSMALL)
                    {
                        triN /= magTri;

                        if ((prevN & triN) < minTwist)
                        {
                            ++nWarped;

                            if (setPtr)
                            {
                                setPtr->insert(facei);
                            }

                            break;
                        }

                        prevN = triN;
                    }
                    else if (minTwist > 0)
                    {
                        ++nWarped;

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        break;
                    }
                }
                while (fp != startFp);
            }
        }
    }


    reduce(nWarped, sumOp<label>());

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between consecutive triangle normals less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between consecutive triangle normals is less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between consecutive triangle normals"
                << " < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool report,
    const scalar minFlatness,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minFlatness < -SMALL || minFlatness > 1+SMALL)
    {
        FatalErrorInFunction
            << "minFlatness should be [0..1] but is now " << minFlatness
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    for (const label facei : checkFaces)
    {
        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Sum triangle areas
            scalar sumArea = 0.0;

            forAll(f, fp)
            {
                sumArea += triPointRef
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                ).mag();
            }

            if (sumArea/mag(faceAreas[facei]) < minFlatness)
            {
                ++nWarped;
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with area of individual triangles"
                << " compared to overall area less than "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the area of individual triangles"
                << " compared to overall area is less than "
                << minFlatness << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " non-flat faces "
                << "(area of individual triangles"
                << " compared to overall area"
                << " < " << minFlatness << ") found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool report,
    const scalar minArea,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    label nZeroArea = 0;

    for (const label facei : checkFaces)
    {
        if (mag(faceAreas[facei]) < minArea)
        {
            ++nZeroArea;
            if (setPtr)
            {
                setPtr->insert(facei);
            }
        }
    }


    reduce(nZeroArea, sumOp<label>());

    if (report)
    {
        if (nZeroArea > 0)
        {
            Info<< "There are " << nZeroArea
                << " faces with area < " << minArea << '.' << nl << endl;
        }
        else
        {
            Info<< "All faces have area > " << minArea << '.' << nl << endl;
        }
    }

    if (nZeroArea > 0)
    {
        if (report)
        {
            WarningInFunction
                << nZeroArea  << " faces with area < " << minArea
                << " found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool report,
    const scalar warnDet,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
)
{
    const cellList& cells = mesh.cells();

    scalar minDet = GREAT;
    scalar sumDet = 0.0;
    label nSumDet = 0;
    label nWarnDet = 0;

    for (const label celli : affectedCells)
    {
        const cell& cFaces = cells[celli];

        tensor areaSum(Zero);
        scalar magAreaSum = 0;

        for (const label facei : cFaces)
        {
            scalar magArea = mag(faceAreas[facei]);

            magAreaSum += magArea;
            areaSum += faceAreas[facei]*(faceAreas[facei]/(magArea+VSMALL));
        }

        scalar scaledDet = det(areaSum/(magAreaSum+VSMALL))/0.037037037037037;

        minDet = min(minDet, scaledDet);
        sumDet += scaledDet;
        ++nSumDet;

        if (scaledDet < warnDet)
        {
            ++nWarnDet;
            if (setPtr)
            {
                setPtr->insert(cFaces);  // Insert all faces of the cell
            }
        }
    }

    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSumDet, sumOp<label>());
    reduce(nWarnDet, sumOp<label>());

    if (report)
    {
        if (nSumDet > 0)
        {
            Info<< "Cell determinant (1 = uniform cube) : average = "
                << sumDet / nSumDet << "  min = " << minDet << endl;
        }

        if (nWarnDet > 0)
        {
            Info<< "There are " << nWarnDet
                << " cells with determinant < " << warnDet << '.' << nl
                << endl;
        }
        else
        {
            Info<< "All faces have determinant > " << warnDet << '.' << nl
                << endl;
        }
    }

    if (nWarnDet > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarnDet << " cells with determinant < " << warnDet
                << " found.\n"
                << endl;
        }

        return true;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkEdgeLength
(
    const bool report,
    const scalar minEdgeLength,
    const polyMesh& mesh,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    const scalar reportLenSqr(Foam::sqr(minEdgeLength));

    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();

    scalar minLenSqr = sqr(GREAT);
    scalar maxLenSqr = -sqr(GREAT);

    label nSmall = 0;

    for (const label facei : checkFaces)
    {
        const face& f = faces[facei];

        forAll(f, fp)
        {
            label fp1 = f.fcIndex(fp);

            scalar magSqrE = magSqr(points[f[fp]] - points[f[fp1]]);

            if (setPtr && magSqrE < reportLenSqr)
            {
                if (setPtr->insert(f[fp]) || setPtr->insert(f[fp1]))
                {
                    nSmall++;
                }
            }

            minLenSqr = min(minLenSqr, magSqrE);
            maxLenSqr = max(maxLenSqr, magSqrE);
        }
    }

    reduce(minLenSqr, minOp<scalar>());
    reduce(maxLenSqr, maxOp<scalar>());

    reduce(nSmall, sumOp<label>());

    if (nSmall > 0)
    {
        if (report)
        {
            Info<< "   *Edges too small, min/max edge length = "
                << sqrt(minLenSqr) << " " << sqrt(maxLenSqr)
                << ", number too small: " << nSmall << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "    Min/max edge length = "
            << sqrt(minLenSqr) << " " << sqrt(maxLenSqr) << " OK." << endl;
    }

    return false;
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool report,
    const scalar orthWarn,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceDotProduct
    (
        report,
        orthWarn,
        mesh_,
        cellCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFacePyramids
    (
        report,
        minPyrVol,
        mesh_,
        cellCentres_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool report,
    const scalar minTetQuality,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceTets
    (
        report,
        minTetQuality,
        mesh_,
        cellCentres_,
        faceCentres_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


//bool Foam::polyMeshGeometry::checkFaceSkewness
//(
//    const bool report,
//    const scalar internalSkew,
//    const scalar boundarySkew,
//    const labelList& checkFaces,
//    const List<labelPair>& baffles,
//    labelHashSet* setPtr
//) const
//{
//    return checkFaceSkewness
//    (
//        report,
//        internalSkew,
//        boundarySkew,
//        mesh_,
//        cellCentres_,
//        faceCentres_,
//        faceAreas_,
//        checkFaces,
//        baffles,
//        setPtr
//    );
//}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool report,
    const scalar warnWeight,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceWeights
    (
        report,
        warnWeight,
        mesh_,
        cellCentres_,
        faceCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool report,
    const scalar warnRatio,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkVolRatio
    (
        report,
        warnRatio,
        mesh_,
        cellVolumes_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceAngles
    (
        report,
        maxDeg,
        mesh_,
        faceAreas_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool report,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceTwist
    (
        report,
        minTwist,
        mesh_,
        cellCentres_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool report,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkTriangleTwist
    (
        report,
        minTwist,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool report,
    const scalar minFlatness,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceFlatness
    (
        report,
        minFlatness,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool report,
    const scalar minArea,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceArea
    (
        report,
        minArea,
        mesh_,
        faceAreas_,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool report,
    const scalar warnDet,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
) const
{
    return checkCellDeterminant
    (
        report,
        warnDet,
        mesh_,
        faceAreas_,
        checkFaces,
        affectedCells,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkEdgeLength
(
    const bool report,
    const scalar minEdgeLength,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkEdgeLength
    (
        report,
        minEdgeLength,
        mesh_,
        checkFaces,
        setPtr
    );
}


// ************************************************************************* //
