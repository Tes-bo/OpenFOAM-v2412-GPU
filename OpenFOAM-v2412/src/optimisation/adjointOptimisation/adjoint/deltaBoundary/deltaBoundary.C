/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "deltaBoundary.H"
#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tensor deltaBoundary::tensorCrossVector(const tensor& T, const vector& v)
{
    // The correct approach when T is not a diagonal tensor
    tensor res(Zero);
    vector vec1(T.xx(), T.yx(), T.zx());
    vector res1(vec1 ^ v);
    res.xx() = res1.x(); res.yx() = res1.y(); res.zx() = res1.z();

    vector vec2(T.xy(), T.yy(), T.zy());
    vector res2(vec2 ^ v);
    res.xy() = res2.x(); res.yy() = res2.y(); res.zy() = res2.z();

    vector vec3(T.xz(), T.yz(), T.zz());
    vector res3(vec3 ^ v);
    res.xz() = res3.x(); res.yz() = res3.y(); res.zz() = res3.z();

    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

deltaBoundary::deltaBoundary(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vectorField deltaBoundary::makeFaceCentresAndAreas_d
(
    const pointField& p,
    const pointField& p_d
)
{
    vector fCtrs_d(Zero);
    vector fAreas_d(Zero);
    vector unitVector_d(Zero);

    // Container field to return results
    vectorField deltaVecs(3, Zero);

    label nPoints = p.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        //fCtrs[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
        vector fAreas = 0.5*((p[1] - p[0])^(p[2] - p[0]));

        fCtrs_d  = (1.0/3.0)*(p_d[0] + p_d[1] + p_d[2]);
        fAreas_d =
            0.5*((p_d[1] - p_d[0])^(p[2] - p[0]))
          + 0.5*((p[1] - p[0])^(p_d[2] - p_d[0]));
        scalar ds = mag(fAreas);
        unitVector_d = fAreas_d/ds - (fAreas*(fAreas&fAreas_d))/ds/ds/ds;

        deltaVecs[0] = fCtrs_d;
        deltaVecs[1] = fAreas_d;
        deltaVecs[2] = unitVector_d;
    }
    else
    {
        vector sumN(Zero);
        vector sumN_d(Zero);
        scalar sumA(0);
        scalar sumA_d(0);
        vector sumAc(Zero);
        vector sumAc_d(Zero);

        point fCentre = p[0];
        point fCentre_d = p_d[0];
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += p[pi];
            fCentre_d += p_d[pi];
        }

        fCentre /= nPoints;
        fCentre_d /= nPoints;

        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = p[(pi + 1) % nPoints];
            const point& nextPoint_d = p_d[(pi + 1) % nPoints];

            vector c = p[pi] + nextPoint + fCentre;
            vector c_d = p_d[pi] + nextPoint_d + fCentre_d;

            vector n = (nextPoint - p[pi])^(fCentre - p[pi]);
            vector n_d =
                ((nextPoint_d - p_d[pi])^(fCentre - p[pi]))
              + ((nextPoint - p[pi])^(fCentre_d - p_d[pi]));

            scalar a = mag(n);
            if (a < ROOTVSMALL)
            {
                // This shouldn't happen in general.
                // Manually zero contribution from zero area face for now
                WarningInFunction
                    << "Zero area face sub triangle found " << endl
                    << p[pi] << " " << nextPoint <<  " " << fCentre << endl
                    << "Neglecting contributions of this element " << endl;
            }
            else
            {
                scalar a_d = (n&n_d)/mag(n);

                sumN += n;
                sumN_d += n_d;

                sumA += a;
                sumA_d += a_d;

                sumAc += a*c;
                sumAc_d += a_d*c + a*c_d;
            }
        }

        // fCtrs[facei] = (1.0/3.0)*sumAc/(sumA + VSMALL);
        vector fAreas = 0.5*sumN;
        fCtrs_d = (1.0/3.0)*(sumAc_d*sumA - sumAc*sumA_d)/sumA/sumA;
        fAreas_d = 0.5*sumN_d;
        scalar ds = mag(fAreas);
        unitVector_d = fAreas_d/ds - (fAreas*(fAreas&fAreas_d))/ds/ds/ds;

        deltaVecs[0] = fCtrs_d;
        deltaVecs[1] = fAreas_d;
        deltaVecs[2] = unitVector_d;
    }

    return deltaVecs;
}


tensorField deltaBoundary::makeFaceCentresAndAreas_d
(
    const pointField& p,
    const tensorField& p_d
)
{
    label nPoints = p.size();
    tensor fCtrs_d(Zero);
    tensor fAreas_d(Zero);
    tensor unitVector_d(Zero);

    // Container field to return results
    tensorField deltaVecs(3, Zero);

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        vector fAreas = 0.5*((p[1] - p[0])^(p[2] - p[0]));

        fCtrs_d  = (1.0/3.0)*(p_d[0] + p_d[1] + p_d[2]);
        fAreas_d =
            0.5*tensorCrossVector(p_d[1] - p_d[0], p[2] - p[0])
            //minus sign since it is vector ^ tensor
          - 0.5*tensorCrossVector(p_d[2] - p_d[0], p[1] - p[0]);
        scalar ds = mag(fAreas);
        unitVector_d = fAreas_d/ds - (fAreas*(fAreas&fAreas_d))/ds/ds/ds;

        deltaVecs[0] = fCtrs_d;
        deltaVecs[1] = fAreas_d;
        deltaVecs[2] = unitVector_d;
    }
    else
    {
        vector sumN(Zero);
        tensor sumN_d(Zero);
        scalar sumA(0);
        vector sumA_d(Zero);
        vector sumAc(Zero);
        tensor sumAc_d(Zero);

        point fCentre = p[0];
        tensor fCentre_d = p_d[0];
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += p[pi];
            fCentre_d += p_d[pi];
        }

        fCentre /= nPoints;
        fCentre_d /= nPoints;

        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = p[(pi + 1) % nPoints];
            const tensor& nextPoint_d = p_d[(pi + 1) % nPoints];

            vector c = p[pi] + nextPoint + fCentre;
            tensor c_d = p_d[pi] + nextPoint_d + fCentre_d;

            vector n = (nextPoint - p[pi])^(fCentre - p[pi]);
            tensor n_d =
                tensorCrossVector(nextPoint_d - p_d[pi], fCentre - p[pi])
                //minus sign since it is vector ^ tensor
              - tensorCrossVector(fCentre_d - p_d[pi], nextPoint - p[pi]);

            scalar a = mag(n);
            if (a < ROOTVSMALL)
            {
                // This shouldn't happen in general.
                // Manually zero contribution from zero area face for now
                WarningInFunction
                    << "Zero area face sub triangle found " << nl
                    << p[pi] << " " << nextPoint <<  " " << fCentre << nl
                    << "Neglecting contributions of this element " << endl;
            }
            else
            {
                vector a_d = (n & n_d)/a;

                sumN += n;
                sumN_d += n_d;

                sumA += a;
                sumA_d += a_d;

                sumAc += a*c;
                // c*a_d since we need to get the correct outer product
                sumAc_d += (c*a_d) + a*c_d;
            }
        }

        vector fAreas = 0.5*sumN;
        fCtrs_d = (1.0/3.0)*(sumAc_d/sumA - (sumAc*sumA_d)/sqr(sumA));
        fAreas_d = 0.5*sumN_d;
        scalar ds = mag(fAreas);
        unitVector_d = fAreas_d/ds - (fAreas*(fAreas&fAreas_d))/ds/ds/ds;

        deltaVecs[0] = fCtrs_d;
        deltaVecs[1] = fAreas_d;
        deltaVecs[2] = unitVector_d;
    }

    return deltaVecs;
}


tmp<tensorField> deltaBoundary::cellCenters_d(const label pointI)
{
    const labelListList& pointCells(mesh_.pointCells());
    const labelList& pointCellsI(pointCells[pointI]);
    const pointField& points(mesh_.points());
    auto tC_d = tmp<tensorField>::New(pointCellsI.size(), Zero);
    auto& C_d = tC_d.ref();

    const labelList& pointFaces(mesh_.pointFaces()[pointI]);
    tensorField Cf_d(pointFaces.size(), Zero);
    tensorField Sf_d(pointFaces.size(), Zero);

    forAll(pointFaces, pfI)
    {
        const label pointFaceI = pointFaces[pfI];
        const face& faceI = mesh_.faces()[pointFaceI];
        tensorField p_d(faceI.size(), Zero);
        forAll(faceI, pI)
        {
            if (faceI[pI] == pointI)
            {
                p_d[pI] = tensor::I;
                break;
            }
        }

        pointField facePoints(faceI.points(points));

        // Compute changes in the face
        tensorField dFace(makeFaceCentresAndAreas_d(facePoints, p_d));
        Cf_d[pfI] = dFace[0];
        Sf_d[pfI] = dFace[1];
    }

    // Face variations have now been computed. Now, compute cell contributions
    forAll(pointCellsI, pcI)
    {
        const label pointCellI = pointCellsI[pcI];
        const cell& cellI(mesh_.cells()[pointCellI]);
        vectorField fAreas(cellI.size(), Zero);
        vectorField fCtrs(cellI.size(), Zero);
        tensorField fAreas_d(cellI.size(), Zero);
        tensorField fCtrs_d(cellI.size(), Zero);
        forAll(cellI, fI)
        {
            const label globalFaceI = cellI[fI];

            // Assign values to faceAreas and faceCtrs
            if (globalFaceI < mesh_.nInternalFaces())
            {
                fAreas[fI] = mesh_.Sf()[globalFaceI];
                fCtrs[fI] = mesh_.Cf()[globalFaceI];
            }
            else
            {
                const label whichPatch =
                    mesh_.boundaryMesh().whichPatch(globalFaceI);
                const fvPatch& patch = mesh_.boundary()[whichPatch];
                const label patchStart = patch.patch().start();
                const label localFace = globalFaceI - patchStart;
                fAreas[fI] = patch.Sf()[localFace];
                fCtrs[fI] = patch.Cf()[localFace];
            }

            // Assign values to differentiated face areas and centres
            forAll(pointFaces, pfI)
            {
                if (pointFaces[pfI] == globalFaceI)
                {
                    fAreas_d[fI] = Sf_d[pfI];
                    fCtrs_d[fI] = Cf_d[pfI];
                }
            }
        }
        C_d[pcI] = makeCellCentres_d(fAreas, fCtrs, fAreas_d, fCtrs_d);
    }

    return tC_d;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
