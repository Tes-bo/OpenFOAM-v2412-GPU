/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "isoSurfacePoint.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::VolumeSliceField<Type>>
Foam::isoSurfacePoint::adaptPatchFields
(
    const VolumeField<Type>& fld
) const
{
    auto tslice = tmp<VolumeSliceField<Type>>::New
    (
        IOobject
        (
            fld.name(),
            fld.instance(),
            fld.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        fld,        // internal field
        true        // preserveCouples
    );
    auto& sliceFld = tslice.ref();

    const fvMesh& mesh = fld.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    auto& sliceFldBf = sliceFld.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isA<emptyPolyPatch>(pp)
         && pp.size() != sliceFldBf[patchi].size()
        )
        {
            // Clear old value. Cannot resize it since is a slice.
            // Set new value we can change

            sliceFldBf.release(patchi);
            sliceFldBf.set
            (
                patchi,
                new calculatedFvPatchField<Type>
                (
                    mesh.boundary()[patchi],
                    sliceFld
                )
            );

            // Note: cannot use patchInternalField since uses emptyFvPatch::size
            // Do our own internalField instead.
            const labelUList& faceCells =
                mesh.boundary()[patchi].patch().faceCells();

            Field<Type>& pfld = sliceFldBf[patchi];
            pfld.setSize(faceCells.size());
            forAll(faceCells, i)
            {
                pfld[i] = sliceFld[faceCells[i]];
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            // Already has interpolate as value
        }
        else if (isA<processorPolyPatch>(pp))
        {
            fvPatchField<Type>& pfld = const_cast<fvPatchField<Type>&>
            (
                sliceFldBf[patchi]
            );

            tmp<Field<Type>> f = lerp
            (
                pfld.patchNeighbourField(),
                pfld.patchInternalField(),
                mesh.weights().boundaryField()[patchi]
            );

            bitSet isCollocated
            (
                collocatedFaces(refCast<const processorPolyPatch>(pp))
            );

            forAll(isCollocated, i)
            {
                if (!isCollocated[i])
                {
                    pfld[i] = f()[i];
                }
            }
        }
    }
    return tslice;
}


template<class Type>
Type Foam::isoSurfacePoint::generatePoint
(
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1
) const
{
    const scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        const scalar s = (iso_-s0)/d;

        if (hasSnap1 && s >= 0.5 && s <= 1)
        {
            return snapP1;
        }
        else if (hasSnap0 && s >= 0.0 && s <= 0.5)
        {
            return snapP0;
        }
        else
        {
            return s*p1 + (1.0-s)*p0;
        }
    }
    else
    {
        constexpr scalar s = 0.4999;

        return s*p1 + (1.0-s)*p0;
    }
}


template<class Type>
void Foam::isoSurfacePoint::generateTriPoints
(
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1,

    const scalar s2,
    const Type& p2,
    const bool hasSnap2,
    const Type& snapP2,

    const scalar s3,
    const Type& p3,
    const bool hasSnap3,
    const Type& snapP3,

    DynamicList<Type>& pts
) const
{
    // Note: cannot use simpler isoSurfaceCell::generateTriPoints since
    //       the need here to sometimes pass in remote 'snappedPoints'

    int triIndex = 0;
    if (s0 < iso_)
    {
        triIndex |= 1;
    }
    if (s1 < iso_)
    {
        triIndex |= 2;
    }
    if (s2 < iso_)
    {
        triIndex |= 4;
    }
    if (s3 < iso_)
    {
        triIndex |= 8;
    }

    // Form the vertices of the triangles for each case
    switch (triIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x01:
        case 0x0E:
        {
            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1)
            );
            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
            );
            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
            if (triIndex == 0x0E)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s0,p0,hasSnap0,snapP0)
            );
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
            );
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            if (triIndex == 0x0D)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            Type p0p2 =
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2);
            Type p1p3 =
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3);

            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
            pts.append(p1p3);
            pts.append(p0p2);

            pts.append(p1p3);
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            pts.append(p0p2);

            if (triIndex == 0x0C)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-5], pts[sz-4]);
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            pts.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s0,p0,hasSnap0,snapP0)
            );
            pts.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s1,p1,hasSnap1,snapP1)
            );
            pts.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3)
            );

            if (triIndex == 0x0B)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            Type p0p1 =
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
            Type p2p3 =
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

            pts.append(p0p1);
            pts.append(p2p3);
            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );

            pts.append(p0p1);
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            pts.append(p2p3);

            if (triIndex == 0x0A)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-5], pts[sz-4]);
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            Type p0p1 =
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
            Type p2p3 =
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

            pts.append(p0p1);
            pts.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
            );
            pts.append(p2p3);

            pts.append(p0p1);
            pts.append(p2p3);
            pts.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
            );

            if (triIndex == 0x09)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-5], pts[sz-4]);
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            pts.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s0,p0,hasSnap0,snapP0)
            );
            pts.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s2,p2,hasSnap2,snapP2)
            );
            pts.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s1,p1,hasSnap1,snapP1)
            );

            if (triIndex == 0x07)
            {
                // Flip normals
                const label sz = pts.size();
                std::swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;
    }
}


template<class Type>
Foam::label Foam::isoSurfacePoint::generateFaceTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const VolumeField<Type>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,
    const label facei,

    const scalar neiVal,
    const Type& neiPt,
    const bool hasNeiSnap,
    const Type& neiSnapPt,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    const label own = mesh_.faceOwner()[facei];

    label oldNPoints = triPoints.size();

    const face& f = mesh_.faces()[facei];

    forAll(f, fp)
    {
        label pointi = f[fp];
        label nextPointi = f[f.fcIndex(fp)];

        generateTriPoints
        (
            pVals[pointi],
            pCoords[pointi],
            snappedPoint[pointi] != -1,
            (
                snappedPoint[pointi] != -1
              ? snappedPoints[snappedPoint[pointi]]
              : Type(Zero)
            ),

            pVals[nextPointi],
            pCoords[nextPointi],
            snappedPoint[nextPointi] != -1,
            (
                snappedPoint[nextPointi] != -1
              ? snappedPoints[snappedPoint[nextPointi]]
              : Type(Zero)
            ),

            cVals[own],
            cCoords[own],
            snappedCc[own] != -1,
            (
                snappedCc[own] != -1
              ? snappedPoints[snappedCc[own]]
              : Type(Zero)
            ),

            neiVal,
            neiPt,
            hasNeiSnap,
            neiSnapPt,

            triPoints
        );
    }

    // Every three triPoints is a triangle
    label nTris = (triPoints.size()-oldNPoints)/3;
    for (label i = 0; i < nTris; i++)
    {
        triMeshCells.append(own);
    }

    return nTris;
}


template<class Type>
void Foam::isoSurfacePoint::generateTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const VolumeField<Type>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    if
    (
        (cVals.size() != mesh_.nCells())
     || (pVals.size() != mesh_.nPoints())
     || (cCoords.size() != mesh_.nCells())
     || (pCoords.size() != mesh_.nPoints())
     || (snappedCc.size() != mesh_.nCells())
     || (snappedPoint.size() != mesh_.nPoints())
    )
    {
        FatalErrorInFunction
            << "Incorrect size." << endl
            << "mesh: nCells:" << mesh_.nCells()
            << " points:" << mesh_.nPoints() << endl
            << "cVals:" << cVals.size() << endl
            << "cCoords:" << cCoords.size() << endl
            << "snappedCc:" << snappedCc.size() << endl
            << "pVals:" << pVals.size() << endl
            << "pCoords:" << pCoords.size() << endl
            << "snappedPoint:" << snappedPoint.size() << endl
            << abort(FatalError);
    }


    // Generate triangle points

    triPoints.clear();
    triMeshCells.clear();

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if ((faceCutType_[facei] & cutType::ANYCUT) != 0)
        {
            generateFaceTriPoints
            (
                cVals,
                pVals,

                cCoords,
                pCoords,

                snappedPoints,
                snappedCc,
                snappedPoint,
                facei,

                cVals[nei[facei]],
                cCoords[nei[facei]],
                snappedCc[nei[facei]] != -1,
                (
                    snappedCc[nei[facei]] != -1
                  ? snappedPoints[snappedCc[nei[facei]]]
                  : Type(Zero)
                ),

                triPoints,
                triMeshCells
            );
        }
    }


    // Determine neighbouring snap status
    boolList neiSnapped(mesh_.nBoundaryFaces(), false);
    List<Type> neiSnappedPoint(neiSnapped.size(), Type(Zero));
    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                label bFacei = facei-mesh_.nInternalFaces();
                label snappedIndex = snappedCc[own[facei]];

                if (snappedIndex != -1)
                {
                    neiSnapped[bFacei] = true;
                    neiSnappedPoint[bFacei] = snappedPoints[snappedIndex];
                }
                facei++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiSnapped);
    syncTools::swapBoundaryFaceList(mesh_, neiSnappedPoint);


    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& cpp =
                refCast<const processorPolyPatch>(pp);

            bitSet isCollocated(collocatedFaces(cpp));

            forAll(isCollocated, i)
            {
                const label facei = pp.start()+i;

                if ((faceCutType_[facei] & cutType::ANYCUT) != 0)
                {
                    if (isCollocated[i])
                    {
                        generateFaceTriPoints
                        (
                            cVals,
                            pVals,

                            cCoords,
                            pCoords,

                            snappedPoints,
                            snappedCc,
                            snappedPoint,
                            facei,

                            cVals.boundaryField()[patchi][i],
                            cCoords.boundaryField()[patchi][i],
                            neiSnapped[facei-mesh_.nInternalFaces()],
                            neiSnappedPoint[facei-mesh_.nInternalFaces()],

                            triPoints,
                            triMeshCells
                        );
                    }
                    else
                    {
                        generateFaceTriPoints
                        (
                            cVals,
                            pVals,

                            cCoords,
                            pCoords,

                            snappedPoints,
                            snappedCc,
                            snappedPoint,
                            facei,

                            cVals.boundaryField()[patchi][i],
                            cCoords.boundaryField()[patchi][i],
                            false,
                            Type(Zero),

                            triPoints,
                            triMeshCells
                        );
                    }
                }
            }
        }
        else
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if ((faceCutType_[facei] & cutType::ANYCUT) != 0)
                {
                    generateFaceTriPoints
                    (
                        cVals,
                        pVals,

                        cCoords,
                        pCoords,

                        snappedPoints,
                        snappedCc,
                        snappedPoint,
                        facei,

                        cVals.boundaryField()[patchi][i],
                        cCoords.boundaryField()[patchi][i],
                        false,  // fc not snapped
                        Type(Zero),

                        triPoints,
                        triMeshCells
                    );
                }
                facei++;
            }
        }
    }

    triPoints.shrink();
    triMeshCells.shrink();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfacePoint::interpolate
(
    const label nPoints,
    const labelList& triPointMergeMap,
    const labelList& interpolatedPoints,
    const List<FixedList<label, 3>>& interpolatedOldPoints,
    const List<FixedList<scalar, 3>>& interpolationWeights,
    const DynamicList<Type>& unmergedValues
)
{
    // One value per point
    auto tvalues = tmp<Field<Type>>::New(nPoints, Type(Zero));
    auto& values = tvalues.ref();


    // Pass1: unweighted average of merged point values
    {
        labelList nValues(values.size(), Zero);

        forAll(unmergedValues, i)
        {
            label mergedPointi = triPointMergeMap[i];

            if (mergedPointi >= 0)
            {
                values[mergedPointi] += unmergedValues[i];
                nValues[mergedPointi]++;
            }
        }

        forAll(values, i)
        {
            if (nValues[i] > 0)
            {
                values[i] /= scalar(nValues[i]);
            }
        }
    }


    // Pass2: weighted average for remaining values (from clipped triangles)

    forAll(interpolatedPoints, i)
    {
        label pointi = interpolatedPoints[i];
        const FixedList<label, 3>& oldPoints = interpolatedOldPoints[i];
        const FixedList<scalar, 3>& w = interpolationWeights[i];

        // Note: zeroing should not be necessary if interpolation only done
        //       for newly introduced points (i.e. not in triPointMergeMap)
        values[pointi] = Type(Zero);
        forAll(oldPoints, j)
        {
            values[pointi] = w[j]*unmergedValues[oldPoints[j]];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfacePoint::interpolateTemplate
(
    const VolumeField<Type>& cCoords,
    const Field<Type>& pCoords
) const
{
    // Recalculate boundary values
    tmp<VolumeSliceField<Type>> c2(adaptPatchFields(cCoords));


    DynamicList<Type> triPoints(3*nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);

    generateTriPoints
    (
        cValsPtr_(),
        pVals_,

        c2(),
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells
    );

    return interpolate
    (
        this->points().size(),
        triPointMergeMap_,
        interpolatedPoints_,
        interpolatedOldPoints_,
        interpolationWeights_,
        triPoints
    );
}


// ************************************************************************* //
