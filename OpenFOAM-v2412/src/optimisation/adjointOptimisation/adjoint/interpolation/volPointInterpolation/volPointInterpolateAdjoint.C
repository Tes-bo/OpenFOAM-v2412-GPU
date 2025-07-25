/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "volPointInterpolationAdjoint.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "coupledPointPatchField.H"
#include "pointConstraints.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::volPointInterpolationAdjoint::pushUntransformedData
(
    List<Type>& pointData
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::volPointInterpolationAdjoint::flatBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvBoundaryMesh& bm = mesh.boundary();

    auto tboundaryVals = tmp<Field<Type>>::New(mesh.nBoundaryFaces());
    auto& boundaryVals = tboundaryVals.ref();

    forAll(vf.boundaryField(), patchi)
    {
        label bFacei = bm[patchi].patch().start() - mesh.nInternalFaces();

        if
        (
           !isA<emptyFvPatch>(bm[patchi])
        && !vf.boundaryField()[patchi].coupled()
        )
        {
            SubList<Type>
            (
                boundaryVals,
                vf.boundaryField()[patchi].size(),
                bFacei
            ) = vf.boundaryField()[patchi];
        }
        else
        {
            const polyPatch& pp = bm[patchi].patch();

            forAll(pp, i)
            {
                boundaryVals[bFacei++] = Zero;
            }
        }
    }

    return tboundaryVals;
}


template<class Type>
void Foam::volPointInterpolationAdjoint::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::addSeparated" << endl;
    }

    auto& pfi = pf.ref();
    auto& pfbf = pf.boundaryFieldRef();

    const label startOfRequests = UPstream::nRequests();

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).initSwapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }

    // Wait for outstanding requests
    UPstream::waitRequests(startOfRequests);

    forAll(pfbf, patchi)
    {
        if (pfbf[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type>>
                (pfbf[patchi]).swapAddSeparated
                (
                    Pstream::commsTypes::nonBlocking,
                    pfi
                );
        }
    }
}


template<class Type>
void Foam::volPointInterpolationAdjoint::interpolateSensitivitiesField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    typename GeometricField<Type, fvPatchField, volMesh>::Boundary& vf,
    const labelHashSet& patchIDs
) const
{
    const Field<Type>& pfi = pf.primitiveField();

    // Get face data in flat list
    const fvMesh& Mesh = mesh();
    const fvBoundaryMesh& bm = Mesh.boundary();

    auto tboundaryVals = tmp<Field<Type>>::New(Mesh.nBoundaryFaces(), Zero);
    auto& boundaryVals = tboundaryVals.ref();

    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const primitivePatch& boundary = boundaryPtr_();
    const labelList& mp = boundary.meshPoints();

    forAll(mp, i)
    {
        label pointi = mp[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];
            const Type& val = pfi[pointi];
            // Face-to-point weights should, in general, have half the weight
            // of what they actually do in volPointInterpolation since, in
            // a complete case, a face laying on the opposite side of the
            // symmetry plane would also contribute to a point laying on
            // the symmetry plane.
            // For face-to-point interpolation this is not a problem, but for
            // the adjoint point-to-face interpolation, the correct value of
            // the weight should be taken into consideration
            scalar mod(isSymmetryPoint_[pointi] ? 0.5 : 1);

            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    boundaryVals[pFaces[j]] += mod*pWeights[j]*val;
                }
            }
        }
    }

    // Transfer values to face-based sensitivity field
    for (const label patchi : patchIDs)
    {
        label bFacei = bm[patchi].patch().start() - Mesh.nInternalFaces();
        if (!isA<emptyFvPatch>(bm[patchi]) && !vf[patchi].coupled())
        {
            vf[patchi] =
                SubList<Type>
                (
                    boundaryVals,
                    vf[patchi].size(),
                    bFacei
                );
        }
    }
}


template<class Type>
void Foam::volPointInterpolationAdjoint::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    const primitivePatch& boundary = boundaryPtr_();

    Field<Type>& pfi = pf.primitiveFieldRef();

    // Get face data in flat list
    tmp<Field<Type>> tboundaryVals(flatBoundaryField(vf));
    const Field<Type>& boundaryVals = tboundaryVals();


    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& mp = boundary.meshPoints();

    forAll(mp, i)
    {
        label pointi = mp[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];

            Type& val = pfi[pointi];

            val = Zero;
            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    scalar mod(1);
                    if (isSymmetryPoint_[pointi])
                    {
                        const label globalFaceI =
                            mesh().nInternalFaces() + pFaces[j];
                        const label facePatchID =
                            mesh().boundaryMesh().whichPatch(globalFaceI);
                        const polyPatch& pp = mesh().boundaryMesh()[facePatchID];
                        if
                        (
                             isA<symmetryPolyPatch>(pp)
                          || isA<symmetryPlanePolyPatch>(pp)
                        )
                        {
                            mod = 0;
                        }
                      //else
                      //{
                      //    mod = 2;
                      //}
                    }
                    val += mod*pWeights[j]*boundaryVals[pFaces[j]];
                }
            }
        }
    }

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh(), pfi, plusEqOp<Type>());

    // And add separated contributions
    addSeparated(pf);

    // Optionally normalise
    /*
    if (normalisationPtr_)
    {
        const scalarField& normalisation = normalisationPtr_();
        forAll(mp, i)
        {
            pfi[mp[i]] *= normalisation[i];
        }
    }
    */


    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves.
    pushUntransformedData(pfi);

    // Apply displacement constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrain(pf, false);
}


// ************************************************************************* //
