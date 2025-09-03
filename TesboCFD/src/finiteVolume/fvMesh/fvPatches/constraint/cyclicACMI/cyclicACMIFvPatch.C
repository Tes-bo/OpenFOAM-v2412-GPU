/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "cyclicACMIFvPatch.H"
#include "fvMesh.H"
#include "transform.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::cyclicACMIFvPatch::updateAreas() const
{
    // Give AMI chance to update itself
    bool updated = cyclicACMIPolyPatch_.updateAreas();

    if (!cyclicACMIPolyPatch_.owner())
    {
        return updated;
    }

    if (updated || !cyclicACMIPolyPatch_.upToDate(areaTime_))
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::updateAreas() : updating fv areas for "
                << name() << " and " << this->nonOverlapPatch().name()
                << endl;
        }

        const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
        const cyclicACMIFvPatch& nbrACMI = neighbPatch();
        const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();

        resetPatchAreas(*this);
        resetPatchAreas(nonOverlapPatch);
        resetPatchAreas(nbrACMI);
        resetPatchAreas(nbrNonOverlapPatch);

        updated = true;

        // Mark my data to be up to date with ACMI polyPatch level
        cyclicACMIPolyPatch_.setUpToDate(areaTime_);
    }
    return updated;
}


void Foam::cyclicACMIFvPatch::resetPatchAreas(const fvPatch& fvp) const
{
    const_cast<vectorField&>(fvp.Sf()) = fvp.patch().faceAreas();
    const_cast<vectorField&>(fvp.Cf()) = fvp.patch().faceCentres();
    const_cast<scalarField&>(fvp.magSf()) = mag(fvp.patch().faceAreas());

    DebugPout
        << fvp.patch().name() << " area:" << sum(fvp.magSf()) << endl;
}


void Foam::cyclicACMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        const scalarField deltas(nf() & coupledFvPatch::delta());

        // These deltas are of the cyclic part alone - they are
        // not affected by the amount of overlap with the nonOverlapPatch
        scalarField nbrDeltas
        (
            interpolate
            (
                nbrPatch.nf() & nbrPatch.coupledFvPatch::delta()
            )
        );

        const scalar tol = cyclicACMIPolyPatch::tolerance();


        forAll(deltas, facei)
        {
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            if (dni < tol)
            {
                // Avoid zero weights on disconnected faces. This value
                // will be weighted with the (zero) face area so will not
                // influence calculations.
                w[facei] = 1.0;
            }
            else
            {
                w[facei] = dni/(di + dni);
            }
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIFvPatch::cyclicACMIFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    coupledFvPatch(patch, bm),
    cyclicACMILduInterface(),
    cyclicACMIPolyPatch_(refCast<const cyclicACMIPolyPatch>(patch)),
    areaTime_
    (
        IOobject
        (
            "areaTime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    )
{
    areaTime_.eventNo() = -1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::cyclicACMIFvPatch::newInternalProcFaces
// (
//     label& newFaces,
//     label& newProcFaces
// ) const
// {
//     const List<labelList>& addSourceFaces =
//         cyclicACMIPolyPatch_.AMI().srcAddress();
//
//     const scalarField& fMask = cyclicACMIPolyPatch_.srcMask();
//
//     // Add new faces as many weights for AMI
//     forAll (addSourceFaces, faceI)
//     {
//         if (fMask[faceI] > cyclicACMIPolyPatch_.tolerance_)
//         {
//             const labelList& nbrFaceIs = addSourceFaces[faceI];
//
//             forAll (nbrFaceIs, j)
//             {
//                 label nbrFaceI = nbrFaceIs[j];
//
//                 if (nbrFaceI < neighbPatch().size())
//                 {
//                     // local faces
//                     newFaces++;
//                 }
//                 else
//                 {
//                     // Proc faces
//                     newProcFaces++;
//                 }
//             }
//         }
//     }
// }


// Foam::refPtr<Foam::labelListList>
// Foam::cyclicACMIFvPatch::mapCollocatedFaces() const
// {
//     const scalarField& fMask = cyclicACMIPolyPatch_.srcMask();
//     const labelListList& srcFaces = cyclicACMIPolyPatch_.AMI().srcAddress();
//     labelListList dOverFaces;
//
//     dOverFaces.setSize(srcFaces.size());
//     forAll (dOverFaces, faceI)
//     {
//         if (fMask[faceI] > cyclicACMIPolyPatch_.tolerance_)
//         {
//             dOverFaces[faceI].setSize(srcFaces[faceI].size());
//
//             forAll (dOverFaces[faceI], subFaceI)
//             {
//                 dOverFaces[faceI][subFaceI] = srcFaces[faceI][subFaceI];
//             }
//         }
//     }
//     return refPtr<labelListList>(new labelListList(dOverFaces));
// }


bool Foam::cyclicACMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


Foam::tmp<Foam::vectorField> Foam::cyclicACMIFvPatch::delta() const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();

        const vectorField patchD(coupledFvPatch::delta());

        vectorField nbrPatchD(interpolate(nbrPatch.coupledFvPatch::delta()));

        auto tpdv = tmp<vectorField>::New(patchD.size());
        vectorField& pdv = tpdv.ref();

        // do the transformation if necessary
        if (parallel())
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }
        else
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform(forwardT()[0], dni);
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, faceCells, tpfld.ref());
    return tpfld;
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


void Foam::cyclicACMIFvPatch::movePoints()
{
    // Update local and higher level areas
    const bool updated = updateAreas();

    // If anything changed update the mesh flux
    if (cyclicACMIPolyPatch_.owner() && updated)
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::movePoints() : areas updated for "
                << name() << "; updating mesh flux now"
                << endl;
        }

        // Scale the mesh flux

        const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
        const cyclicACMIFvPatch& nbrACMI = neighbPatch();
        const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();

        const labelListList& newSrcAddr = AMI().srcAddress();
        const labelListList& newTgtAddr = AMI().tgtAddress();

        const fvMesh& mesh = boundaryMesh().mesh();
        surfaceScalarField& meshPhi = const_cast<fvMesh&>(mesh).setPhi().ref();
        surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();

        // Note: phip and phiNonOverlap will be different sizes if new faces
        // have been added
        scalarField& phip = meshPhiBf[cyclicACMIPolyPatch_.index()];
        scalarField& phiNonOverlapp =
            meshPhiBf[nonOverlapPatch.patch().index()];

        const auto& points = mesh.points();

        forAll(phip, facei)
        {
            if (newSrcAddr[facei].empty())
            {
                // AMI patch with no connection to other coupled faces
                phip[facei] = 0.0;
            }
            else
            {
                // Scale the mesh flux according to the area fraction
                const face& fAMI = cyclicACMIPolyPatch_[facei];

                // Note: using raw point locations to calculate the geometric
                // area - faces areas are currently scaled (decoupled from
                // mesh points)
                const scalar geomArea = fAMI.mag(points);
                phip[facei] *= magSf()[facei]/geomArea;
            }
        }

        forAll(phiNonOverlapp, facei)
        {
            const scalar w = 1.0 - cyclicACMIPolyPatch_.srcMask()[facei];
            phiNonOverlapp[facei] *= w;
        }

        const cyclicACMIPolyPatch& nbrPatch = nbrACMI.cyclicACMIPatch();
        scalarField& nbrPhip = meshPhiBf[nbrPatch.index()];
        scalarField& nbrPhiNonOverlapp =
            meshPhiBf[nbrNonOverlapPatch.patch().index()];

        forAll(nbrPhip, facei)
        {
            if (newTgtAddr[facei].empty())
            {
                nbrPhip[facei] = 0.0;
            }
            else
            {
                const face& fAMI = nbrPatch[facei];

                // Note: using raw point locations to calculate the geometric
                // area - faces areas are currently scaled (decoupled from
                // mesh points)
                const scalar geomArea = fAMI.mag(points);
                nbrPhip[facei] *= nbrACMI.magSf()[facei]/geomArea;
            }
        }

        forAll(nbrPhiNonOverlapp, facei)
        {
            const scalar w = 1.0 - cyclicACMIPolyPatch_.tgtMask()[facei];
            nbrPhiNonOverlapp[facei] *= w;
        }
    }
}


// ************************************************************************* //
