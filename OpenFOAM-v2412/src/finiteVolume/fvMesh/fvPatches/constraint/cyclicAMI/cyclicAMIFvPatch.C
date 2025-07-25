/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "cyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "transform.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicAMIFvPatch, polyPatch);
    addNamedToRunTimeSelectionTable
    (
        fvPatch,
        cyclicAMIFvPatch,
        polyPatch,
        cyclicPeriodicAMI
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::cyclicAMIFvPatch::newInternalProcFaces
// (
//     label& newFaces,
//     label& newProcFaces
// ) const
// {
//     const labelListList& addSourceFaces = AMI().srcAddress();
//
//     // Add new faces as many weights for AMI
//     forAll (addSourceFaces, faceI)
//     {
//         const labelList& nbrFaceIs = addSourceFaces[faceI];
//
//         forAll (nbrFaceIs, j)
//         {
//             label nbrFaceI = nbrFaceIs[j];
//
//             if (nbrFaceI < neighbPatch().size())
//             {
//                 // local faces
//                 newFaces++;
//             }
//             else
//             {
//                 // Proc faces
//                 newProcFaces++;
//             }
//         }
//     }
// }


bool Foam::cyclicAMIFvPatch::coupled() const
{
    return
        Pstream::parRun()
     || !this->boundaryMesh().mesh().time().processorCase();
}


void Foam::cyclicAMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

        const scalarField deltas(nf() & coupledFvPatch::delta());

        tmp<scalarField> tnbrDeltas;
        if (applyLowWeightCorrection())
        {
            tnbrDeltas =
                interpolate
                (
                    nbrPatch.nf() & nbrPatch.coupledFvPatch::delta(),
                    scalarField(this->size(), 1.0)
                );
        }
        else
        {
            tnbrDeltas =
                interpolate(nbrPatch.nf() & nbrPatch.coupledFvPatch::delta());
        }

        const scalarField& nbrDeltas = tnbrDeltas();

        forAll(deltas, facei)
        {
            // Note use of mag
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            w[facei] = dni/(di + dni);
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


void Foam::cyclicAMIFvPatch::makeDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
}


void Foam::cyclicAMIFvPatch::makeNonOrthoDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
    //coeffs = Zero;
}


void Foam::cyclicAMIFvPatch::makeNonOrthoCorrVectors(vectorField& vecs) const
{
    // Apply correction to default vectors
    //vecs = Zero;
}


Foam::tmp<Foam::vectorField> Foam::cyclicAMIFvPatch::delta() const
{
    const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        tmp<vectorField> tnbrPatchD;
        if (applyLowWeightCorrection())
        {
            tnbrPatchD =
                interpolate
                (
                    nbrPatch.coupledFvPatch::delta(),
                    vectorField(this->size(), Zero)
                );
        }
        else
        {
            tnbrPatchD = interpolate(nbrPatch.coupledFvPatch::delta());
        }

        const vectorField& nbrPatchD = tnbrPatchD();

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


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, faceCells, tpfld.ref());
    return tpfld;
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


void Foam::cyclicAMIFvPatch::movePoints()
{
    if (!owner() || !cyclicAMIPolyPatch_.createAMIFaces())
    {
        // Only manipulating patch face areas and mesh motion flux if the AMI
        // creates additional faces
        return;
    }

    // Update face data based on values set by the AMI manipulations
    const_cast<vectorField&>(Sf()) = cyclicAMIPolyPatch_.faceAreas();
    const_cast<vectorField&>(Cf()) = cyclicAMIPolyPatch_.faceCentres();
    const_cast<scalarField&>(magSf()) = mag(Sf());

    const cyclicAMIFvPatch& nbr = neighbPatch();
    const_cast<vectorField&>(nbr.Sf()) = nbr.cyclicAMIPatch().faceAreas();
    const_cast<vectorField&>(nbr.Cf()) = nbr.cyclicAMIPatch().faceCentres();
    const_cast<scalarField&>(nbr.magSf()) = mag(nbr.Sf());


    // Set consitent mesh motion flux
    // TODO: currently maps src mesh flux to tgt - update to
    // src = src + mapped(tgt) and tgt = tgt + mapped(src)?

    const fvMesh& mesh = boundaryMesh().mesh();
    surfaceScalarField& meshPhi = const_cast<fvMesh&>(mesh).setPhi().ref();
    surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();

    if (cyclicAMIPolyPatch_.owner())
    {
        scalarField& phip = meshPhiBf[patch().index()];
        forAll(phip, facei)
        {
            const face& f = cyclicAMIPolyPatch_.localFaces()[facei];

            // Note: using raw point locations to calculate the geometric
            // area - faces areas are currently scaled by the AMI weights
            // (decoupled from mesh points)
            const scalar geomArea = f.mag(cyclicAMIPolyPatch_.localPoints());

            const scalar scaledArea = magSf()[facei];
            phip[facei] *= scaledArea/geomArea;
        }

        scalarField srcMeshPhi(phip);
        if (AMI().distributed())
        {
            AMI().srcMap().distribute(srcMeshPhi);
        }

        const labelListList& tgtToSrcAddr = AMI().tgtAddress();
        scalarField& nbrPhip = meshPhiBf[nbr.index()];

        forAll(tgtToSrcAddr, tgtFacei)
        {
            // Note: now have 1-to-1 mapping so tgtToSrcAddr[tgtFacei] is size 1
            const label srcFacei = tgtToSrcAddr[tgtFacei][0];
            nbrPhip[tgtFacei] = -srcMeshPhi[srcFacei];
        }

        DebugInfo
            << "patch:" << patch().name()
            << " sum(area):" << gSum(magSf())
            << " min(mag(faceAreas):" << gMin(magSf())
            << " sum(meshPhi):" << gSum(phip) << nl
            << " sum(nbrMeshPhi):" << gSum(nbrPhip) << nl
            << endl;
    }
}

// ************************************************************************* //
