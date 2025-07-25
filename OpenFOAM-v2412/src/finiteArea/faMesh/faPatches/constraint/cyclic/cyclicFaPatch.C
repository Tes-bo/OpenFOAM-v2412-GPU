/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "cyclicFaPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicFaPatch, 0);
    addToRunTimeSelectionTable(faPatch, cyclicFaPatch, dictionary);
}

const Foam::scalar Foam::cyclicFaPatch::matchTol_ = 1e-3;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicFaPatch::cyclicFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm,
    const word& patchType
)
:
    coupledFaPatch(name, dict, index, bm, patchType)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicFaPatch::calcTransforms()
{
    const label sizeby2 = this->size()/2;
    if (size() > 0)
    {
        pointField half0Ctrs(sizeby2);
        pointField half1Ctrs(sizeby2);
        for (label edgei=0; edgei < sizeby2; ++edgei)
        {
            half0Ctrs[edgei] = this->edgeCentres()[edgei];
            half1Ctrs[edgei] = this->edgeCentres()[edgei + sizeby2];
        }

        vectorField half0Normals(sizeby2);
        vectorField half1Normals(sizeby2);

        const vectorField eN(edgeNormals()*magEdgeLengths());

        scalar maxMatchError = 0;
        label errorEdge = -1;

        for (label edgei = 0; edgei < sizeby2; ++edgei)
        {
            half0Normals[edgei] = eN[edgei];
            half1Normals[edgei] = eN[edgei + sizeby2];

            scalar magLe = mag(half0Normals[edgei]);
            scalar nbrMagLe = mag(half1Normals[edgei]);
            scalar avLe = 0.5*(magLe + nbrMagLe);

            if (magLe < ROOTVSMALL && nbrMagLe < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                half0Normals[edgei] = point(1, 0, 0);
                half1Normals[edgei] = half0Normals[edgei];
            }
            else if (mag(magLe - nbrMagLe)/avLe > matchTol_)
            {
                // Error in area matching.  Find largest error
                maxMatchError =
                    Foam::max(maxMatchError, mag(magLe - nbrMagLe)/avLe);
                errorEdge = edgei;
            }
            else
            {
                half0Normals[edgei] /= magLe;
                half1Normals[edgei] /= nbrMagLe;
            }
        }

        // Check for error in edge matching
        if (maxMatchError > matchTol_)
        {
            label nbrEdgei = errorEdge + sizeby2;
            scalar magLe = mag(half0Normals[errorEdge]);
            scalar nbrMagLe = mag(half1Normals[errorEdge]);
            scalar avLe = 0.5*(magLe + nbrMagLe);

            FatalErrorInFunction
                << "edge " << errorEdge
                << " area does not match neighbour "
                << nbrEdgei << " by "
                << 100*mag(magLe - nbrMagLe)/avLe
                << "% -- possible edge ordering problem." << endl
                << "patch:" << name()
                << " my area:" << magLe
                << " neighbour area:" << nbrMagLe
                << " matching tolerance:" << matchTol_
                << endl
                << "Mesh edge:" << start() + errorEdge
                << endl
                << "Neighbour edge:" << start() + nbrEdgei
                << endl
                << "Other errors also exist, only the largest is reported. "
                << "Please rerun with cyclic debug flag set"
                << " for more information." << exit(FatalError);
        }

        // Calculate transformation tensors
        calcTransformTensors
        (
            half0Ctrs,
            half1Ctrs,
            half0Normals,
            half1Normals
        );

        // Check transformation tensors
        if (!parallel())
        {
            if (forwardT().size() > 1 || reverseT().size() > 1)
            {
                SeriousErrorInFunction
                    << "Transformation tensor is not constant for the cyclic "
                    << "patch.  Please reconsider your setup and definition of "
                    << "cyclic boundaries." << endl;
            }
        }
    }
}


void Foam::cyclicFaPatch::makeWeights(scalarField& w) const
{
    const scalarField& magL = magEdgeLengths();

    const scalarField deltas(edgeNormals() & coupledFaPatch::delta());
    const label sizeby2 = deltas.size()/2;

    scalar maxMatchError = 0;
    label errorEdge = -1;

    for (label edgei = 0; edgei < sizeby2; ++edgei)
    {
        scalar avL = 0.5*(magL[edgei] + magL[edgei + sizeby2]);

        if
        (
            mag(magL[edgei] - magL[edgei + sizeby2])/avL
          > matchTol_
        )
        {
            // Found error.  Look for largest matching error
            maxMatchError =
                Foam::max
                (
                    maxMatchError,
                    mag(magL[edgei] - magL[edgei + sizeby2])/avL
                );

            errorEdge = edgei;
        }

        scalar di = deltas[edgei];
        scalar dni = deltas[edgei + sizeby2];

        w[edgei] = dni/(di + dni);
        w[edgei + sizeby2] = 1 - w[edgei];
    }

    // Check for error in matching
    if (maxMatchError > matchTol_)
    {
        scalar avL = 0.5*(magL[errorEdge] + magL[errorEdge + sizeby2]);

        FatalErrorInFunction
            << "edge " << errorEdge << " and " << errorEdge + sizeby2
            <<  " areas do not match by "
            << 100*mag(magL[errorEdge] - magL[errorEdge + sizeby2])/avL
            << "% -- possible edge ordering problem." << nl
            << "Cyclic area match tolerance = "
            << matchTol_ << " patch: " << name()
            << abort(FatalError);
    }
}


void Foam::cyclicFaPatch::makeLPN(scalarField& lPN) const
{
    makeDeltaCoeffs(lPN);
    lPN = scalar(1)/lPN;
}


void Foam::cyclicFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    const scalarField deltas(edgeNormals() & coupledFaPatch::delta());
    const label sizeby2 = deltas.size()/2;

    for (label edgei = 0; edgei < sizeby2; ++edgei)
    {
        scalar di = deltas[edgei];
        scalar dni = deltas[edgei + sizeby2];

        dc[edgei] = 1.0/(di + dni);
        dc[edgei + sizeby2] = dc[edgei];
    }
}


void Foam::cyclicFaPatch::initGeometry(PstreamBuffers& pBufs)
{
    faPatch::initGeometry(pBufs);
}


void Foam::cyclicFaPatch::calcGeometry(PstreamBuffers& pBufs)
{
    faPatch::calcGeometry(pBufs);
    calcTransforms();
}


void Foam::cyclicFaPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    faPatch::initMovePoints(pBufs, p);
}


void Foam::cyclicFaPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    faPatch::movePoints(pBufs, p);
    calcTransforms();
}


Foam::tmp<Foam::vectorField> Foam::cyclicFaPatch::delta() const
{
    const vectorField patchD(coupledFaPatch::delta());
    const label sizeby2 = patchD.size()/2;

    auto tpdv = tmp<vectorField>::New(patchD.size());
    auto& pdv = tpdv.ref();

    // Do the transformation if necessary
    if (parallel())
    {
        for (label edgei = 0; edgei < sizeby2; ++edgei)
        {
            const vector& ddi = patchD[edgei];
            const vector& dni = patchD[edgei + sizeby2];

            pdv[edgei] = ddi - dni;
            pdv[edgei + sizeby2] = -pdv[edgei];
        }
    }
    else
    {
        for (label edgei = 0; edgei < sizeby2; ++edgei)
        {
            const vector& ddi = patchD[edgei];
            const vector& dni = patchD[edgei + sizeby2];

            pdv[edgei] = ddi - transform(forwardT()[0], dni);
            pdv[edgei + sizeby2] = -transform(reverseT()[0], pdv[edgei]);
        }
    }

    return tpdv;
}


Foam::tmp<Foam::labelField> Foam::cyclicFaPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicFaPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& edgeFaces
) const
{
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, edgeFaces, tpfld.ref());
    return tpfld;
}


Foam::tmp<Foam::labelField> Foam::cyclicFaPatch::transfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    auto tpnf = tmp<labelField>::New(this->size());
    auto& pnf = tpnf.ref();

    const label sizeby2 = this->size()/2;

    for (label edgei=0; edgei < sizeby2; ++edgei)
    {
        pnf[edgei] = interfaceData[edgei + sizeby2];
        pnf[edgei + sizeby2] = interfaceData[edgei];
    }

    return tpnf;
}


Foam::tmp<Foam::labelField> Foam::cyclicFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return internalFieldTransfer(commsType, iF, this->faceCells());
}


Foam::tmp<Foam::labelField> Foam::cyclicFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF,
    const labelUList& edgeCells
) const
{
    auto tpnf = tmp<labelField>::New(this->size());
    auto& pnf = tpnf.ref();

    const label sizeby2 = this->size()/2;

    for (label edgei=0; edgei < sizeby2; ++edgei)
    {
        pnf[edgei] = iF[edgeCells[edgei + sizeby2]];
        pnf[edgei + sizeby2] = iF[edgeCells[edgei]];
    }

    return tpnf;
}


// ************************************************************************* //
