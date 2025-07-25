/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "fvMatrix.H"
#include "cyclicACMIFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p)),
    patchNeighbourFieldPtr_(nullptr)
{}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p, dict)),
    patchNeighbourFieldPtr_(nullptr)
{
    if (!isA<cyclicACMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (cacheNeighbourField())
    {
        // Handle neighbour value first, before any evaluate()
        const auto* hasNeighbValue =
            dict.findEntry("neighbourValue", keyType::LITERAL);

        if (hasNeighbValue)
        {
            patchNeighbourFieldPtr_.reset
            (
                new Field<Type>(*hasNeighbValue, p.size())
            );
        }
    }

    // Use 'value' supplied, or set to coupled field
    if (!this->readValueEntry(dict) && this->coupled())
    {
        // Extra check: make sure that the non-overlap patch is before
        // this so it has actually been read - evaluate will crash otherwise

        const GeometricField<Type, fvPatchField, volMesh>& fld =
            static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
            (
                this->primitiveField()
            );
        if (!fld.boundaryField().set(cyclicACMIPatch_.nonOverlapPatchID()))
        {
            FatalIOErrorInFunction(dict)
                << "    patch " << p.name()
                << " of field " << this->internalField().name()
                << " refers to non-overlap patch "
                << cyclicACMIPatch_.cyclicACMIPatch().nonOverlapPatchName()
                << " which is not constructed yet." << nl
                << "    Either supply an initial value or change the ordering"
                << " in the file"
                << exit(FatalIOError);
        }

        // Tricky: avoid call to evaluate without call to initEvaluate.
        // For now just disable the localConsistency to make it use the
        // old logic (ultimately calls the fully self contained
        // patchNeighbourField)

        int& consistency =
            GeometricField<Type, fvPatchField, volMesh>::
            Boundary::localConsistency;

        const int oldConsistency = consistency;
        consistency = 0;

        this->evaluate(Pstream::commsTypes::buffered);

        consistency = oldConsistency;
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p)),
    patchNeighbourFieldPtr_(nullptr)
{
    //if (ptf.patchNeighbourFieldPtr_ && cacheNeighbourField())
    //{
    //    patchNeighbourFieldPtr_.reset
    //    (
    //        new Field<Type>(ptf.patchNeighbourFieldPtr_(), mapper)
    //    );
    //}

    if (!isA<cyclicACMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << cyclicACMIPatch_.name()
            << abort(FatalError);
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_),
    patchNeighbourFieldPtr_(nullptr)
{
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << cyclicACMIPatch_.name()
            << abort(FatalError);
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_),
    patchNeighbourFieldPtr_(nullptr)
{
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << cyclicACMIPatch_.name()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::coupled() const
{
    return cyclicACMIPatch_.coupled();
}


template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::all_ready() const
{
    int done = 0;

    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
    )
    {
        recvRequests_.clear();
        ++done;
    }

    if
    (
        UPstream::finishedRequests
        (
            sendRequests_.start(),
            sendRequests_.size()
        )
    )
    {
        sendRequests_.clear();
        ++done;
    }

    return (done == 2);
}


template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::ready() const
{
    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
    )
    {
        recvRequests_.clear();

        if
        (
            UPstream::finishedRequests
            (
                sendRequests_.start(),
                sendRequests_.size()
            )
        )
        {
            sendRequests_.clear();
        }

        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicACMIFvPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    DebugPout
        << "cyclicACMIFvPatchField::patchNeighbourField(const Field<Type>&) :"
        << " field:" << this->internalField().name()
        << " patch:" << this->patch().name()
        << endl;

    // By pass polyPatch to get nbrId. Instead use cyclicACMIFvPatch virtual
    // neighbPatch()
    const cyclicACMIFvPatch& neighbPatch = cyclicACMIPatch_.neighbPatch();
    const labelUList& nbrFaceCells = neighbPatch.faceCells();

    tmp<Field<Type>> tpnf
    (
        cyclicACMIPatch_.interpolate
        (
            Field<Type>
            (
                iField,
                nbrFaceCells
                //cpp.neighbPatch().faceCells()
            )
        )
    );

    if (doTransform())
    {
        transform(tpnf.ref(), forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::cacheNeighbourField()
{
    /*
    return
    (
        GeometricField<Type, fvPatchField, volMesh>::Boundary::localConsistency
     != 0
    );
    */
    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicACMIFvPatchField<Type>::patchNeighbourField() const
{
    if (this->ownerAMI().distributed() && cacheNeighbourField())
    {
        if (!this->ready())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Initialise if not done in construct-from-dictionary
        if (!patchNeighbourFieldPtr_)
        {
            DebugPout
                << "cyclicACMIFvPatchField::patchNeighbourField() :"
                << " field:" << this->internalField().name()
                << " patch:" << this->patch().name()
                << " calculating&caching patchNeighbourField"
                << endl;

            // Do interpolation and store result
            patchNeighbourFieldPtr_.reset
            (
                patchNeighbourField(this->primitiveField()).ptr()
            );
        }
        else
        {
            DebugPout
                << "cyclicACMIFvPatchField::patchNeighbourField() :"
                << " field:" << this->internalField().name()
                << " patch:" << this->patch().name()
                << " returning cached patchNeighbourField"
                << endl;
        }
        return patchNeighbourFieldPtr_();
    }
    else
    {
        // Do interpolation
        DebugPout
            << "cyclicACMIFvPatchField::evaluate() :"
            << " field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << " calculating up-to-date patchNeighbourField"
            << endl;

        return patchNeighbourField(this->primitiveField());
    }
}


template<class Type>
const Foam::cyclicACMIFvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicACMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicACMIPatch_.neighbPatchID()]
    );
}


template<class Type>
const Foam::fvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::nonOverlapPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    // WIP: Needs to re-direct nonOverlapPatchID to new patchId for assembly?
    return fld.boundaryField()[cyclicACMIPatch_.nonOverlapPatchID()];
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (this->ownerAMI().distributed() && cacheNeighbourField())
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            // Invalidate old field - or flag as fatal?
            patchNeighbourFieldPtr_.reset(nullptr);
            return;
        }

        // Start sending
        DebugPout
            << "cyclicACMIFvPatchField::initEvaluate() :"
            << " field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << " starting send&receive"
            << endl;

        // Bypass polyPatch to get nbrId.
        // - use cyclicACMIFvPatch::neighbPatch() virtual instead
        const cyclicACMIFvPatch& neighbPatch = cyclicACMIPatch_.neighbPatch();
        const labelUList& nbrFaceCells = neighbPatch.faceCells();
        const Field<Type> pnf(this->primitiveField(), nbrFaceCells);

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();

        cyclicACMIPatch_.initInterpolate
        (
            pnf,
            sendRequests_,
            sendBufs_,
            recvRequests_,
            recvBufs_
        );
    }
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const auto& AMI = this->ownerAMI();

    if (AMI.distributed() && cacheNeighbourField())
    {
        // Calculate patchNeighbourField
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        patchNeighbourFieldPtr_.reset(nullptr);

        if (!this->ready())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        patchNeighbourFieldPtr_.reset
        (
            cyclicACMIPatch_.interpolate
            (
                Field<Type>::null(),    // Not used for distributed
                recvRequests_,
                recvBufs_
            ).ptr()
        );

        // Receive requests all handled by last function call
        recvRequests_.clear();


        auto& patchNeighbourField = patchNeighbourFieldPtr_.ref();

        if (doTransform())
        {
            // In-place transform
            transform(patchNeighbourField, forwardT(), patchNeighbourField);
        }
     }

    // Use patchNeighbourField() and patchInternalField() to obtain face value
    coupledFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->ownerAMI().distributed())
    {
        // Start sending
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        DebugPout
            << "cyclicACMIFvPatchField::initInterfaceMatrixUpdate() :"
            << " field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << " starting send&receive"
            << endl;

        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

        solveScalarField pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();

        cyclicACMIPatch_.initInterpolate
        (
            pnf,
            sendRequests_,
            scalarSendBufs_,
            recvRequests_,
            scalarRecvBufs_
        );
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    DebugPout<< "cyclicACMIFvPatchField::updateInterfaceMatrix() :"
        << " field:" << this->internalField().name()
        << " patch:" << this->patch().name()
        << endl;

    // note: only applying coupled contribution

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    solveScalarField pnf;

    if (this->ownerAMI().distributed())
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        DebugPout
            << "cyclicACMIFvPatchField::evaluate() :"
            << " field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << " consuming received coupled neighbourfield"
            << endl;

        pnf =
            cyclicACMIPatch_.interpolate
            (
                solveScalarField::null(),    // Not used for distributed
                recvRequests_,
                scalarRecvBufs_
            );

        // Receive requests all handled by last function call
        recvRequests_.clear();
    }
    else
    {
        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

        pnf = solveScalarField(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        pnf = cyclicACMIPatch_.interpolate(pnf);
    }

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);

    this->updatedMatrix(true);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    const auto& AMI = this->ownerAMI();

    if (AMI.distributed())
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

        Field<Type> pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf);

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();

        cyclicACMIPatch_.initInterpolate
        (
            pnf,
            sendRequests_,
            sendBufs_,
            recvRequests_,
            recvBufs_
        );
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    // note: only applying coupled contribution

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    const auto& AMI = this->ownerAMI();

    Field<Type> pnf;

    if (AMI.distributed())
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        pnf =
            cyclicACMIPatch_.interpolate
            (
                Field<Type>::null(),    // Not used for distributed
                recvRequests_,
                recvBufs_
            );

        // Receive requests all handled by last function call
        recvRequests_.clear();
    }
    else
    {
        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

        pnf = Field<Type>(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf);

        pnf = cyclicACMIPatch_.interpolate(pnf);
    }

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);

    this->updatedMatrix(true);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

    // Nothing to be done by the AMI, but re-direct to non-overlap patch
    // with non-overlap patch weights
    const fvPatchField<Type>& npf = nonOverlapPatchField();

    const_cast<fvPatchField<Type>&>(npf).manipulateMatrix(matrix, 1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label mat,
    const direction cmpt
)
{
    if (this->cyclicACMIPatch().owner())
    {
        label index = this->patch().index();

        const label globalPatchID =
            matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

        const Field<scalar> intCoeffsCmpt
        (
            matrix.internalCoeffs()[globalPatchID].component(cmpt)
        );

        const Field<scalar> boundCoeffsCmpt
        (
            matrix.boundaryCoeffs()[globalPatchID].component(cmpt)
        );

        tmp<Field<scalar>> tintCoeffs(coeffs(matrix, intCoeffsCmpt, mat));
        tmp<Field<scalar>> tbndCoeffs(coeffs(matrix, boundCoeffsCmpt, mat));
        const Field<scalar>& intCoeffs = tintCoeffs.ref();
        const Field<scalar>& bndCoeffs = tbndCoeffs.ref();

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        label subFaceI = 0;

        const labelList& faceMap =
            matrix.lduMeshAssembly().faceBoundMap()[mat][index];

        forAll (faceMap, j)
        {
            label globalFaceI = faceMap[j];

            const scalar boundCorr = -bndCoeffs[subFaceI];
            const scalar intCorr = -intCoeffs[subFaceI];

            matrix.upper()[globalFaceI] += boundCorr;
            matrix.diag()[u[globalFaceI]] -= intCorr;
            matrix.diag()[l[globalFaceI]] -= boundCorr;

            if (matrix.asymmetric())
            {
                matrix.lower()[globalFaceI] += intCorr;
            }
            subFaceI++;
        }

        // Set internalCoeffs and boundaryCoeffs in the assembly matrix
        // on cyclicACMI patches to be used in the individual matrix by
        // matrix.flux()
        if (matrix.psi(mat).mesh().fluxRequired(this->internalField().name()))
        {
            matrix.internalCoeffs().set
            (
                globalPatchID, intCoeffs*pTraits<Type>::one
            );
            matrix.boundaryCoeffs().set
            (
                globalPatchID, bndCoeffs*pTraits<Type>::one
            );

            const label nbrPathID =
                cyclicACMIPatch_.cyclicACMIPatch().neighbPatchID();

            const label nbrGlobalPatchID =
                matrix.lduMeshAssembly().patchLocalToGlobalMap()
                [mat][nbrPathID];

            matrix.internalCoeffs().set
            (
                nbrGlobalPatchID, intCoeffs*pTraits<Type>::one
            );
            matrix.boundaryCoeffs().set
            (
                nbrGlobalPatchID, bndCoeffs*pTraits<Type>::one
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::cyclicACMIFvPatchField<Type>::coeffs
(
    fvMatrix<Type>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    const label index(this->patch().index());

    const label nSubFaces
    (
        matrix.lduMeshAssembly().cellBoundMap()[mat][index].size()
    );

    auto tmapCoeffs = tmp<Field<scalar>>::New(nSubFaces, Zero);
    auto& mapCoeffs = tmapCoeffs.ref();

    const scalarListList& srcWeight =
        cyclicACMIPatch_.cyclicACMIPatch().AMI().srcWeights();

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

    const scalar tol = cyclicACMIPolyPatch::tolerance();
    label subFaceI = 0;
    forAll(mask, faceI)
    {
        const scalarList& w = srcWeight[faceI];
        for(label i=0; i<w.size(); i++)
        {
            if (mask[faceI] > tol)
            {
                const label localFaceId =
                    matrix.lduMeshAssembly().facePatchFaceMap()
                    [mat][index][subFaceI];
                mapCoeffs[subFaceI] = w[i]*coeffs[localFaceId];
            }
            subFaceI++;
        }
    }

    return tmapCoeffs;
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateCoeffs()
{
    // Update non-overlap patch - some will implement updateCoeffs, and
    // others will implement evaluate

    // Pass in (1 - mask) to give non-overlap patch the chance to do
    // manipulation of non-face based data

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();
    const fvPatchField<Type>& npf = nonOverlapPatchField();
    const_cast<fvPatchField<Type>&>(npf).updateWeightedCoeffs(1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    fvPatchField<Type>::writeValueEntry(os);

    if (patchNeighbourFieldPtr_)
    {
        patchNeighbourFieldPtr_->writeEntry("neighbourValue", os);
    }
}


// ************************************************************************* //
