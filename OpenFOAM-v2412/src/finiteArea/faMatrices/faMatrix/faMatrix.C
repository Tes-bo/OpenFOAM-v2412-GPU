/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "areaFields.H"
#include "edgeFields.H"
#include "extrapolatedCalculatedFaPatchFields.H"
#include "IndirectList.H"
#include "UniformList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void Foam::faMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
        else if (couples)
        {
            tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const labelUList& addr = lduAddr().patchAddr(patchI);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::faMatrix<Type>::faMatrix
(
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), Zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size())
{
    DebugInFunction
        << "constructing faMatrix<Type> for field " << psi_.name()
        << endl;

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchi)
    {
        internalCoeffs_.set
        (
            patchi,
            new Field<Type>(psi.mesh().boundary()[patchi].size(), Zero)
        );

        boundaryCoeffs_.set
        (
            patchi,
            new Field<Type>(psi.mesh().boundary()[patchi].size(), Zero)
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    auto& psiRef =
        const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    const label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::faMatrix<Type>::faMatrix(const faMatrix<Type>& fam)
:
    lduMatrix(fam),
    psi_(fam.psi_),
    dimensions_(fam.dimensions_),
    source_(fam.source_),
    internalCoeffs_(fam.internalCoeffs_),
    boundaryCoeffs_(fam.boundaryCoeffs_)
{
    DebugInFunction
        << "Copying faMatrix<Type> for field " << psi_.name() << endl;

    if (fam.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *(fam.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
Foam::faMatrix<Type>::faMatrix(const tmp<faMatrix<Type>>& tmat)
:
    lduMatrix(tmat.constCast(), tmat.movable()),
    psi_(tmat().psi_),
    dimensions_(tmat().dimensions_),
    source_(tmat.constCast().source_, tmat.movable()),
    internalCoeffs_(tmat.constCast().internalCoeffs_, tmat.movable()),
    boundaryCoeffs_(tmat.constCast().boundaryCoeffs_, tmat.movable())
{
    DebugInFunction
        << "Copy/Move faMatrix<Type> for field " << psi_.name() << endl;

    if (tmat().faceFluxCorrectionPtr_)
    {
        if (tmat.movable())
        {
            faceFluxCorrectionPtr_ =
                std::move(tmat.constCast().faceFluxCorrectionPtr_);
        }
        else if (tmat().faceFluxCorrectionPtr_)
        {
            faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
            (
                *(tmat().faceFluxCorrectionPtr_)
            );
        }
    }

    tmat.clear();
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::faMatrix<Type>::~faMatrix()
{
    DebugInFunction
        << "Destroying faMatrix<Type> for field " << psi_.name() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class ListType>
void Foam::faMatrix<Type>::setValuesFromList
(
    const labelUList& faceLabels,
    const ListType<Type>& values
)
{
#if 0  /* Specific to foam-extend */
    // Record face labels of eliminated equations
    for (const label i : faceLabels)
    {
        this->eliminatedEqns().insert(i);
    }
#endif

    const faMesh& mesh = psi_.mesh();

    const labelListList& edges = mesh.patch().faceEdges();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField& Diag = diag();
    Field<Type>& psi =
        const_cast
        <
            GeometricField<Type, faPatchField, areaMesh>&
        >(psi_).primitiveFieldRef();

    // Following actions:
    // - adjust local field psi
    // - set local matrix to be diagonal (so adjust source)
    //      - cut connections to neighbours
    // - make (on non-adjusted cells) contribution explicit

    if (symmetric() || asymmetric())
    {
        forAll(faceLabels, i)
        {
            const label facei = faceLabels[i];
            const Type& value = values[i];

            for (const label edgei : edges[facei])
            {
                if (mesh.isInternalEdge(edgei))
                {
                    if (symmetric())
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= upper()[edgei]*value;
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*value;
                        }

                        upper()[edgei] = 0.0;
                    }
                    else
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= lower()[edgei]*value;
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*value;
                        }

                        upper()[edgei] = 0.0;
                        lower()[edgei] = 0.0;
                    }
                }
                else
                {
                    const label patchi = mesh.boundary().whichPatch(edgei);

                    if (internalCoeffs_[patchi].size())
                    {
                        const label patchEdgei =
                            mesh.boundary()[patchi].whichEdge(edgei);

                        internalCoeffs_[patchi][patchEdgei] = Zero;
                        boundaryCoeffs_[patchi][patchEdgei] = Zero;
                    }
                }
            }
        }
    }

    // Note: above loop might have affected source terms on adjusted cells
    // so make sure to adjust them afterwards
    forAll(faceLabels, i)
    {
        const label facei = faceLabels[i];
        const Type& value = values[i];

        psi[facei] = value;
        source_[facei] = value*Diag[facei];
    }
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const Type& value
)
{
    this->setValuesFromList(faceLabels, UniformList<Type>(value));
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const UList<Type>& values
)
{
    this->setValuesFromList(faceLabels, values);
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(faceLabels, values);
}


template<class Type>
void Foam::faMatrix<Type>::setReference
(
    const label facei,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && facei >= 0)
    {
        if (Pstream::master())
        {
            source()[facei] += diag()[facei]*value;
            diag()[facei] += diag()[facei];
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::setReferences
(
    const labelUList& faceLabels,
    const Type& value,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        forAll(faceLabels, facei)
        {
            const label faceId = faceLabels[facei];
            if (faceId >= 0)
            {
                source()[faceId] += diag()[faceId]*value;
                diag()[faceId] += diag()[faceId];
            }
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::setReferences
(
    const labelUList& faceLabels,
    const UList<Type>& values,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        forAll(faceLabels, facei)
        {
            const label faceId = faceLabels[facei];
            if (faceId >= 0)
            {
                source()[faceId] += diag()[faceId]*values[facei];
                diag()[faceId] += diag()[faceId];
            }
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), Zero);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries subtract the diagonal
                // contribution off-diagonal sum which avoids having to remove
                // it from the diagonal later.
                // Also add the source contribution from the relaxation
                forAll(pa, face)
                {
                    Type iCoeff0 = iCoeffs[face];
                    iCoeffs[face] = cmptMag(iCoeffs[face]);
                    sumOff[pa[face]] -= cmptMin(iCoeffs[face]);
                    iCoeffs[face] /= alpha;
                    S[pa[face]] +=
                        cmptMultiply(iCoeffs[face] - iCoeff0, psi_[pa[face]]);
                }
            }
        }
    }

    // Ensure the matrix is diagonally dominant...
    max(D, D, sumOff);

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.primitiveField();
}


template<class Type>
void Foam::faMatrix<Type>::relax()
{
    scalar relaxCoeff = 0;

    if (psi_.mesh().relaxEquation(psi_.name(), relaxCoeff))
    {
        relax(relaxCoeff);
    }
    else
    {
        DebugInFunction
            << "No relaxation specified for field " << psi_.name() << nl;
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::faMatrix<Type>::D() const
{
    auto tdiag = tmp<scalarField>::New(diag());
    addCmptAvBoundaryDiag(tdiag.ref());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::areaScalarField> Foam::faMatrix<Type>::A() const
{
    auto tAphi = areaScalarField::New
    (
        "A(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/psi_.dimensions()/dimArea,
        faPatchFieldBase::extrapolatedCalculatedType()
    );

    tAphi.ref().primitiveFieldRef() = D()/psi_.mesh().S();
    tAphi.ref().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::faMatrix<Type>::H() const
{
    auto tHphi = GeometricField<Type, faPatchField, areaMesh>::New
    (
        "H(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/dimArea,
        faPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Hphi = tHphi.ref();

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; ++cmpt)
    {
        scalarField psiCmpt(psi_.primitiveField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), Zero);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.primitiveFieldRef().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.primitiveFieldRef() += lduMatrix::H(psi_.primitiveField()) + source_;
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.primitiveFieldRef() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::faMatrix<Type>::flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary of faSchemes"
            << abort(FatalError);
    }

    auto tfieldFlux = GeometricField<Type, faePatchField, edgeMesh>::New
    (
        "flux(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions()
    );
    auto& fieldFlux = tfieldFlux.ref();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
    {
        fieldFlux.primitiveFieldRef().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.primitiveField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryFieldRef()[patchI] =
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


template<class Type>
const Foam::dictionary& Foam::faMatrix<Type>::solverDict
(
    const word& name
) const
{
    return psi_.mesh().solverDict(name);
}


template<class Type>
const Foam::dictionary& Foam::faMatrix<Type>::solverDict() const
{
    return psi_.mesh().solverDict
    (
        psi_.name()
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::faMatrix<Type>::operator=(const faMatrix<Type>& famv)
{
    if (this == &famv)
    {
        return;  // Self-assignment is a no-op
    }

    if (&psi_ != &(famv.psi_))
    {
        FatalErrorInFunction
            << "different fields"
            << abort(FatalError);
    }

    lduMatrix::operator=(famv);
    source_ = famv.source_;
    internalCoeffs_ = famv.internalCoeffs_;
    boundaryCoeffs_ = famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator=(const tmp<faMatrix<Type>>& tfamv)
{
    operator=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator+=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ += famv.dimensions_;
    lduMatrix::operator+=(famv);
    source_ += famv.source_;
    internalCoeffs_ += famv.internalCoeffs_;
    boundaryCoeffs_ += famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator+=(const tmp<faMatrix<Type>>& tfamv)
{
    operator+=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ -= famv.dimensions_;
    lduMatrix::operator-=(famv);
    source_ -= famv.source_;
    internalCoeffs_ -= famv.internalCoeffs_;
    boundaryCoeffs_ -= famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            -*famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator-=(const tmp<faMatrix<Type>>& tfamv)
{
    operator-=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().S()*su.field();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().S()*su.field();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= psi().mesh().S()*su;
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += psi().mesh().S()*su;
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const areaScalarField::Internal& dsf
)
{
    dimensions_ *= dsf.dimensions();
    lduMatrix::operator*=(dsf.field());
    source_ *= dsf.field();

    forAll(boundaryCoeffs_, patchi)
    {
        const scalarField pisf
        (
            dsf.mesh().boundary()[patchi].patchInternalField(dsf.field())
        );
        internalCoeffs_[patchi] *= pisf;
        boundaryCoeffs_[patchi] *= pisf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const tmp<areaScalarField::Internal>& tfld
)
{
    operator*=(tfld());
    tfld.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const tmp<areaScalarField>& tfld
)
{
    operator*=(tfld());
    tfld.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& mat1,
    const faMatrix<Type>& mat2,
    const char* op
)
{
    if (&mat1.psi() != &mat2.psi())
    {
        FatalErrorInFunction
            << "Incompatible fields for operation\n    "
            << "[" << mat1.psi().name() << "] "
            << op
            << " [" << mat2.psi().name() << "]"
            << abort(FatalError);
    }

    if
    (
        dimensionSet::checking()
     && mat1.dimensions() != mat2.dimensions()
    )
    {
        FatalErrorInFunction
            << "Incompatible dimensions for operation\n    "
            << "[" << mat1.psi().name() << mat1.dimensions()/dimArea << " ] "
            << op
            << " [" << mat2.psi().name() << mat2.dimensions()/dimArea << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& mat,
    const DimensionedField<Type, areaMesh>& fld,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && mat.dimensions()/dimArea != fld.dimensions()
    )
    {
        FatalErrorInFunction
            <<  "Incompatible dimensions for operation\n    "
            << "[" << mat.psi().name() << mat.dimensions()/dimArea << " ] "
            << op
            << " [" << fld.name() << fld.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& mat,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && mat.dimensions()/dimArea != dt.dimensions()
    )
    {
        FatalErrorInFunction
            << "Incompatible dimensions for operation\n    "
            << "[" << mat.psi().name() << mat.dimensions()/dimArea << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    faMatrix<Type>& mat,
    const dictionary& solverControls
)
{
    return mat.solve(solverControls);
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<faMatrix<Type>>& tmat,
    const dictionary& solverControls
)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve(solverControls));

    tmat.clear();

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    faMatrix<Type>& mat,
    const word& name
)
{
    return mat.solve(name);
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<faMatrix<Type>>& tmat,
    const word& name
)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve(name));

    tmat.clear();

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(faMatrix<Type>& mat)
{
    return mat.solve();
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(const tmp<faMatrix<Type>>& tmat)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve());

    tmat.clear();

    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    auto tC = tmp<faMatrix<Type>>::New(A);
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A
)
{
    auto tC = tmp<faMatrix<Type>>::New(A);
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    auto tC = tmp<faMatrix<Type>>::New(A);
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(A, su, "+");
    auto tC(A.clone());
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    auto tC(A.clone());
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    auto tC(A.clone());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    auto tC(A.clone());
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    auto tC(A.clone());
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    auto tC(A.clone());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(A, su, "-");
    auto tC(A.clone());
    tC.ref().source() += su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    auto tC(A.clone());
    tC.ref().source() += tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    auto tC(A.clone());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    auto tC(A.clone());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    auto tC(A.clone());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    auto tC(A.clone());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    auto tC(A.clone());
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    auto tC(A.clone());
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    auto tC(A.clone());
    tC.ref().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    auto tC(A.clone());
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(A, su, "==");
    auto tC(A.clone());
    tC.ref().source() += su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    auto tC(A.clone());
    tC.ref().source() += tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    auto tC(A.clone());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const DimensionedField<Type, areaMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.field();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().field();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    auto tC(A.clone());
    tC.ref().source() += A.psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const Foam::zero
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const Foam::zero
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const areaScalarField::Internal& dsf,
    const faMatrix<Type>& A
)
{
    auto tC(A.clone());
    tC.ref() *= dsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField::Internal>& tdsf,
    const faMatrix<Type>& A
)
{
    auto tC(A.clone());
    tC.ref() *= tdsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField>& tvsf,
    const faMatrix<Type>& A
)
{
    auto tC(A.clone());
    tC.ref() *= tvsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const areaScalarField::Internal& dsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= dsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField::Internal>& tdsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tdsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField>& tvsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const faMatrix<Type>& A
)
{
    auto tC(A.clone());
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const faMatrix<Type>& M,
    const DimensionedField<Type, areaMesh>& psi
)
{
    auto tMphi = GeometricField<Type, faPatchField, areaMesh>::New
    (
        "M&" + psi.name(),
        psi.mesh(),
        M.dimensions()/dimArea,
        faPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Mphi = tMphi.ref();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalarField psiCmpt(psi.field().component(cmpt));
            scalarField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.primitiveFieldRef().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.primitiveFieldRef() = Zero;
    }

    Mphi.primitiveFieldRef() += M.lduMatrix::H(psi.field()) + M.source();
    M.addBoundarySource(Mphi.primitiveFieldRef());

    Mphi.primitiveFieldRef() /= -psi.mesh().S();
    Mphi.correctBoundaryConditions();

    return tMphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const faMatrix<Type>& M,
    const tmp<DimensionedField<Type, areaMesh>>& tpsi
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const faMatrix<Type>& M,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tpsi
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const tmp<faMatrix<Type>>& tM,
    const DimensionedField<Type, areaMesh>& psi
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const tmp<faMatrix<Type>>& tM,
    const tmp<DimensionedField<Type, areaMesh>>& tpsi
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::operator&
(
    const tmp<faMatrix<Type>>& tM,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tpsi
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const faMatrix<Type>& fam)
{
    os  << static_cast<const lduMatrix&>(fam) << nl
        << fam.dimensions_ << nl
        << fam.source_ << nl
        << fam.internalCoeffs_ << nl
        << fam.boundaryCoeffs_ << endl;

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "faMatrixSolve.C"

// ************************************************************************* //
