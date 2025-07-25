/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::Amul
(
    Field<Type>& Apsi,
    const tmp<Field<Type>>& tpsi
) const
{
    Type* __restrict__ ApsiPtr = Apsi.begin();

    const Field<Type>& psi = tpsi();
    const Type* const __restrict__ psiPtr = psi.begin();

    const DType* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ upperPtr = upper().begin();
    const LUType* const __restrict__ lowerPtr = lower().begin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfacesUpper_,
        psi,
        Apsi
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = dot(diagPtr[cell], psiPtr[cell]);
    }


    const label nFaces = upper().size();
    for (label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += dot(lowerPtr[face], psiPtr[lPtr[face]]);
        ApsiPtr[lPtr[face]] += dot(upperPtr[face], psiPtr[uPtr[face]]);
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfacesUpper_,
        psi,
        Apsi,
        startRequest
    );

    tpsi.clear();
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::Tmul
(
    Field<Type>& Tpsi,
    const tmp<Field<Type>>& tpsi
) const
{
    Type* __restrict__ TpsiPtr = Tpsi.begin();

    const Field<Type>& psi = tpsi();
    const Type* const __restrict__ psiPtr = psi.begin();

    const DType* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ lowerPtr = lower().begin();
    const LUType* const __restrict__ upperPtr = upper().begin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfacesLower_,
        psi,
        Tpsi
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = dot(diagPtr[cell], psiPtr[cell]);
    }

    const label nFaces = upper().size();
    for (label face=0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += dot(upperPtr[face], psiPtr[lPtr[face]]);
        TpsiPtr[lPtr[face]] += dot(lowerPtr[face], psiPtr[uPtr[face]]);
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfacesLower_,
        psi,
        Tpsi,
        startRequest
    );

    tpsi.clear();
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumA
(
    Field<Type>& sumA
) const
{
    Type* __restrict__ sumAPtr = sumA.begin();

    const DType* __restrict__ diagPtr = diag().begin();

    const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* __restrict__ lowerPtr = lower().begin();
    const LUType* __restrict__ upperPtr = upper().begin();

    const label nCells = diag().size();
    const label nFaces = upper().size();

    for (label cell=0; cell<nCells; cell++)
    {
        sumAPtr[cell] = dot(diagPtr[cell], pTraits<Type>::one);
    }

    for (label face=0; face<nFaces; face++)
    {
        sumAPtr[uPtr[face]] += dot(lowerPtr[face], pTraits<Type>::one);
        sumAPtr[lPtr[face]] += dot(upperPtr[face], pTraits<Type>::one);
    }

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces_, patchi)
    {
        if (interfaces_.set(patchi))
        {
            const labelUList& pa = lduAddr().patchAddr(patchi);
            const Field<LUType>& pCoeffs = interfacesUpper_[patchi];

            forAll(pa, face)
            {
                sumAPtr[pa[face]] -= dot(pCoeffs[face], pTraits<Type>::one);
            }
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::residual
(
    Field<Type>& rA,
    const Field<Type>& psi
) const
{
    Type* __restrict__ rAPtr = rA.begin();

    const Type* const __restrict__ psiPtr = psi.begin();
    const DType* const __restrict__ diagPtr = diag().begin();
    const Type* const __restrict__ sourcePtr = source().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ upperPtr = upper().begin();
    const LUType* const __restrict__ lowerPtr = lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update to add the contibution to the r.h.s.

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        false,          // negate interface contributions
        interfacesUpper_,
        psi,
        rA
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - dot(diagPtr[cell], psiPtr[cell]);
    }


    const label nFaces = upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= dot(lowerPtr[face], psiPtr[lPtr[face]]);
        rAPtr[lPtr[face]] -= dot(upperPtr[face], psiPtr[uPtr[face]]);
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        false,
        interfacesUpper_,
        psi,
        rA,
        startRequest
    );
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::Field<Type>> Foam::LduMatrix<Type, DType, LUType>::residual
(
    const Field<Type>& psi
) const
{
    auto trA = tmp<Field<Type>>::New(psi.size());
    residual(trA.ref(), psi);
    return trA;
}


// ************************************************************************* //
