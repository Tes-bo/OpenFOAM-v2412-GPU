/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "calculatedFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    fvPatchField<Type>(p, iF, dict, requireValue)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::fvPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    auto* patchTypeCtor = patchConstructorTable(p.type());

    if (patchTypeCtor)
    {
        return patchTypeCtor
        (
            p,
            DimensionedField<Type, volMesh>::null()
        );
    }
    else
    {
        return tmp<fvPatchField<Type>>
        (
            new calculatedFvPatchField<Type>
            (
                p,
                DimensionedField<Type, volMesh>::null()
            )
        );
    }
}


template<class Type>
template<class AnyType>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::fvPatchField<Type>::NewCalculatedType
(
    const fvPatchField<AnyType>& pf
)
{
    return NewCalculatedType(pf.patch());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
void Foam::calculatedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
