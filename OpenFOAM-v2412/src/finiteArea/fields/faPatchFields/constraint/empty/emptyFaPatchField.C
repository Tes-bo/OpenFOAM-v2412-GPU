/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "emptyFaPatchField.H"
#include "faPatchFieldMapper.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{}


template<class Type>
Foam::emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>&,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper&
)
:
    faPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    if (!isType<emptyFaPatch>(p))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    // Empty means empty, so no patchType override
    // with faPatchFieldBase::readDict(dict);

    if (!isA<emptyFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf.patch(), iF, Field<Type>())  // zero-sized
{}


template<class Type>
Foam::emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>& ptf
)
:
    emptyFaPatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::emptyFaPatchField<Type>::updateCoeffs()
{
    if (this->internalField().mesh().nFaces())
    {
        if
        (
            this->patch().faPatch::size()
          % this->internalField().mesh().nFaces()
        )
        {
            FatalErrorInFunction
                << "This mesh contains patches of type empty but is "
                << "not 1D or 2D by virtue of the fact that the number of "
                << "faces of this empty patch is not divisible by the number "
                << "of cells."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::emptyFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    // Never write "value"
}


// ************************************************************************* //
