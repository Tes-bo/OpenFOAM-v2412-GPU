/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::basicSymmetryPointPatchField

Description
    A Symmetry boundary condition for pointField

SourceFiles
    basicsymmetryPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_basicSymmetryPointPatchField_H
#define Foam_basicSymmetryPointPatchField_H

#include "pointPatchField.H"
#include "symmetryPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class basicSymmetryPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class basicSymmetryPointPatchField
:
    public pointPatchField<Type>
{
public:

    // Constructors

        //- Construct from patch and internal field
        basicSymmetryPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        basicSymmetryPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        basicSymmetryPointPatchField
        (
            const basicSymmetryPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Default copy construct
        basicSymmetryPointPatchField
        (
            const basicSymmetryPointPatchField<Type>&
        ) = default;

        //- Construct as copy setting internal field reference
        basicSymmetryPointPatchField
        (
            const basicSymmetryPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return pointPatchField<Type>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return pointPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
        );


    // Member Operators

        //- Inherit assignment operators
        using pointPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicSymmetryPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
