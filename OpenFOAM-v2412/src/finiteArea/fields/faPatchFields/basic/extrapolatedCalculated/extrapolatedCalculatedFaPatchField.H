/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::extrapolatedCalculatedFaPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition applies a zero-gradient condition from the patch
    internal field onto the patch faces when \c evaluated but may also be
    assigned.  \c snGrad returns the patch gradient evaluated from the current
    internal and patch field values rather than returning zero.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            extrapolatedCalculated;
    }
    \endverbatim

SourceFiles
    extrapolatedCalculatedFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_extrapolatedCalculatedFaPatchField_H
#define Foam_extrapolatedCalculatedFaPatchField_H

#include "calculatedFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class extrapolatedCalculatedFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extrapolatedCalculatedFaPatchField
:
    public calculatedFaPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("extrapolatedCalculated");
    // fieldTypes::extrapolatedCalculatedTypeName_()


    // Constructors

        //- Construct from patch and internal field
        extrapolatedCalculatedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedCalculatedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField onto a new patch
        extrapolatedCalculatedFaPatchField
        (
            const extrapolatedCalculatedFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedCalculatedFaPatchField
        (
            const extrapolatedCalculatedFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        extrapolatedCalculatedFaPatchField
        (
            const extrapolatedCalculatedFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return faPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::buffered
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extrapolatedCalculatedFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
