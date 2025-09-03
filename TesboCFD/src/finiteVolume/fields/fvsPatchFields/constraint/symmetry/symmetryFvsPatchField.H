/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::symmetryFvsPatchField

Description
    This boundary condition enforces a symmetry constraint

    The "value" entry is MUST_READ.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type    symmetry;
        value   ...;
    }
    \endverbatim

SourceFiles
    symmetryFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_symmetryFvsPatchField_H
#define Foam_symmetryFvsPatchField_H

#include "fvsPatchField.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class symmetryFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryFvsPatchField
:
    public fvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(symmetryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        symmetryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        //- The "value" entry is MUST_READ.
        symmetryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given symmetryFvsPatchField onto a new patch
        symmetryFvsPatchField
        (
            const symmetryFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetryFvsPatchField
        (
            const symmetryFvsPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        symmetryFvsPatchField
        (
            const symmetryFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return fvsPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return fvsPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Write includes "value" entry
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmetryFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
