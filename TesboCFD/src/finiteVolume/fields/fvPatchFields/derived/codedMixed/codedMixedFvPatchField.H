/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::codedMixedFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Constructs on-the-fly a new boundary condition
    (derived from mixedFvPatchField) which is then used to evaluate.

    The code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions;
       code        | c++; patch value assignment
       codeContext | additional dictionary context for the code
    \endplaintable

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedMixed;

        refValue        uniform (0 0 0);
        refGradient     uniform (0 0 0);
        valueFraction   uniform 1;

        name    rampedMixed;   // name of generated BC

        code
        #{
            this->refValue() =
                vector(1, 0, 0)
               *min(10, 0.1*this->db().time().value());
            this->refGrad() = Zero;
            this->valueFraction() = 1.0;
        #};

        //codeInclude
        //#{
        //    #include "fvCFD.H"
        //#};

        //codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry

    \verbatim
    <patchName>
    {
        code
        #{
            this->refValue() = min(10, 0.1*this->db().time().value());
            this->refGrad() = Zero;
            this->valueFraction() = 1.0;
        #};
    }
    \endverbatim

Note
    The code context dictionary can be supplied separately as the
    \c codeContext entry.

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    codedMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_codedMixedFvPatchField_H
#define Foam_codedMixedFvPatchField_H

#include "mixedFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class codedMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedMixedFvPatchField
:
    public mixedFvPatchField<Type>,
    public codedBase
{
    //- The parent boundary condition type
    typedef mixedFvPatchField<Type> parent_bctype;


    // Private Data

        //- Dictionary contents for the boundary condition
        dictionary dict_;

        const word name_;

        mutable autoPtr<mixedFvPatchField<Type>> redirectPatchFieldPtr_;


protected:

    // Protected Member Functions

        //- Mutable access to the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Description (type + name) for the output
        virtual string description() const;

        //- Clear redirected object(s)
        virtual void clearRedirect() const;

        //- Additional 'codeContext' dictionary to pass through
        virtual const dictionary& codeContext() const;

        //- The code dictionary. Inline "code" or from system/codeDict
        virtual const dictionary& codeDict() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "mixedFvPatchFieldTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "mixedFvPatchFieldTemplate.H";


    //- Runtime type information
    TypeName("codedMixed");


    // Constructors

        //- Construct from patch and internal field
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedMixedFvPatchField
        //  onto a new patch
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        codedMixedFvPatchField
        (
            const codedMixedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        //- Get reference to the underlying patchField
        const mixedFvPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        //  This is only needed to set the updated() flag of the name
        //  to false.
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::buffered
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "codedMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
