/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

Class
    Foam::codedFixedValuePointPatchField

Description
    Constructs on-the-fly a new boundary condition (derived from
    fixedValuePointPatchField) which is then used to evaluate.

    The \c value entry (optional) is used for the initial values.
    Otherwise the code is executed.

    The code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions
       code        | c++; patch value assignment
       codeContext | additional dictionary context for the code
    \endplaintable

    Example:
    \verbatim
    movingWall
    {
        type    codedFixedValue;
        value   uniform 0;
        name    rampedFixedValue;   // name of generated bc

        code
        #{
            operator==
            (
                vector(0,0,1) * min(10, 0.1*this->db().time().value())
            );
        #};

        codeContext
        {
             ...
        }

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

    A special form is if the \c code section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry

    \verbatim
    rampedFixedValue
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

Note
    The code context dictionary can be supplied separately as the
    \c codeContext entry.

See also
    codedFixedValueFvPatchField

SourceFiles
    codedFixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_codedFixedValuePointPatchField_H
#define Foam_codedFixedValuePointPatchField_H

#include "fixedValuePointPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class codedFixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedFixedValuePointPatchField
:
    public fixedValuePointPatchField<Type>,
    protected codedBase
{
    //- The parent boundary condition type
    typedef fixedValuePointPatchField<Type> parent_bctype;


    // Private Data

        //- Dictionary contents for the boundary condition
        dictionary dict_;

        const word name_;

        mutable autoPtr<pointPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

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
            = "fixedValuePointPatchFieldTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "fixedValuePointPatchFieldTemplate.H";


    //- Runtime type information
    TypeName("codedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        codedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedFixedValuePointPatchField
        //  onto a new patch
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        codedFixedValuePointPatchField
        (
            const codedFixedValuePointPatchField<Type>&,
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


    // Member functions

        //- Get reference to the underlying patch
        const pointPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets updated() to false
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
    #include "codedFixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
