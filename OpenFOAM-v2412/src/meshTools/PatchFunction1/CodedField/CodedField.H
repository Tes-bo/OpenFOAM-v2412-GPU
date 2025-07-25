/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Foam::PatchFunction1Types::CodedField

Description
    PatchFunction1 with the code supplied by an on-the-fly compiled C++
    expression.

    The code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions
       code        | c++; return the patch values at (scalar x)
    \endplaintable

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            uniformFixedValue;
        uniformValue
        {
            type  coded;
            name  myExpression;  // Name of generated PatchFunction1

            code
            #{
                const polyPatch& pp = this->patch();
                Pout<< "** Patch size:" << pp.size() << endl;
                return tmp<vectorField>::New(pp.size(), vector(1, 0, 0))
            #};

            //codeInclude
            //#{
            //    #include "volFields.H"
            //#};

            //codeOptions
            //#{
            //    -I$(LIB_SRC)/finiteVolume/lnInclude
            //#};
        }
    }
    \endverbatim

Note
    The code context dictionary is simply the dictionary used to specify
    the PatchFunction1 coefficients.

See also
    Foam::dynamicCode
    Foam::codedFixedValue
    Foam::functionEntries::codeStream

SourceFiles
    CodedField.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_CodedField_H
#define PatchFunction1Types_CodedField_H

#include "PatchFunction1.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                         Class CodedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedField
:
    public PatchFunction1<Type>,
    protected codedBase
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        const word redirectName_;

        mutable autoPtr<PatchFunction1<Type>> redirectFunctionPtr_;


    // Private Member Functions

        //- Get reference to the underlying Function1
        const PatchFunction1<Type>& redirectFunction() const;


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

        // Get the code (sub)dictionary
        virtual const dictionary& codeDict(const dictionary& dict) const;

        // Get the code dictionary
        virtual const dictionary& codeDict() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


        //- No copy assignment
        void operator=(const CodedField<Type>&) = delete;

public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedPatchFunction1Template.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedPatchFunction1Template.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from entry name and dictionary
        CodedField
        (
            const polyPatch& pp,
            const word& redirectType,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct, setting patch
        explicit CodedField
        (
            const CodedField<Type>& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        explicit CodedField(const CodedField<Type>& rhs);

        //- Return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return PatchFunction1<Type>::Clone(*this);
        }

        //- Return a clone, setting the patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return PatchFunction1<Type>::Clone(*this, pp);
        }


    //- Destructor
    virtual ~CodedField() = default;


    // Member Functions

        //- Is value uniform (i.e. independent of coordinate)
        virtual inline bool uniform() const { return false; }


        // Evaluation

            //- Return CodedField value
            virtual tmp<Field<Type>> value(const scalar x) const;

            //- Integrate between two values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& pf1,
                const labelList& addr
            );


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
