/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::Function1Types::CodedFunction1

Description
    Function1 with the code supplied by an on-the-fly compiled C++
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
    CodedFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_CodedFunction1_H
#define Foam_Function1Types_CodedFunction1_H

#include "Function1.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                         Class CodedFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedFunction1
:
    public Function1<Type>,
    protected codedBase
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        const word redirectName_;

        mutable autoPtr<Function1<Type>> redirectFunctionPtr_;


    // Private Member Functions

        //- Get reference to the underlying Function1
        const Function1<Type>& redirectFunction() const;


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


    // Generated Methods

        //- No copy assignment
        void operator=(const CodedFunction1<Type>&) = delete;

public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedFunction1Template.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedFunction1Template.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        CodedFunction1
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit CodedFunction1(const CodedFunction1<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~CodedFunction1() = default;


    // Member Functions

        //- Return value at current time
        virtual inline Type value(const scalar x) const;

        // Integrate etc are not implemented!

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
