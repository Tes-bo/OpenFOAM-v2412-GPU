/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::dlLibraryTable&
Foam::Function1Types::CodedFunction1<Type>::libs() const
{
    return this->time().libs();
}


template<class Type>
Foam::string
Foam::Function1Types::CodedFunction1<Type>::description() const
{
    return "CodedFunction1 " + redirectName_;
}


template<class Type>
void Foam::Function1Types::CodedFunction1<Type>::clearRedirect() const
{
    redirectFunctionPtr_.reset(nullptr);
}


template<class Type>
const Foam::dictionary&
Foam::Function1Types::CodedFunction1<Type>::codeContext() const
{
    // What else would make sense?
    return dict_;
}


template<class Type>
const Foam::dictionary&
Foam::Function1Types::CodedFunction1<Type>::codeDict
(
    const dictionary& dict
) const
{
    // Use named subdictionary if present to provide the code.
    // This allows running with multiple Function1s

    return
    (
        dict.found("code")
      ? dict
      : dict.subDict(redirectName_)
    );
}


template<class Type>
const Foam::dictionary&
Foam::Function1Types::CodedFunction1<Type>::codeDict() const
{
    return codeDict(dict_);
}


template<class Type>
void Foam::Function1Types::CodedFunction1<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    if (context.code().empty())
    {
        FatalIOErrorInFunction(dict_)
            << "No code section in input dictionary for Function1 "
            << " name " << redirectName_
            << exit(FatalIOError);
    }

    // Take no chances - typeName must be identical to redirectName_
    dynCode.setFilterVariable("typeName", redirectName_);

    // Set TemplateType and FieldType filter variables
    dynCode.setFieldTemplates<Type>();

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    #ifdef FULLDEBUG
    dynCode.setFilterVariable("verbose", "true");
    DetailInfo
        <<"compile " << redirectName_ << " sha1: " << context.sha1() << endl;
    #endif

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lOpenFOAM \\\n"
        "    -lmeshTools \\\n"
      + context.libs()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::CodedFunction1<Type>::CodedFunction1
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    codedBase(),
    dict_(dict),
    redirectName_(dict.getOrDefault<word>("name", entryName))
{
    this->codedBase::setCodeContext(dict_);

    // No additional code chunks...

    updateLibrary(redirectName_);
}


template<class Type>
Foam::Function1Types::CodedFunction1<Type>::CodedFunction1
(
    const CodedFunction1<Type>& rhs
)
:
    Function1<Type>(rhs),
    codedBase(),
    dict_(rhs.dict_),
    redirectName_(rhs.redirectName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::Function1<Type>&
Foam::Function1Types::CodedFunction1<Type>::redirectFunction() const
{
    if (!redirectFunctionPtr_)
    {
        dictionary constructDict;
        // Force 'redirectName_' sub-dictionary into existence
        dictionary& coeffs = constructDict.subDictOrAdd(redirectName_);

        coeffs = dict_;  // Copy input code and coefficients
        coeffs.remove("name");      // Redundant
        coeffs.set("type", redirectName_);  // Specify our new (redirect) type

        redirectFunctionPtr_.reset
        (
            Function1<Type>::New
            (
                redirectName_,
                constructDict,
                this->whichDb()
            )
        );

        // Forward copy of codeContext to the code template
        auto* contentPtr =
            dynamic_cast<dictionaryContent*>(redirectFunctionPtr_.get());

        if (contentPtr)
        {
            contentPtr->dict(this->codeContext());
        }
        else
        {
            WarningInFunction
                << redirectName_ << " Did not derive from dictionaryContent"
                << nl << nl;
        }
    }
    return *redirectFunctionPtr_;
}


template<class Type>
Type Foam::Function1Types::CodedFunction1<Type>::value
(
    const scalar x
) const
{
    // Ensure library containing user-defined code is up-to-date
    updateLibrary(redirectName_);

    return redirectFunction().value(x);
}


template<class Type>
void Foam::Function1Types::CodedFunction1<Type>::writeData
(
    Ostream& os
) const
{
    // Should really only output only relevant entries but since using
    // Function1-from-subdict upon construction our dictionary contains
    // only the relevant entries.
    dict_.writeEntry(this->name(), os);
}


// ************************************************************************* //
