/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
Foam::PatchFunction1Types::CodedField<Type>::libs() const
{
    return this->patch_.boundaryMesh().mesh().time().libs();
}


template<class Type>
Foam::string
Foam::PatchFunction1Types::CodedField<Type>::description() const
{
    return "CodedField " + redirectName_;
}


template<class Type>
void Foam::PatchFunction1Types::CodedField<Type>::clearRedirect() const
{
    redirectFunctionPtr_.reset(nullptr);
}


template<class Type>
const Foam::dictionary&
Foam::PatchFunction1Types::CodedField<Type>::codeContext() const
{
    // What else would make sense?
    return dict_;
}


template<class Type>
const Foam::dictionary&
Foam::PatchFunction1Types::CodedField<Type>::codeDict
(
    const dictionary& dict
) const
{
    // Use named subdictionary if present to provide the code.
    // This allows running with multiple PatchFunction1s

    return
    (
        dict.found("code")
      ? dict
      : dict.subDict(redirectName_)
    );
}


template<class Type>
const Foam::dictionary&
Foam::PatchFunction1Types::CodedField<Type>::codeDict() const
{
    return codeDict(dict_);
}


template<class Type>
void Foam::PatchFunction1Types::CodedField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    if (context.code().empty())
    {
        FatalIOErrorInFunction(dict_)
            << "No code section in input dictionary for patch "
            << this->patch_.name()
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
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lOpenFOAM \\\n"
        "    -lfiniteVolume \\\n"
        "    -lmeshTools \\\n"
      + context.libs()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::CodedField<Type>::CodedField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    codedBase(),
    dict_(dict),
    redirectName_(dict.getOrDefault<word>("name", entryName))
{
    this->codedBase::setCodeContext(dict_);

    // No additional code chunks...

    updateLibrary(redirectName_);
}


template<class Type>
Foam::PatchFunction1Types::CodedField<Type>::CodedField
(
    const CodedField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    codedBase(),
    dict_(rhs.dict_),
    redirectName_(rhs.redirectName_)
{}


template<class Type>
Foam::PatchFunction1Types::CodedField<Type>::CodedField
(
    const CodedField<Type>& rhs
)
:
    CodedField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::PatchFunction1<Type>&
Foam::PatchFunction1Types::CodedField<Type>::redirectFunction() const
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
            PatchFunction1<Type>::New
            (
                this->patch(),
                redirectName_,
                constructDict,
                this->faceValues()
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
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::CodedField<Type>::value
(
    const scalar x
) const
{
    // Ensure library containing user-defined code is up-to-date
    updateLibrary(redirectName_);

    return redirectFunction().value(x);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::CodedField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // Ensure library containing user-defined code is up-to-date
    updateLibrary(redirectName_);

    return redirectFunction().integrate(x1, x2);
}


template<class Type>
void Foam::PatchFunction1Types::CodedField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    PatchFunction1<Type>::autoMap(mapper);
    if (redirectFunctionPtr_)
    {
        redirectFunctionPtr_->autoMap(mapper);
    }
}


template<class Type>
void Foam::PatchFunction1Types::CodedField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    PatchFunction1<Type>::rmap(pf1, addr);
    if (redirectFunctionPtr_)
    {
        redirectFunctionPtr_->rmap(pf1, addr);
    }
}


template<class Type>
void Foam::PatchFunction1Types::CodedField<Type>::writeData
(
    Ostream& os
) const
{
    // Should really only output only relevant entries but since using
    // PatchFunction1-from-subdict upon construction our dictionary contains
    // only the relevant entries. It would be different if PatchFunction1-from
    // primitiveEntry when the whole 'value' entry would be present
    dict_.writeEntry(this->name(), os);
}


// ************************************************************************* //
