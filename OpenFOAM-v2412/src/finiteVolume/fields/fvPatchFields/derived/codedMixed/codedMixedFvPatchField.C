/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "codedMixedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "dynamicCode.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::dlLibraryTable& Foam::codedMixedFvPatchField<Type>::libs() const
{
    return this->db().time().libs();
}


template<class Type>
Foam::string Foam::codedMixedFvPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::clearRedirect() const
{
    redirectPatchFieldPtr_.reset(nullptr);
}


template<class Type>
const Foam::dictionary&
Foam::codedMixedFvPatchField<Type>::codeContext() const
{
    const dictionary* ptr = dict_.findDict("codeContext", keyType::LITERAL);
    return (ptr ? *ptr : dictionary::null);
}


template<class Type>
const Foam::dictionary&
Foam::codedMixedFvPatchField<Type>::codeDict() const
{
    // Inline "code" or from system/codeDict
    return
    (
        dict_.found("code")
      ? dict_
      : codedBase::codeDict(this->db()).subDict(name_)
    );
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Take no chances - typeName must be identical to name_
    dynCode.setFilterVariable("typeName", name_);

    // Set TemplateType and FieldType filter variables
    dynCode.setFieldTemplates<Type>();

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    #ifdef FULLDEBUG
    dynCode.setFilterVariable("verbose", "true");
    DetailInfo
        <<"compile " << name_ << " sha1: " << context.sha1() << endl;
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
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    codedBase(),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& rhs,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    codedBase(),
    dict_(rhs.dict_),  // Deep copy
    name_(rhs.name_),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    codedBase(),
    dict_
    (
        // Copy dictionary, but without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value", "refValue", "refGradient", "valueFraction"
            })
        )
    ),
    name_(dict.getCompat<word>("name", {{"redirectType", 1706}})),
    redirectPatchFieldPtr_(nullptr)
{
    updateLibrary(name_);
}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& rhs
)
:
    parent_bctype(rhs),
    codedBase(),
    dict_(rhs.dict_),  // Deep copy
    name_(rhs.name_),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& rhs,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(rhs, iF),
    codedBase(),
    dict_(rhs.dict_),  // Deep copy
    name_(rhs.name_),
    redirectPatchFieldPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mixedFvPatchField<Type>&
Foam::codedMixedFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_)
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        // Write the data from the mixed b.c.
        OCharStream os;
        this->parent_bctype::write(os);
        ISpanStream is(os.view());
        dictionary constructDict(is);

        // Override type
        constructDict.set("type", name_);

        redirectPatchFieldPtr_.reset
        (
            dynamic_cast<parent_bctype*>
            (
                fvPatchField<Type>::New
                (
                    this->patch(),
                    this->internalField(),
                    constructDict
                ).ptr()
            )
        );

        // Forward copy of dictionary content to the code template
        auto* contentPtr =
            dynamic_cast<dictionaryContent*>(redirectPatchFieldPtr_.get());

        if (contentPtr)
        {
            contentPtr->dict(this->codeContext());
        }
        else
        {
            WarningInFunction
                << name_ << " Did not derive from dictionaryContent"
                << nl << nl;
        }
    }
    return *redirectPatchFieldPtr_;
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const parent_bctype& fvp = redirectPatchField();
    const_cast<parent_bctype&>(fvp).updateCoeffs();

    // Copy through coefficients
    this->refValue() = fvp.refValue();
    this->refGrad() = fvp.refGrad();
    this->valueFraction() = fvp.valueFraction();

    this->parent_bctype::updateCoeffs();
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const parent_bctype& fvp = redirectPatchField();

    // - updates the value of fvp (though not used)
    // - resets the updated() flag
    const_cast<parent_bctype&>(fvp).evaluate(commsType);

    // Update the value (using the coefficients) locally
    parent_bctype::evaluate(commsType);
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    os.writeEntry("name", name_);

    codedBase::writeCodeDict(os, dict_);
}


// ************************************************************************* //
