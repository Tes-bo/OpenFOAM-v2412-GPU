/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "codedFixedValuePointPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"
#include "dynamicCode.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::dlLibraryTable& Foam::codedFixedValuePointPatchField<Type>::libs() const
{
    return this->db().time().libs();
}


template<class Type>
Foam::string Foam::codedFixedValuePointPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::clearRedirect() const
{
    redirectPatchFieldPtr_.reset(nullptr);
}


template<class Type>
const Foam::dictionary&
Foam::codedFixedValuePointPatchField<Type>::codeContext() const
{
    const dictionary* ptr = dict_.findDict("codeContext", keyType::LITERAL);
    return (ptr ? *ptr : dictionary::null);
}


template<class Type>
const Foam::dictionary&
Foam::codedFixedValuePointPatchField<Type>::codeDict() const
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
void Foam::codedFixedValuePointPatchField<Type>::prepare
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
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    parent_bctype(p, iF),
    codedBase(),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& rhs,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    codedBase(),
    dict_(rhs.dict_),  // Deep copy
    name_(rhs.name_),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    // The 'value' is optional (handled below)
    parent_bctype(p, iF, dict, IOobjectOption::NO_READ),
    codedBase(),
    dict_
    (
        // Copy dictionary without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value"
            })
        )
    ),
    name_(dict.getCompat<word>("name", {{"redirectType", 1706}})),
    redirectPatchFieldPtr_(nullptr)
{
    updateLibrary(name_);

    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        // Evaluate to assign a value
        this->evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& rhs
)
:
    parent_bctype(rhs),
    codedBase(),
    dict_(rhs.dict_),  // Deep copy
    name_(rhs.name_),
    redirectPatchFieldPtr_(nullptr)
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& rhs,
    const DimensionedField<Type, pointMesh>& iF
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
const Foam::pointPatchField<Type>&
Foam::codedFixedValuePointPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_)
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OCharStream os;
        this->writeValueEntry(os);
        ISpanStream is(os.view());
        dictionary constructDict(is);

        constructDict.set("type", name_);

        redirectPatchFieldPtr_.reset
        (
            pointPatchField<Type>::New
            (
                this->patch(),
                this->internalField(),
                constructDict
            ).ptr()
        );


        // Forward copy of codeContext to the code template
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
void Foam::codedFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary(name_);

    const pointPatchField<Type>& fvp = redirectPatchField();

    const_cast<pointPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through value
    this->operator==(fvp);

    parent_bctype::updateCoeffs();
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary(name_);

    const pointPatchField<Type>& fvp = redirectPatchField();

    const_cast<pointPatchField<Type>&>(fvp).evaluate(commsType);

    parent_bctype::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    os.writeEntry("name", name_);

    codedBase::writeCodeDict(os, dict_);
}


// ************************************************************************* //
