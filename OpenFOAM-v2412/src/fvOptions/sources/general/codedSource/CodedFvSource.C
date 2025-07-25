/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "CodedFvSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class Type>
Foam::dlLibraryTable& Foam::fv::CodedSource<Type>::libs() const
{
    return mesh_.time().libs();
}


template<class Type>
Foam::string Foam::fv::CodedSource<Type>::description() const
{
    return "fvOption::" + name_;
}


template<class Type>
void Foam::fv::CodedSource<Type>::clearRedirect() const
{
    redirectOptionPtr_.reset(nullptr);
}


template<class Type>
const Foam::dictionary& Foam::fv::CodedSource<Type>::codeDict() const
{
    return coeffs_;
}


template<class Type>
void Foam::fv::CodedSource<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    const word sourceType(pTraits<Type>::typeName);

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("TemplateType", sourceType);
    dynCode.setFilterVariable("SourceType", sourceType + "Source");

    //dynCode.removeFilterVariable("code");
    dynCode.setFilterVariable("codeCorrect", codeCorrect_);
    dynCode.setFilterVariable("codeConstrain", codeConstrain_);

    // Missing codeAddSup or codeAddSupRho handled on input,
    // but also check for empty() and inject "NotImplemented" to avoid
    // cases where the user has specified or used the wrong one
    // (ie, compresssible vs incompressible)

    if (codeAddSup_.empty())
    {
        dynCode.setFilterVariable("codeAddSup", "NotImplemented");
    }
    else
    {
        dynCode.setFilterVariable("codeAddSup", codeAddSup_);
    }
    if (codeAddSupRho_.empty())
    {
        dynCode.setFilterVariable("codeAddSupRho", "NotImplemented");
    }
    else
    {
        dynCode.setFilterVariable("codeAddSupRho", codeAddSupRho_);
    }

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    #ifdef FULLDEBUG
    dynCode.setFilterVariable("verbose", "true");
    DetailInfo
        <<"compile " << name_ << " sha1: " << context.sha1() << endl;
    #endif

    // define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/fvOptions/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
        "-I$(LIB_SRC)/sampling/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lfiniteVolume \\\n"
        "    -lfvOptions \\\n"
        "    -lmeshTools \\\n"
        "    -lsampling \\\n"
      + context.libs()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::CodedSource<Type>::CodedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    codedBase::setCodeContext(coeffs_);

    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    coeffs_.readEntry("fields", fieldNames_);

    fv::option::resetApplied();

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);


    // Code context chunks

    auto& ctx = codedBase::codeContext();

    ctx.readEntry("codeCorrect", codeCorrect_);

    // Need one or both codeAddSup/codeAddSupRho.
    // - do precheck and let the usual mechanism fail

    const bool mandatory =
    (
        !ctx.findEntry("codeAddSup")
     && !ctx.findEntry("codeAddSupRho")
    );

    ctx.readEntry("codeAddSup", codeAddSup_, mandatory);
    ctx.readEntry("codeAddSupRho", codeAddSupRho_, mandatory);

    // ctx.readEntry("codeConstrain", codeConstrain_);
    ctx.readEntry  // Compatibility
    (
        coeffs_.lookupEntryCompat
        (
            "codeConstrain",
            {{ "codeSetValue", 1812 }},
            keyType::LITERAL
        ).keyword(),
        codeConstrain_
    );

    return true;
}


template<class Type>
Foam::fv::option& Foam::fv::CodedSource<Type>::redirectOption() const
{
    if (!redirectOptionPtr_)
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);
        constructDict.changeKeyword(modelType_ & "Coeffs", name_ & "Coeffs");

        redirectOptionPtr_ = fv::option::New
        (
            name_,
            constructDict,
            mesh_
        );
    }
    return *redirectOptionPtr_;
}


template<class Type>
void Foam::fv::CodedSource<Type>::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    DebugInfo
        << "fv::CodedSource<" << pTraits<Type>::typeName
        << ">::correct for source " << name_ << endl;

    updateLibrary(name_);
    redirectOption().correct(field);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "fv::CodedSource<" << pTraits<Type>::typeName
        << ">::addSup for source " << name_ << endl;

    updateLibrary(name_);
    redirectOption().addSup(eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "fv::CodedSource<" << pTraits<Type>::typeName
        << ">::addSup(rho) for source " << name_ << endl;

    updateLibrary(name_);
    redirectOption().addSup(rho, eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "fv::CodedSource<" << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    updateLibrary(name_);
    redirectOption().constrain(eqn, fieldi);
}


// ************************************************************************* //
