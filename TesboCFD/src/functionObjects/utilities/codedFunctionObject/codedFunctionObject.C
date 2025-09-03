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

\*---------------------------------------------------------------------------*/

#include "codedFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "dictionaryContent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(codedFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        codedFunctionObject,
        dictionary
    );
} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable& Foam::functionObjects::codedFunctionObject::libs() const
{
    return time_.libs();
}


Foam::string Foam::functionObjects::codedFunctionObject::description() const
{
    return "functionObject " + name();
}


void Foam::functionObjects::codedFunctionObject::clearRedirect() const
{
    redirectFunctionObjectPtr_.reset(nullptr);
}


const Foam::dictionary&
Foam::functionObjects::codedFunctionObject::codeContext() const
{
    const dictionary* ptr = dict_.findDict("codeContext", keyType::LITERAL);
    return (ptr ? *ptr : dictionary::null);
}


const Foam::dictionary&
Foam::functionObjects::codedFunctionObject::codeDict() const
{
    return dict_;
}


void Foam::functionObjects::codedFunctionObject::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("codeData", codeData_);
    dynCode.setFilterVariable("codeRead", codeRead_);
    dynCode.setFilterVariable("codeExecute", codeExecute_);
    dynCode.setFilterVariable("codeWrite", codeWrite_);
    dynCode.setFilterVariable("codeEnd", codeEnd_);

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

Foam::functionObjects::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    codedBase(),
    dict_(dict)
{
    read(dict_);

    updateLibrary(name_);
    redirectFunctionObject();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject&
Foam::functionObjects::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_)
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            name_,
            time_,
            constructDict
        );


        // Forward copy of codeContext to the code template
        auto* contentPtr =
            dynamic_cast<dictionaryContent*>(redirectFunctionObjectPtr_.get());

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
    return *redirectFunctionObjectPtr_;
}


bool Foam::functionObjects::codedFunctionObject::execute()
{
    updateLibrary(name_);
    return redirectFunctionObject().execute();
}


bool Foam::functionObjects::codedFunctionObject::write()
{
    updateLibrary(name_);
    return redirectFunctionObject().write();
}


bool Foam::functionObjects::codedFunctionObject::end()
{
    updateLibrary(name_);
    return redirectFunctionObject().end();
}


bool Foam::functionObjects::codedFunctionObject::read(const dictionary& dict)
{
    timeFunctionObject::read(dict);

    codedBase::setCodeContext(dict);

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);

    auto& ctx = codedBase::codeContext();

    // Get code chunks, no short-circuiting
    int nKeywords = 0;
    nKeywords += ctx.readIfPresent("codeData", codeData_);
    nKeywords += ctx.readIfPresent("codeRead", codeRead_);
    nKeywords += ctx.readIfPresent("codeExecute", codeExecute_);
    nKeywords += ctx.readIfPresent("codeWrite", codeWrite_);
    nKeywords += ctx.readIfPresent("codeEnd", codeEnd_);

    if (!nKeywords)
    {
        IOWarningInFunction(dict)
            << "No critical \"code\" prefixed keywords found." << nl
            << "Please check the code documentation for more details." << nl
            << endl;
    }

    updateLibrary(name_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
