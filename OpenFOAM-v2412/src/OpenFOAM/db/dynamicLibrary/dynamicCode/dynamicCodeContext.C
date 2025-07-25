/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "dynamicCodeContext.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::dynamicCodeContext::inplaceExpand
(
    string& str,
    const dictionary& dict
)
{
    stringOps::inplaceTrim(str);
    stringOps::inplaceExpand(str, dict);
}


unsigned Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    label lineNum,
    const string& file
)
{
    ++lineNum;  // Change from 0-based to 1-based

    const auto len = code.length();

    if (lineNum > 0 && len && !file.empty())
    {
        code = "#line " + Foam::name(lineNum) + " \"" + file + "\"\n" + code;

        return (code.length() - len);
    }

    return 0;
}


unsigned Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    label lineNum,
    const dictionary& dict
)
{
    return addLineDirective(code, lineNum, dict.name());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCodeContext::dynamicCodeContext()
:
    dict_(std::cref<dictionary>(dictionary::null))
{}


Foam::dynamicCodeContext::dynamicCodeContext(const dictionary& dict)
:
    dynamicCodeContext()
{
    setCodeContext(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::dynamicCodeContext::good() const noexcept
{
    return &(dict_.get()) != &(dictionary::null);
}


const Foam::entry* Foam::dynamicCodeContext::findEntry(const word& key) const
{
    return this->dict().findEntry(key, keyType::LITERAL);
}


bool Foam::dynamicCodeContext::readEntry
(
    const word& key,
    string& str,
    bool mandatory,
    bool withLineNum
)
{
    str.clear();
    sha1_.append("<" + key + ">");

    const dictionary& dict = this->dict();
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);

    if (!eptr)
    {
        if (mandatory)
        {
            FatalIOErrorInFunction(dict)
                << "Entry '" << key << "' not found in dictionary "
                << dict.name() << nl
                << exit(FatalIOError);
        }

        return false;
    }

    // Expand dictionary entries.
    // Removing any leading/trailing whitespace is necessary for compilation
    // options, but is also convenient for includes and code body.

    eptr->readEntry(str);
    dynamicCodeContext::inplaceExpand(str, dict);
    sha1_.append(str);

    if (withLineNum)
    {
        addLineDirective(str, eptr->startLineNumber(), dict);
    }

    return true;
}


bool Foam::dynamicCodeContext::readIfPresent
(
    const word& key,
    string& str,
    bool withLineNum
)
{
    return readEntry(key, str, false, withLineNum);
}


void Foam::dynamicCodeContext::setCodeContext(const dictionary& dict)
{
    dict_ = std::cref<dictionary>(dict);
    sha1_.clear();

    // No #line for options (Make/options)
    readIfPresent("codeOptions", codeOptions_, false);

    // No #line for libs (LIB_LIBS)
    readIfPresent("codeLibs", codeLibs_, false);

    readIfPresent("codeInclude", codeInclude_);
    readIfPresent("localCode", localCode_);
    readIfPresent("code", code_);
}


// ************************************************************************* //
