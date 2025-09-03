/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ifeqEntry.H"
#include "ifEntry.H"
#include "stringOps.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifeqEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifeqEntry,
        execute,
        dictionaryIstream,
        ifeq
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::ifeqEntry::readToken(token& t, Istream& is)
{
    // Skip dummy tokens - avoids entry::getKeyword consuming #else, #endif
    do
    {
        if (is.read(t).bad() || is.eof() || !t.good())
        {
            return;
        }
    }
    while (t == token::END_STATEMENT);
}


Foam::token Foam::functionEntries::ifeqEntry::expandToken
(
    const dictionary& dict,
    const string& keyword,
    const token& t
)
{
    if (keyword[0] == '$')
    {
        const word varName(keyword.substr(1));

        // Lookup the variable name in the given dictionary
        const entry* ePtr = dict.findScoped(varName, keyType::REGEX_RECURSIVE);
        if (ePtr)
        {
            return token(ePtr->stream());
        }
        else
        {
            // String expansion. Allow unset variables
            string expanded(keyword);
            stringOps::inplaceExpand(expanded, dict, true, true);

            // Re-form as a string token so we can compare to string
            return token(expanded, t.lineNumber());
        }
    }
    else if (!t.isString())
    {
        // Re-form as a string token so we can compare to string
        return token(keyword, t.lineNumber());
    }

    return t;
}


Foam::token Foam::functionEntries::ifeqEntry::expandToken
(
    const dictionary& dict,
    const token& t
)
{
    if (t.isWord())
    {
        return expandToken(dict, t.wordToken(), t);
    }
    else if (t.isVariable())
    {
        return expandToken(dict, t.stringToken(), t);
    }
    else if (t.isString())
    {
        return expandToken(dict, t.stringToken(), t);
    }

    return t;
}


bool Foam::functionEntries::ifeqEntry::equalToken
(
    const token& t1,
    const token& t2
)
{
    const bool eqType = (t1.type() == t2.type());

    switch (t1.type())
    {
        case token::UNDEFINED:
            return eqType;

        case token::BOOL:
            return (eqType && t1.boolToken() == t2.boolToken());

        case token::FLAG:
            return (eqType && t1.flagToken() == t2.flagToken());

        case token::PUNCTUATION:
            return (eqType && t1.pToken() == t2.pToken());

        case token::WORD:
        case token::DIRECTIVE:
            if (t2.isWord())
            {
                return t1.wordToken() == t2.wordToken();
            }
            else if (t2.isString())
            {
                const wordRe w2(t2.stringToken(), wordRe::DETECT);
                return w2.match(t1.wordToken());
            }
            return false;

        case token::STRING:
            if (eqType)
            {
                const wordRe w1(t1.stringToken(), wordRe::DETECT);
                const wordRe w2(t2.stringToken(), wordRe::DETECT);
                return w1.match(w2) || w2.match(w1);
            }
            else if (t2.isWord())
            {
                const wordRe w1(t1.stringToken(), wordRe::DETECT);
                return w1.match(t2.wordToken());
            }
            return false;

        case token::VARIABLE:
        case token::VERBATIM:
        case token::CHAR_DATA:
            if (t2.isStringType())
            {
                return t1.stringToken() == t2.stringToken();
            }
            return false;

        case token::LABEL:
            if (eqType)
            {
                return t1.labelToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.labelToken() == t2.scalarToken();
            }
            return false;

        case token::FLOAT:
            if (eqType)
            {
                return equal(t1.floatToken(), t2.floatToken());
            }
            else if (t2.isLabel())
            {
                return t1.floatToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.scalarToken() == t2.scalarToken();
            }
            return false;

        case token::DOUBLE:
            if (eqType)
            {
                return equal(t1.doubleToken(), t2.doubleToken());
            }
            else if (t2.isLabel())
            {
                return t1.doubleToken() == t2.labelToken();
            }
            else if (t2.isScalar())
            {
                return t1.scalarToken() == t2.scalarToken();
            }
            return false;

        case token::EXPRESSION:
            return false;

        case token::COMPOUND:
            return false;

        case token::ERROR:
            return eqType;
    }

    return false;
}


void Foam::functionEntries::ifeqEntry::skipUntil
(
    DynamicList<filePos>& stack,
    const dictionary& parentDict,
    const word& endDirective,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);

        if (!t.isDirective())
        {
            continue;
        }
        else if (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
        {
            stack.push_back(filePos(is.name(), is.lineNumber()));
            skipUntil(stack, parentDict, "#endif", is);
            stack.pop_back();
        }
        else if (t.wordToken() == endDirective)
        {
            return;
        }
    }

    FatalIOErrorInFunction(parentDict)
        << "Did not find matching " << endDirective << nl
        << exit(FatalIOError);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionEntries::ifeqEntry::evaluate
(
    const bool doIf,
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);
        bool pending = false;

        if (t.isDirective())
        {
            if (t.wordToken() == "#ifeq")
            {
                // Recurse to evaluate
                execute(stack, parentDict, is);
            }
            else if (t.wordToken() == "#if")
            {
                // Recurse to evaluate
                ifEntry::execute(stack, parentDict, is);
            }
            else if
            (
                doIf
             && (t.wordToken() == "#else" || t.wordToken() == "#elif")
            )
            {
                // Now skip until #endif
                skipUntil(stack, parentDict, "#endif", is);
                stack.pop_back();
                break;
            }
            else if (t.wordToken() == "#endif")
            {
                stack.pop_back();
                break;
            }
            else
            {
                pending = true;
            }
        }
        else
        {
            pending = true;
        }

        if (pending)
        {
            is.putBack(t);
            bool ok = entry::New(parentDict, is);
            if (!ok)
            {
                return false;
            }
        }
    }
    return true;
}


bool Foam::functionEntries::ifeqEntry::execute
(
    const bool doIf,
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    if (doIf)
    {
        evaluate(true, stack, parentDict, is);
    }
    else
    {
        // Fast-forward to #else
        token t;
        while (!is.eof())
        {
            readToken(t, is);

            // Only consider directives
            if (!t.isDirective())
            {
                continue;
            }

            if (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
            {
                stack.push_back(filePos(is.name(), is.lineNumber()));
                skipUntil(stack, parentDict, "#endif", is);
                stack.pop_back();
            }
            else if (t.wordToken() == "#else")
            {
                break;
            }
            else if (t.wordToken() == "#elif")
            {
                // const label lineNo = is.lineNumber();

                // Read line
                string line;
                dynamic_cast<ISstream&>(is).getLine(line);
                line += ';';
                IStringStream lineStream(line);
                const primitiveEntry e("ifEntry", parentDict, lineStream);

                if (ifEntry::isTrue(e.stream()))
                {
                    // Info<< "Using #elif " << doIf << " - line " << lineNo
                    //     << " in file " << is.relativeName() << endl;
                    break;
                }
            }
            else if (t.wordToken() == "#endif")
            {
                stack.pop_back();
                break;
            }
        }

        if (t.wordToken() == "#else")
        {
            // Evaluate until we hit #endif
            evaluate(false, stack, parentDict, is);
        }
        else if (t.wordToken() == "#elif")
        {
            // Evaluate until we hit #else or #endif
            evaluate(true, stack, parentDict, is);
        }
    }
    return true;
}


bool Foam::functionEntries::ifeqEntry::execute
(
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    const label nNested = stack.size();

    stack.push_back(filePos(is.name(), is.lineNumber()));

    // Read first token and expand any string
    token cond1(is);
    cond1 = expandToken(parentDict, cond1);

    // Read second token and expand any string
    token cond2(is);
    cond2 = expandToken(parentDict, cond2);

    const bool equal = equalToken(cond1, cond2);

    // Info<< "Using #" << typeName << " " << cond1
    //     << " == " << cond2
    //     << " at line " << stack.back().second()
    //     << " in file " <<  stack.back().first() << endl;

    bool ok = ifeqEntry::execute(equal, stack, parentDict, is);

    if (stack.size() != nNested)
    {
        FatalIOErrorInFunction(parentDict)
            << "Did not find matching #endif for condition starting"
            << " at line " << stack.back().second()
            << " in file " <<  stack.back().first() << exit(FatalIOError);
    }

    return ok;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifeqEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    return execute(stack, parentDict, is);
}


// ************************************************************************* //
