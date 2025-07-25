/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

Description
    Tests for regular expressions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "Switch.H"

#include "stringOps.H"
#include "SubStrings.H"
#include "regExpCxx.H"
#ifndef _WIN32
#include "regExpPosix.H"
#endif

using namespace Foam;


//- Simple structure for holding regex tests
struct regexTest
{
    bool expected;
    string pattern;
    string text;

    friend Istream& operator>>(Istream& is, struct regexTest& obj)
    {
        is.readBegin("struct");
        is >> obj.expected >> obj.pattern >> obj.text;
        is.readEnd("struct");
        is.check(FUNCTION_NAME);
        return is;
    }

    friend Ostream& operator<<(Ostream& os, const struct regexTest& obj)
    {
        os  << "( "
            << obj.expected << " " << obj.pattern << " " << obj.text
            << " )";
        return os;
    }
};


// Needed for list output. Just treat everything as unequal.
bool operator!=(const struct regexTest&, const struct regexTest&)
{
    return true;
}


// Simple output of match groups
static Ostream& operator<<(Ostream& os, const regExpCxx::results_type& sm)
{
    for (std::smatch::size_type i = 1; i < sm.size(); ++i)
    {
        os << " " << sm.str(i);
    }

    return os;
}


// Simple output of match groups
#ifndef _WIN32
static Ostream& operator<<(Ostream& os, const regExpPosix::results_type& sm)
{
    for (std::smatch::size_type i = 1; i < sm.size(); ++i)
    {
        os << " " << sm.str(i);
    }

    return os;
}
#endif


template<class RegexType>
void generalTests()
{
    Info<< nl << "test regExp(const char*) ..." << endl;
    string me("Mark");

    // Expect some failures:
    const bool oldThrowingError = FatalError.throwing(true);

    try
    {
        // Handling of null strings
        if (RegexType(nullptr).match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - null pointer is no expression" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Normal match
        if (RegexType("[Mm]ar[ck]").match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Match ignore case
        if (RegexType("mar[ck]", true).match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Embedded prefix for match ignore case
        if (RegexType("(?i)mar[ck]").match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Handling of empty expression
        if (RegexType("").match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - no match on empty expression" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Embedded prefix - but expression is empty
        if (RegexType("(?i)").match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - no match on empty expression" << endl;
        }
    }
    catch (const Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    FatalError.throwing(oldThrowingError);
}


template<class RegexType>
void testExpressions(const UList<regexTest>& tests)
{
    typename RegexType::results_type match;

    // Expect some failures:
    const bool oldThrowingError = FatalError.throwing(true);

    // Report matches:
    for (const auto& testseq : tests)
    {
        const bool expected = testseq.expected;
        const string& pat = testseq.pattern;
        const string& str = testseq.text;

        Info<< "Test " << Switch(expected) << ": "
            << str << " =~ m/" << pat.c_str() << "/ == ";

        RegexType re;

        try
        {
            re = pat;

            if (re.match(str, match))
            {
                Info<< "true";
                if (re.ngroups())
                {
                    Info<< " (" << re.ngroups() << " groups):" << match;
                }
            }
            else if (re.search(str))
            {
                Info<< "partial";
            }
            else
            {
                Info<< "false";
            }
            if (re.negated())
            {
                Info<< " (negated)";
            }
            Info<< endl;
        }
        catch (const Foam::error& err)
        {
            Info<< "Caught FatalError " << err << nl << endl;
            continue;
        }

        if (false)
        {
            RegexType re2(std::move(re));
            Info<<"move construct: " << re.exists() << "/" << re2.exists()
                << endl;

            re = std::move(re2);
            Info<<"move assign: " << re.exists() << "/" << re2.exists()
                << endl;

            re.swap(re2);
            Info<<"swap: " << re.exists() << "/" << re2.exists()
                << endl;
        }
    }

    FatalError.throwing(oldThrowingError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noParallel();

    argList::addBoolOption("cxx", "Test C++11 regular expressions");
    #ifndef _WIN32
    argList::addBoolOption("posix", "Test POSIX regular expressions");
    #endif
    argList::addOption
    (
        "regex",
        "expression",
        "regular expression to test"
    );

    argList::addArgument("file");
    argList::addArgument("...");
    argList::addArgument("fileN");
    argList::noMandatoryArgs();

    #include "setRootCase.H"

    // Newer compilers support regex directly
    #ifdef _GLIBCXX_RELEASE
    Info<< "_GLIBCXX_RELEASE = " << (_GLIBCXX_RELEASE) << nl;
    #endif

    if (std::is_same<regExp, regExpCxx>::value)
    {
        Info<< "Foam::regExp uses C++11 regex" << nl;
    }
    #ifndef _WIN32
    if (std::is_same<regExp, regExpPosix>::value)
    {
        Info<< "Foam::regExp uses POSIX regex" << nl;
    }
    #endif

    Info<< "sizeof std::regex:  " << sizeof(std::regex) << nl;
    Info<< "sizeof regex C++11: " << sizeof(regExpCxx) << nl;
    #ifndef _WIN32
    Info<< "sizeof regex POSIX: " << sizeof(regExpPosix) << nl;
    #endif
    Info<< "sizeof word: " << sizeof(Foam::word) << nl;
    Info<< "sizeof wordRe: " << sizeof(Foam::wordRe) << nl;
    Info<< "sizeof keyType: " << sizeof(Foam::keyType) << nl;

    if (!args.count({"cxx", "posix"}))
    {
        args.setOption("cxx");
        Info<< "Assuming -cxx as default" << nl;
    }
    Info<< nl;

    if (args.found("regex"))
    {
        std::string expr(args["regex"]);
        Info<< "regex: " << expr << nl;

        Info<< "(cxx)" << nl
            << "meta : " << Switch(regExpCxx::is_meta(expr)) << nl
            << "quotemeta: "
            << stringOps::quotemeta(expr, regExpCxx::meta()) << nl
            << nl;

        #ifndef _WIN32
        Info<< "(posix):" << nl
            << "meta : " << Switch(regExpPosix::is_meta(expr)) << nl
            << "quotemeta: "
            << stringOps::quotemeta(expr, regExpPosix::meta()) << nl
            << nl;
        #endif
        Info<< nl;
    }
    else if (args.size() < 2)
    {
        Info<< "No test files specified .. restrict to general tests" << nl;

        if (args.found("cxx"))
        {
            generalTests<regExpCxx>();
        }

        #ifndef _WIN32
        if (args.found("posix"))
        {
            generalTests<regExpPosix>();
        }
        #endif
    }

    for (label argi = 1; argi < args.size(); ++argi)
    {
        IFstream is(args.get<fileName>(argi));
        List<regexTest> tests(is);

        Info<< "Test expressions:" << tests << endl;
        IOobject::writeDivider(Info) << endl;

        if (args.found("cxx"))
        {
            testExpressions<regExpCxx>(tests);
        }

        #ifndef _WIN32
        if (args.found("posix"))
        {
            testExpressions<regExpPosix>(tests);
        }
        #endif
    }

    Info<< "\nDone" << nl << endl;

    return 0;
}


// ************************************************************************* //
