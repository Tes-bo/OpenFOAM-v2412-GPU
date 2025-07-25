/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Test word/regex

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "ITstream.H"
#include "List.H"
#include "Tuple2.H"
#include "keyType.H"
#include "wordRes.H"
#include "predicates.H"
#include "Random.H"

using namespace Foam;


word typeOf(wordRe::compOption retval)
{
    if (wordRe::LITERAL == retval) return "(literal)";
    if (wordRe::UNKNOWN == retval) return "(unknown)";
    if (wordRe::REGEX == retval)   return "(regex)";
    return "";
}


Ostream& printInfo(const wordRe& wre)
{
    if (wre.isPattern())
    {
        Info<< "wordRe(regex) ";
    }
    else
    {
        Info<< "wordRe(plain) ";
    }
    Info<< wre;
    return Info;
}


// Could use something like this for reading wordRes
void exptl_reading(Istream& is, wordRes& list)
{
    token tok(is);

    bool ok = ((tok.isWord() || tok.isQuotedString()) && !tok.isCompound());
    if (ok)
    {
        list.resize(1);
        ok = list[0].assign(tok);
    }
    if (!ok)
    {
        if (tok.good())
        {
            is.putBack(tok);
        }
        list.readList(is);
    }
}


bool testReadList_wordRes(const std::string& input)
{
    ITstream is(input);
    wordRes list;

    exptl_reading(is, list);

    const label nTrailing = is.nRemainingTokens();

    Info<< "input:<<<<" << nl << input.c_str() << nl
        << ">>>> with " << nTrailing << " tokens remaining" << nl
        << "list: " << flatOutput(list) << nl;

    return !nTrailing;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    if (true)
    {
        Info<< "Test input for wordRes" << nl << nl;

        testReadList_wordRes
        (
            "("
            "this \"x.*\" \"file[a-b]\" xvalues \"xv.*\""
            ")"
        );

        testReadList_wordRes
        (
            "\".*\""
        );

        Info<< nl << nl;
    }

    wordRe wre;
    std::string s1("this .* file");
    Foam::string s2("this .* file");
    const char * s3 = "this .* file";

    keyType keyre("x.*", keyType::REGEX);

    wordRes wordrelist
    ({
        {"this", wordRe::LITERAL},
        {"x.*", wordRe::REGEX},
        {"file[a-b]", wordRe::REGEX},
        {"xvalues", wordRe::LITERAL},
        {"xv.*", wordRe::REGEX},
    });

    if (false)
    {
        Info<<"keyType: " << keyre << nl;

        keyType key2(std::move(keyre));

        Info<<"move construct: <" << keyre << "> <" << key2 << ">" << nl;

        keyre = std::move(key2);

        Info<<"move assign: <" << keyre << "> <" << key2 << ">" << nl;

        keyType key3;

        keyre.swap(key3);

        Info<<"swap: <" << keyre << "> <" << key3 << ">" << nl;

        keyre = std::move(key3);
        Info<<"move assign: <" << keyre << "> <" << key3 << ">" << nl;

        return 0;
    }

    if (false)
    {
        wordRe keyre("y.*", wordRe::REGEX);

        Info<<"wordRe: " << keyre << nl;

        wordRe key2(std::move(keyre));

        Info<<"keyTypes: " << keyre << " " << key2 << nl;

        keyre = std::move(key2);

        Info<<"keyTypes: " << keyre << " " << key2 << nl;

        wordRe key3;

        keyre.swap(key3);

        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << nl;

        keyre = std::move(key3);
        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << nl;

        return 0;
    }

    wordRes wres1(wordrelist);

    Info<< "re-list:" << wres1 << nl;
    Info<< "match this: " << wres1("this") << nl;
    Info<< "match xyz: "  << wres1("xyz") << nl;
    Info<< "match zyx: "  << wres1("zyx") << nl;
    Info<< "match xyz: "  << wres1.match("xyz") << nl;
    Info<< "match any: "  << predicates::always()("any junk") << nl;

    Info<< "match xvalues: "  << wres1.match("xvalues") << nl;
    Info<< "matched xvalues = "  << typeOf(wres1.matched("xvalues")) << nl;
    Info<< "matched xval = "  << typeOf(wres1.matched("xval")) << nl;
    Info<< "matched zyx = "  << typeOf(wres1.matched("zyx")) << nl;

    Info<< "keyre match: "  << keyre("xyz") << nl;
    Info<< "string match: "  << string("this").match("xyz") << nl;
    Info<< "string match: "  << string("x.*")("xyz") << nl;
    Info<< "string match: "  << string("x.*")(keyre) << nl;


    // Test uniq
    {
        Random rnd;
        const label last = wres1.size()-1;

        for (label i = 0; i < 8; ++i)
        {
            // Make a copy
            wordRe wre(wres1[rnd.position<label>(0,last)]);

            // Append
            wres1.push_back(wre);
        }

        // Add some entropy
        Foam::shuffle(wres1);

        Info<< nl
            << "Test uniq on " << wres1
            << "  ==  " << wordRes::uniq(wres1) << nl;

        // Inplace
        wres1.uniq();
        Info<< nl << "Inplace: " << wres1 << nl;
    }
    Info<< nl;

    printInfo(wordRe(s1, wordRe::DETECT)) << nl;
    printInfo(wordRe(s2)) << nl;
    printInfo(wordRe(s2, wordRe::DETECT)) << nl;
    printInfo(wordRe(s3, wordRe::REGEX)) << nl;

    wre = "this .* file";

    Info<<"substring: " << wre.substr(4) << nl;

    printInfo(wre) << nl;
    wre = s1;
    printInfo(wre) << nl;
    wre.uncompile();
    printInfo(wre) << nl;

    wre = "something";
    printInfo(wre) << " before" << nl;
    wre.uncompile();
    printInfo(wre) << " uncompiled" << nl;
    wre.compile(wordRe::DETECT);
    printInfo(wre) << " after DETECT" << nl;
    wre.compile(wordRe::ICASE);
    printInfo(wre) << " after ICASE" << nl;
    wre.compile(wordRe::DETECT_ICASE);
    printInfo(wre) << " after DETECT_ICASE" << nl;

    wre = "something .* value";
    printInfo(wre) << " before" << nl;
    wre.uncompile();
    printInfo(wre) << " uncompiled" << nl;
    wre.compile(wordRe::DETECT);
    printInfo(wre) << " after DETECT" << nl;
    wre.uncompile();
    printInfo(wre) << " uncompiled" << nl;
    wre.compile();
    printInfo(wre) << " re-compiled" << nl;

    wre.set("something .* value", wordRe::LITERAL);
    printInfo(wre) << " set as LITERAL" << nl;

    IOobject::writeDivider(Info);

    List<Tuple2<wordRe, string>> rawList(IFstream("testRegexps")());
    Info<< "input list:" << rawList << nl;
    IOobject::writeDivider(Info) << nl;

    forAll(rawList, elemI)
    {
        const wordRe& wre = rawList[elemI].first();
        const string& str = rawList[elemI].second();

        printInfo(wre)
            << " equals:" << (wre == str)
            << "(" << wre.match(str, true) << ")"
            << " match:" << wre.match(str)
            << "  str=" << str
            << nl;

        wordRe wre2;
        wre2.set(wre, wordRe::ICASE);

        printInfo(wre2)
            << " match:" << wre2.match(str)
            << "  str=" << str
            << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
