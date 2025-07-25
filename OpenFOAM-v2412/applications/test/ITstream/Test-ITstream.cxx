/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Basic test for ITstream, which is the input token stream used by
    primitiveEntry (dictionary).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"
#include "ITstream.H"
#include "ListOps.H"
#include "stringListOps.H"

using namespace Foam;

template<class T>
Ostream& toString(Ostream& os, const T& str)
{
    os << str;
    return os;
}


template<>
Ostream& toString(Ostream& os, const UList<char>& list)
{
    for (const char c : list)
    {
        os << c;
    }

    return os;
}


template<>
Ostream& toString(Ostream& os, const List<char>& list)
{
    for (const char c : list)
    {
        os << c;
    }

    return os;
}


void printTokens(Istream& is)
{
    label count = 0;
    token t;
    while (is.good())
    {
        is >> t;
        if (t.good())
        {
            ++count;
            Info<< "token: " << t << endl;
        }
    }

    Info<< count << " tokens" << endl;
}


Ostream& reportPeek(const ITstream& is)
{
    Info<< "  index : " << is.tokenIndex() << nl
        << "  peek  : " << is.peek().info() << nl;
    return Info;
}


template<class BUF>
void doTest
(
    const string& name,
    const BUF& input,
    bool verbose = false,
    bool testskip = false
)
{
    Info<< "test " << name.c_str() << ":" << nl
        << "====" << nl;
    toString(Info, input)
        << nl
        << "====" << nl << endl;

    ITstream its(input);
    Info<< "got " << its.size() << " tokens - index at "
        << its.tokenIndex() << endl;

    if (verbose)
    {
        for (const token& tok : its)
        {
            Info<< "    " << tok.info() << nl;
        }
        Info<< nl;
    }

    if (testskip)
    {
        Info<< " front : " << its.front().info() << nl
            << " back  : " << its.back().info() << nl;

        Info<< "rewind():" << nl;
        reportPeek(its);

        its.skip(3);
        Info<< "skip(3):" << nl;
        reportPeek(its);

        its.skip(2);
        Info<< "skip(2):" << nl;
        reportPeek(its);

        its.skip(-2);
        Info<< "skip(-2):" << nl;
        reportPeek(its);

        its.skip(100);
        Info<< "skip(100):" << nl;
        reportPeek(its);

        its.skip(-1000);
        Info<< "skip(-1000):" << nl;
        reportPeek(its);
    }
}


void printToken(const label index, const token& tok)
{
    Info<< "  " << index << "  " << tok.name();
    if (tok.good())
    {
        Info<< " : " << tok;
    }
    Info<< nl;
}


template<class BUF>
void testWalk1
(
    const std::string& name,
    const BUF& input,
    const int verbose
)
{
    Info<< "tokenized " << name.c_str() << ":" << nl
        << "====" << nl;
    toString(Info, input)
        << nl
        << "====" << endl;

    ITstream is(input);
    Info<< is.size() << " tokens" << endl;
    for (is.rewind(); !is.eof(); is.skip())
    {
        printToken(is.tokenIndex(), is.currentToken());
    }
    Info<< nl;

    Info<< "every other token:" << nl;
    for (is.seek(1); is.nRemainingTokens(); is.skip(2))
    {
        printToken(is.tokenIndex(), is.currentToken());
    }

    for (int i : { 3, 7, 11, 20 })
    {
        Info<< "peekToken: ";
        printToken(i, is.peekToken(i));
    }

    labelRange range(is.size()-2, 2);
    Info<< nl
        << "remove: " << range << " of 0/" << is.size() << " tokens" << endl;
    is.remove(range);

    Info<< "Now " << is.size() << " tokens" << endl;
    for (is.rewind(); !is.eof(); is.skip())
    {
        printToken(is.tokenIndex(), is.currentToken());
    }

    range.reset(10, 3);
    Info<< nl
        << "remove: " << range << " of 0/" << is.size() << " tokens" << endl;
    is.remove(range);

    Info<< "Now " << is.size() << " tokens" << endl;
    for (is.rewind(); !is.eof(); is.skip())
    {
        printToken(is.tokenIndex(), is.currentToken());
    }

    Info<< nl;
}


void testRewrite(const std::string& input, const int verbose)
{
    Info<< "tokens" << nl
        << "====" << nl;
    toString(Info, input)
        << nl
        << "====" << endl;

    ITstream is(input);
    Info<< is.size() << " tokens" << endl;

    if (verbose)
    {
        for (is.rewind(); !is.eof(); is.skip())
        {
            printToken(is.tokenIndex(), is.currentToken());
        }
        Info<< nl;
    }
    else
    {
        Info<< "==>";
        for (const token& tok : is)
        {
            Info<< ' ' << tok;
        }
        Info<< nl;
    }

    Info<< nl
        << "removing sub-dictionary tokens" << nl;

    for (is.rewind(); !is.eof(); is.skip())
    {
        if (is.currentToken().isPunctuation(token::BEGIN_BLOCK))
        {
            labelRange slice(is.tokenIndex(), 0);

            #if 0
            // This is a bad way to remove things since we lose the parse
            // point!
            for (/*nil*/; !is.eof(); is.skip())
            {
                if (is.currentToken().isPunctuation(token::END_BLOCK))
                {
                    slice.size() = (is.tokenIndex() - slice.start()) + 1;
                    break;
                }
            }
            #else
            for (label toki = is.tokenIndex()+1; toki < is.size(); ++toki)
            {
                if (is.peekToken(toki).isPunctuation(token::END_BLOCK))
                {
                    slice.size() = (toki - slice.start()) + 1;
                    break;
                }
            }
            #endif

            Info<< "remove range: " << slice
                << " currentIndex: " << is.tokenIndex() << '/' << is.size()
                // NB peekToken handles out-of-range
                << " token: " << is.peekToken(is.tokenIndex()) << nl;

            const label nRemoved = is.remove(slice);

            Info<< "remove " << nRemoved
                << " new current: " << is.tokenIndex() << '/' << is.size()
                // NB peekToken handles out-of-range
                << " token: " << is.peekToken(is.tokenIndex()) << nl;

            Info<< "==>";
            for (const token& tok : is)
            {
                Info<< ' ' << tok;
            }
            Info<< nl << nl;
        }
    }
    Info<< nl;
}


void testRemoveDict(const std::string& input, const int verbose)
{
    Info<< "tokens" << nl
        << "====" << nl;
    toString(Info, input)
        << nl
        << "====" << endl;

    ITstream is(input);
    Info<< is.size() << " tokens" << endl;

    if (verbose)
    {
        for (is.rewind(); !is.eof(); is.skip())
        {
            printToken(is.tokenIndex(), is.currentToken());
        }
        Info<< nl;
    }
    else
    {
        Info<< "==>";
        for (const token& tok : is)
        {
            Info<< ' ' << tok;
        }
        Info<< nl;
    }

    for (label pos = 0; pos < is.size(); /*nil*/)
    {
        labelRange slice
        (
            is.find(token::BEGIN_BLOCK, token::END_BLOCK, pos)
        );

        if (slice.good())
        {
            pos = slice.end_value();

            tokenList::subList substream(is.slice(slice));

            Info<< "  dict " << slice << " ==>";
            for (const token& tok : substream)
            {
                Info<< ' ' << tok;
            }
            Info<< nl;
        }
        else
        {
            break;
        }
    }


    Info<< nl
        << "removing sub-dictionary tokens" << nl;

    for (is.rewind(); !is.eof(); is.skip())
    {
        if (is.currentToken().isPunctuation(token::BEGIN_BLOCK))
        {
            labelRange slice
            (
                is.find(token::BEGIN_BLOCK, token::END_BLOCK, is.tokenIndex())
            );

            if (slice.good())
            {
                ITstream substream(is.extract(slice));

                Info<< "got " << slice << " ==>";
                for (const token& tok : substream)
                {
                    Info<< ' ' << tok;
                }
                Info<< nl;

                dictionary dict(substream);

                Info<< "tokenIndex: " << is.tokenIndex() << nl;
                Info<< "sub-dict " << dict << nl;

                Info<< "remove range: " << slice
                    << " currentIndex: " << is.tokenIndex() << '/' << is.size()
                    << " token: " << is.peekToken(is.tokenIndex()) << nl;

                const label nRemoved = is.remove(slice);

                Info<< "remove " << nRemoved
                    << " new current: " << is.tokenIndex() << '/' << is.size()
                    << " token: " << is.peekToken(is.tokenIndex()) << nl;

                Info<< "==>";
                for (const token& tok : is)
                {
                    Info<< ' ' << tok;
                }
                Info<< nl << nl;

                // Reposition the parse point
                is.seek(slice.start());
                is.skip(-1);

                Info<< "continue after " << is.tokenIndex()
                    << " : " << is.peekToken(is.tokenIndex()) << nl;
            }
        }
    }
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addVerboseOption("additional verbosity");
    argList::addBoolOption("basic", "basic tests");
    argList::addBoolOption("rewrite", "test rewriting only");
    argList::addBoolOption("remove-dict", "test rewriting only");

    argList args(argc, argv);

    if
    (
        !args.found("basic")
     && !args.found("rewrite")
     && !args.found("remove-dict")
    )
    {
        Info<< "No test options specified!" << nl << nl;
    }

    if (args.found("basic"))
    {
        const char* charInput =
            "( const char input \"string\" to tokenize )\n"
            "List<label> 5(0 1 2 3 4);";

        string stringInput("( string     ;    input \"string\" to tokenize )");

        List<char> listInput
        (
            ListOps::create<char>
            (
                stringInput.cbegin(),
                stringInput.cend(),
                Foam::identityOp{}
            )
        );

        doTest("empty", "", true, true);

        doTest("char*", charInput, true, true);
        doTest("string", stringInput, true);
        doTest("List<char>", listInput, true);

        reverse(listInput);
        doTest("List<char>", listInput, true);


        {
            auto& empty = ITstream::empty_stream();
            Info<< "The empty stream:" << nl
                << "    name: " << empty.name()
                << " good:" << empty.good()
                << " content:" << empty << nl;

            const token& tok0 = empty.peek();
            Info<< "First token from empty:" << tok0.info() << nl;
        }
    }

    if (args.found("rewrite"))
    {
        testWalk1
        (
            "std::string",
            "( string ;    input \"string\" to tokenize )"
            "{ other entry; value 100; value2 200; }"
            , args.verbose()
        );

        testRewrite
        (
            "some entry ( string1 ; )"
            "{ sub dict1; value 100; value2 200; }"
            "other entry ( string2 ; )"
            "{ sub dict2; value 100; value2 200; }"
            "{ sub dict3; value 100; value2 200; }"
            "trailing entry"
            , args.verbose()
        );
    }

    if (args.found("remove-dict"))
    {
        testRemoveDict
        (
            "some entry ( string1 ; )"
            "{ sub dict1; value 100; value2 200; }"
            "other entry ( string2 ; )"
            "{ sub dict2; value 100; value2 200; }"
            "{ sub dict3; value 100; value2 200; }"
            "trailing entry"
            , args.verbose()
        );

        testRemoveDict
        (
            "some entry no dictionary"
            , args.verbose()
        );
        testRemoveDict
        (
            "{ leading dict; } last-stuff"
            , args.verbose()
        );
        testRemoveDict
        (
            "first-stuff { trailing dict; }"
            , args.verbose()
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
