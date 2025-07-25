/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Test some string functionality

\*---------------------------------------------------------------------------*/

#include "fieldExprFwd.H"
#include "int.H"
#include "uint.H"
#include "scalar.H"
#include "Switch.H"
#include "fileName.H"
#include "stringList.H"
#include "stringOps.H"
#include "dictionary.H"
#include "IOstreams.H"
#include "OSspecific.H"


using namespace Foam;

void testCommentStripping(const std::string& s)
{
    Info<< "input" << nl
        << "========" << nl
        << s << nl
        << "========" << nl;

    Info<< "output" << nl
        << "========" << nl
        << stringOps::removeComments(s) << nl
        << "========" << nl << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    string test
    (
        "  $HOME kjhkjhkjh \" \\$HOME/tyetyery $; ${FOAM_RUN} \n $; hkjh;"
        " $(DONOTSUBST) some other <${USER}> with '${__UNKNOWN:-some default}'"
        " value "
        " or with '${HOME:+Home was set}' via :+ alternative"
        " or with '${__UNKNOWN:+unknown}' empty  "
    );

    setEnv("FOAM_CASE", cwd(), true);

    dictionary dict;
    dict.add("HOME", "myHome");

    dictionary subDict;
    subDict.add("value1", "test1");
    subDict.add("value2", "test2");
    dict.add("FOAM_RUN", subDict);

    if (false)
    {
        typedef std::string inputType;
        typedef string outputType;

        inputType in1("move-construct-from");

        Info<<"move construct from " << in1.length() << nl;

        outputType out1(std::move(in1));

        Info<<"after move construct "
            << out1.size() << ", " << in1.size() << nl;

        in1 = "move-assign-from";
        out1 = "some-text-rubbish";
        out1.resize(10);

        Info<<"move assign from " << in1.length() << nl;

        out1 = std::move(in1);

        Info<<"after move assign "
            << out1.size() << ", " << in1.size() << nl;

        return 0;
    }


    // basic expansions
    {
        for
        (
            const auto& cstr
          :
            {
                "~OpenFOAM/controlDict",    "<etc>/controlDict",
                "<etc:ugo>/controlDict",
                "<etc:u>/controlDict",
                "<etc:ug>/controlDict",
                "<etc:go>/controlDict",
                "<etc:o>/controlDict",
                "<etc:JUNK>/controlDict",   // rubbish input
                "<etc:>/controlDict",       // rubbish input
                "$FOAM_CASE/xyz",           "<case>/xyz",
                "$FOAM_CASE/constant/xyz",  "<constant>/xyz",
                "$FOAM_CASE/system/xyz",    "<system>/xyz",

                // corner cases
                "~OpenFOAM",                "<etc>",
                "~OpenFOAM/",               "<etc>/",
                "$FOAM_CASE",               "<case>",
                "$FOAM_CASE/constant",      "<constant>",
                "$FOAM_CASE/system",        "<system>",

                "$FOAM_CASE/",              "<case>/",
                "$FOAM_CASE/constant/",     "<constant>/",
                "$FOAM_CASE/system/",       "<system>/",
            }
        )
        {
            string input(cstr);
            string output(stringOps::expand(input));

            Info<< "input:  " << input  << nl
                << "expand: " << output << nl
                << "split: " << stringOps::split(output, '/') << nl << nl;
        }
    }


    // Test Foam::name with formatting string
    {
        word formatted = word::printf("formatted=<%X>", 0xdeadbeef);
        Info<<"formatted: " << formatted << nl;
    }

    Info<<"formatted: "
        << word::printf("formatted not checked for validity=<%X>", 0xdeadbeef)
        << nl
        << endl;


    Info<< "string:" << test << nl << "hash:"
        << unsigned(string::hasher()(test)) << endl;

    Info<<"trimLeft: " << stringOps::trimLeft(test) << endl;
    Info<<"trimRight: " << stringOps::trimRight(test) << endl;
    Info<<"trim: " << stringOps::trim(test) << endl;

    // With replace, replaceAll
    {
        string test2(test);
        Info<<"replaceAny: (\"abcj\", '?')" << nl
            << test2.replaceAny("abcj", '?') << endl;
        Info<<"replaceAll: (\"k\", \"?\")" << nl
            << test2.replaceAll("k", "?") << endl;
        Info<<"replaceAll: (\"h\", null) - ie, remove them" << nl
            << test2.replaceAll("h", "") << endl;
        Info<<"replaceAny: (\"it${}?\", null) - ie, remove them" << nl
            << test2.replaceAny("it${}?", '\0') << endl;
    }

    if (false)
    {
        Info<<"test move construct - string size:" << test.size() << nl;
        string test2(std::move(test));

        Info<<"input size:" << test.size() << nl;
        Info<<"moved size:" << test2.size() << nl;

        Info<<"test move assign - string sizes:"
            << test.size() << "/" << test2.size() << nl;

        test = std::move(test2);

        Info<<"input size:" << test.size() << nl;
        Info<<"moved size:" << test2.size() << nl;
    }

    if (false)
    {
        std::string str("some text");

        Info<<"test move construct to string:" << str.size() << nl;

        Foam::string test2(std::move(str));

        Info<<"input/moved sizes:" << str.size() << "/" << test2.size() << nl;

        str = std::move(test2);

        Info<<"test move assign - sizes:"
            << str.size() << "/" << test2.size() << nl;
    }

    if (false)
    {
        Foam::string str("thisIsAWord");

        Info<<"test move construct to word:" << str.size() << nl;

        word test2(std::move(str));

        Info<<"input/moved sizes:" << str.size() << "/" << test2.size() << nl;

        str = std::move(test2);

        Info<<"test move assign - sizes:"
            << str.size() << "/" << test2.size() << nl;

        // move back
        test2.swap(str);

        Info<<"test move assign - sizes:"
            << str.size() << "/" << test2.size() << nl;

        string str2(std::move(test2));
        Info<<"input/moved sizes:" << test2.size() << "/" << str2.size() << nl;

    }

    {
        fileName test1("libFooBar.so");

        Info<< nl;
        Info<< "trim filename: " << test1 << nl;

        test1.removeStart("lib");
        Info<<"without leading 'lib': " << test1 << nl;
    }

    Info<< nl;
    Info<<"camel-case => " << (word("camel") & "case") << nl;
    for (const auto& s : { " text with \"spaces'", "08/15 value" })
    {
        // Character sequence

        Info<<"validated 5 chars from \" => "
            << word::validate(s, s+5, true) << nl;

        Info<<"validated (via string convert) \"" << s << "\" => "
            << word::validate(s, true) << nl;
    }
    Info<< nl;

    // test sub-strings via iterators
    string::const_iterator iter  = test.end();
    string::const_iterator iter2 = test.end();
    string::size_type fnd = test.find('\\');

    if (fnd != string::npos)
    {
        iter  = test.begin() + fnd;
        iter2 = iter + 6;
    }

    Info<< "sub-string via iterators : >";
    while (iter != iter2)
    {
        Info<< *iter;
        ++iter;
    }
    Info<< "<\n";

    Info<< string(test).replace("kj", "zzz") << endl;
    Info<< string(test).replace("kj", "") << endl;
    Info<< string(test).replaceAll("kj", "zzz") << endl;
    Info<< string(test).replaceAll("kj", "z") << endl;

    Info<< "expanded: " << string(test).expand() << endl;

    Info<<"dictionary-based substitution: " << dict << endl;
    Info<< "expand dict: " << stringOps::expand(test, dict) << endl;

    string test2("~OpenFOAM/controlDict");
    Info<< test2 << " => " << test2.expand() << endl;

    // replace controlDict with "newName"
    {
        string::size_type i = test2.rfind('/');

        if (i == string::npos)
        {
            test2 = "newName";
        }
        else
        {
            // this uses the std::string::replace
            test2.replace(i+1, string::npos, "newName");
        }
        Info<< "after replace: " << test2 << endl;

        // do another replace
        // this uses the Foam::string::replace
        test2.replace("OpenFOAM", "openfoam");

        Info<< "after replace: " << test2 << endl;
    }

    cout<< "\nEnter some string to test:\n";

    string s;
    Sin.getLine(s);

    string s2(s.expand());

    cout<< "output string with " << s2.length() << " characters\n";
    cout<< "ostream<<  >" << s2 << "<\n";
    Info<< "Ostream<<  >" << s2 << "<\n";
    Info<< "hash:" << hex << string::hasher()(s2) << dec << endl;

    cout<< "\ntest Foam::name()\n";

    Info<< "hash: = " << word::printf("0x%012X", string::hasher()(s2)) << endl;

    // Test formatting on int
    {
        label val = 25;
        Info<< "int " << val << " nameOp='"
            << nameOp<label>()(val) << "', name='"
            << Foam::name(val) << "' or "
            << word::printf("formatted '%08d'", val) << "\n";
    }

    // Test formatting on scalar
    {
        scalar val = 3.1415926535897931;
        Info<< "scalar " << val << " nameOp='"
            << nameOp<scalar>()(val) << "', name='"
            << Foam::name(val) << "' or "
            << word::printf("formatted '%.9f'", val) << "\n";
    }

    // Test formatting on uint
    {
        uint64_t val = 25000000ul;
        Info<< "uint64 " << val << " nameOp='"
            << nameOp<uint64_t>()(val) << "', name='"
            << Foam::name(val) << "' or "
            << word::printf("formatted '%08d'", val) << "\n";
    }

    // test starts_with, ends_with methods
    {
        string empty; //;
        string input1 = "shorter input";
        string input2 = "longer input text";

        stringList checks{"match", "long", "", "short", "text", "s", "l", "t"};

        Info<< nl;
        Info<< "check starts_with:" << nl
            << "~~~~~~~~~~~~~~~~~" << nl;

        Info<<"input: " << empty << nl;
        for (const string& test : checks)
        {
            Info<< "    starts_with(" << test << ") = "
                << Switch(empty.starts_with(test)) << nl;
        }
        Info<<"input: " << input1 << nl;
        for (const string& test : checks)
        {
            Info<< "    starts_with(" << test << ") = "
                << Switch(input1.starts_with(test)) << nl;
        }
        Info<<"input: " << input2 << nl;
        for (const string& test : checks)
        {
            Info<< "    starts_with(" << test << ") = "
                << Switch(input2.starts_with(test)) << nl;
        }


        Info<< nl;
        Info<< "check ends_with:" << nl
            << "~~~~~~~~~~~~~~~~~" << nl;

        Info<<"input: " << empty << nl;
        for (const string& test : checks)
        {
            Info<< "    ends_with(" << test << ") = "
                << Switch(empty.ends_with(test)) << nl;
        }
        Info<<"input: " << input1 << nl;
        for (const string& test : checks)
        {
            Info<< "    ends_with(" << test << ") = "
                << Switch(input1.ends_with(test)) << nl;
        }
        Info<<"input: " << input2 << nl;
        for (const string& test : checks)
        {
            Info<< "    ends_with(" << test << ") = "
                << Switch(input2.ends_with(test)) << nl;
        }

        Info<< nl;
        Info<< "check ends_with as applied to field names:" << nl
            << "~~~~~~~~~~~~~~~~~" << nl;

        string input3 = "field_0";
        string input4 = "_0";

        Info<<input3 << " ends_with(\"_0\") = "
            << Switch(input3.ends_with("_0")) << nl;

        Info<<input4 << " ends_with(\"_0\") = "
            << Switch(input4.ends_with("_0")) << nl;
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
