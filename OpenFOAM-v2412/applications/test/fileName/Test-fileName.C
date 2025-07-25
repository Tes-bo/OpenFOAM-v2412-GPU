/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Application
    Test-fileName

Description
    Test some basic fileName functionality

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "SubList.H"
#include "DynamicList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "POSIX.H"
#include "Switch.H"
#include "etcFiles.H"
#include "Pair.H"
#include "Tuple2.H"
#include <fstream>

using namespace Foam;

// Create named file with some dummy content
void touchFileContent(const fileName& file)
{
    std::ofstream os(file);
    os << "file=<" << file << ">" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

unsigned testClean(std::initializer_list<Pair<std::string>> tests)
{
    unsigned nFail = 0;

    for (const Pair<std::string>& test : tests)
    {
        const std::string& input = test.first();
        const std::string& expected = test.second();

        fileName cleaned(test.first());
        cleaned.clean();  // Remove unneeded ".."

        if (cleaned == expected)
        {
            Info<< "(pass)"
                << "  clean " << input << " -> " << cleaned << nl;
        }
        else
        {
            Info<< "(fail)"
                << "  clean " << input << " -> " << cleaned
                << "  expected=" << expected
                << nl;

            ++nFail;
        }
    }

    return nFail;
}


unsigned testStrip
(
    const bool doClean,
    std::initializer_list
    <
        Tuple2<bool, std::string>
    > tests
)
{
    Info<< nl << "Checking with clean=" << Switch(doClean) << nl << endl;

    unsigned nFail = 0;

    for (const Tuple2<bool, std::string>& test : tests)
    {
        const bool expected = test.first();
        const std::string& input = test.second();

        fileName output(fileName::validate(input, doClean));

        // Check for real failure (invalid chars) vs.
        // spurious failure (removed double slashes with 'doClean').

        const bool same =
        (
            doClean
          ? fileName::equals(input, output)
          : (input == output)
        );

        if (same)
        {
            if (expected)
            {
                Info<< "(pass) validated " << input << " = " << output << nl;
            }
            else
            {
                ++nFail;
                Info<< "(fail) unexpected success for " << input << nl;
            }
        }
        else
        {
            if (expected)
            {
                ++nFail;
                Info<< "(fail) unexpected";
            }
            else
            {
                Info<< "(pass) expected";
            }

            Info<< " failure for " << input << nl;
        }
    }

    return nFail;
}


unsigned testEquals
(
    std::initializer_list
    <
        Tuple2<bool, Pair<std::string>>
    > tests
)
{
    Info<< nl << "Checking fileName::equals()" << nl << endl;

    unsigned nFail = 0;

    for (const Tuple2<bool, Pair<std::string>>& test : tests)
    {
        const bool expected = test.first();
        const std::string& s1 = test.second().first();
        const std::string& s2 = test.second().second();

        const bool same = fileName::equals(s1, s2);

        if (same)
        {
            if (expected)
            {
                Info<< "(pass) success";
            }
            else
            {
                ++nFail;
                Info<< "(fail) unexpected success";
            }
        }
        else
        {
            if (expected)
            {
                ++nFail;
                Info<< "(fail) unexpected failure";
            }
            else
            {
                Info<< "(pass) expected failure";
            }

        }

        Info<< " for " << s1 << " == " << s2 << nl;
    }
    return nFail;
}


unsigned testRelative(std::initializer_list<Pair<std::string>> tests)
{
    Info<< nl << "Checking fileName::relative()" << nl << endl;

    unsigned nFail = 0;

    for (const Pair<std::string>& test : tests)
    {
        const std::string& dir    = test.first();
        const std::string& parent = test.second();

        Info<< "directory: " << dir << nl
            << "parent   : " << parent << nl
            << "relative = " << fileName(dir).relative(parent) << nl
            << "case-rel = " << fileName(dir).relative(parent, true) << nl
            << endl;
    }

    return nFail;
}


void testDirname(const fileName& input)
{
    Info<< "input:" << input
        << "   path:" << input.path()
        << "   name:\"" << input.name() << '"'
        << "   ext:\"" << input.ext()  << '"'
        << "   components: " << flatOutput(input.components())
        << "   last: " << input.component(string::npos) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("validate", "test fileName::validate");
    argList::addBoolOption("ext", "test handing of file extensions");
    argList::addBoolOption("construct", "test constructors");
    argList::addBoolOption("relative", "test relative operations");
    argList::addBoolOption("system", "test filesystem operations");
    argList::addBoolOption("default", "reinstate default tests");
    argList::addBoolOption("clean", "clean()");
    argList::addBoolOption("dirname", "basename/dirname tables");
    argList::addNote("runs default tests or specified ones only");

    #include "setRootCase.H"

    // Run default tests, unless specific tests are requested
    const bool defaultTests =
        args.found("default") || args.options().empty();

    if (args.found("construct"))
    {
        Info<< "From initializer_list<word> = ";
        fileName file1
        {
            "hello",
            "hello1",
            "hello2",
            "hello3",
            "hello4.ext"
        };

        Info<< file1 << nl;

        Info<< "From initializer_list<fileName> = ";
        fileName file2
        {
            file1,
            "some",
            "more/things.ext"
        };

        Info<< file2 << nl;


        Info<< "From initializer_list<fileName> with nesting = ";
        fileName file3
        {
            std::string("ffO"),
            "some",
            "more/things.ext"
        };
        Info<< file3 << nl;

        DynamicList<word> base
        {
            "hello",
            "hello1"
        };

        fileName file4
        {
            "some",
            file3,
            "more/things.ext",
            file1
        };
        Info<< "All ==> " << file4 << nl;
    }

    if (args.found("dirname"))
    {
        testDirname("");
        testDirname(".");
        testDirname("abc");
        testDirname("/");
        testDirname("/abc");
        testDirname("abc/def");
        testDirname("/abc/def");
        testDirname("/abc/def/");
        testDirname("/abc///def///");
        testDirname("/abc/../def");
    }


    // Test various ext() methods
    if (args.found("ext"))
    {
        Info<<nl << nl << "handling of fileName extension" << nl;

        fileName empty;
        fileName endWithDot("some.path/name.");
        fileName endWithSlash("some.path/");
        fileName input0("some.file/with.out/extension");
        fileName input1("path.to/media/image.png");

        Info<<"File : " << input0 << " ext: "
            << Switch(input0.has_ext())
            <<  " = " << input0.ext() << nl;
        Info<<"File : " << input1 << " ext: "
            << Switch(input1.has_ext())
            <<  " = " << input1.ext() << nl;
        Info<<"File : " << endWithDot << " ext: "
            << Switch(endWithDot.has_ext())
            <<  " = " << endWithDot.ext() << " <-- perhaps return false?" << nl;
        Info<<"File : " << endWithSlash << " ext: "
            << Switch(endWithSlash.has_ext())
            <<  " = " << endWithSlash.ext() << nl;


        Info<<"Remove extension " << (input0.remove_ext());
        Info<< "  now: " << input0 << nl;

        Info<<"Remove extension " << (input1.removeExt());
        Info<< "  now: " << input1 << nl;

        Info<<"Remove extension " << (endWithSlash.remove_ext());
        Info<< "  now: " << endWithSlash << nl;

        wordList exts{ "jpg", "png", "txt", word::null };
        Info<<"Add extension(s): " << input1 << nl;
        for (const word& e : exts)
        {
            Info<<"<" << e << ">  -> " << input1.ext(e) << nl;
        }
        Info<< nl;


        Info<<"Test has_ext(word)" << nl
            <<"~~~~~~~~~~~~~~~~~" << nl;
        Info<<"Has extension(s):" << nl
            << "input: " << input1 << nl;
        for (const word& e : exts)
        {
            Info<<"  '" << e << "'  -> "
                << Switch(input1.has_ext(e)) << nl;
        }
        Info<< nl;

        Info<<"Has extension(s):" << nl
            << "input: " << endWithDot << nl;
        for (const word& e : exts)
        {
            Info<<"  '" << e << "'  -> "
                << Switch(endWithDot.has_ext(e)) << nl;
        }
        Info<< nl;


        Info<<"Test has_ext(wordRe)" << nl
            <<"~~~~~~~~~~~~~~~~~~~" << nl;

        // A regex with a zero length matcher doesn't work at all:
        // eg "(png|jpg|txt|)" regex matcher itself

        wordRe matcher0("()", wordRe::REGEX);
        wordRe matcher1("(png|jpg|txt)", wordRe::REGEX);
        wordRe matcher2("(png|txt)", wordRe::REGEX);

        Info<<"Has extension(s):" << nl
            << "input: " << endWithDot << nl;
        Info<<"    " << matcher0 << "  -> "
            << Switch(endWithDot.has_ext(matcher0)) << nl;
        Info<<"    " << matcher1 << "  -> "
            << Switch(endWithDot.has_ext(matcher1)) << nl;
        Info<<"    " << matcher2 << "  -> "
            << Switch(endWithDot.has_ext(matcher2)) << nl;

        Info<< "input: " << input1 << nl;
        Info<<"    " << matcher0 << "  -> "
            << Switch(input1.has_ext(matcher0)) << nl;
        Info<<"    " << matcher1 << "  -> "
            << Switch(input1.has_ext(matcher1)) << nl;
        Info<<"    " << matcher2 << "  -> "
            << Switch(input1.has_ext(matcher2)) << nl;
        Info<< nl;

        Info<<"Remove extension(s):" << nl << "input: " << input1 << nl;
        while (!input1.empty())
        {
            if (input1.remove_ext())
            {
                Info<< "   -> " << input1 << nl;
            }
            else
            {
                Info<< "stop> " << input1 << nl;
                break;
            }
        }
        Info<< nl;

        input0.clear();
        Info<<"test with zero-sized: " << input0 << nl;
        Info<<"add extension: " << input0.ext("abc") << nl;
        Info<< nl;

        input0 = "this/";
        Info<<"test add after slash: " << input0 << nl;
        Info<<"add extension: " << input0.ext("abc")
            << " <-- avoids accidentally creating hidden files" << nl;
        Info<< nl;

        input0 = "this.file.";
        Info<<"test after dot: " << input0 << nl;
        Info<<"add extension: " << input0.ext("abc")
            << " <-- No check for repeated dots (user error!)" << nl;
        Info<< nl;
    }


    if (args.found("clean"))
    {
        Info<< nl << "Test fileName::clean()" << nl << nl;

        unsigned nFail = testClean
        ({
            { "/", "/" },
            { "/abc/", "/abc" },
            { "/abc////def", "/abc/def" },
            { "/abc/def/./ghi/.", "/abc/def/ghi" },
            { "abc/def/./",   "abc/def" },
            { "./abc/", "./abc" },
            { "/abc/def/../ghi/jkl/nmo/..", "/abc/ghi/jkl" },
            { "abc/../def/ghi/../jkl", "abc/../def/jkl" },
        });

        Info<< nl;
        if (nFail)
        {
            Info<< "failed " << nFail;
        }
        else
        {
            Info<< "passed all";
        }
        Info<< " fileName::clean tests" << nl;
    }


    if (args.found("validate"))
    {
        unsigned nFail = 0;
        Info<< nl << "Test fileName::validate" << nl;

        // Without clean
        nFail += testEquals
        (
            {
                { true,  { "abc", "abc/" } },
                { true,  { "///abc/", "//abc///" } },
                { false, { " ab //c/", "ab/c" } },
            }
        );

        Info<< nl << "Test fileName::validate" << nl;

        // Without clean
        nFail += testStrip
        (
            false,
            {
                { true,  "abc/" },
                { true,  "/", },
                { true,  "//", },
                { true,  "/abc/def", },
                { true,  "/abc/def/", },
                { false, "/abc  def" },
                { true,  "/abc////def///", },
                { false,  "/abc//// def///" },
            }
        );

        // With clean
        nFail += testStrip
        (
            true,
            {
                { true,  "abc/" },
                { true,  "/" },
                { true,  "//" },
                { true,  "/abc/def" },
                { true,  "/abc/def/" },
                { false, "/abc  def" },
                { true,  "/abc////def///" },
                { false, "/abc//// def///" },
            }
        );

        Info<< nl;
        if (nFail)
        {
            Info<< "failed " << nFail;
        }
        else
        {
            Info<< "passed all";
        }
        Info<< " fileName::validate tests" << nl;
    }


    if (args.found("relative"))
    {
        unsigned nFail = 0;

        nFail += testRelative
        (
            {
                { "",  "" },
                { "",  "/" },
                { "",  "/some" },

                { "/some/dir/subdir/name",  "" },
                { "/some/dir/subdir/name",  "/" },
                { "/some/dir/subdir/name",  "/some" },
                { "/some/dir/subdir/name",  "/some/" },
                { "/some/dir/subdir/name",  "/some/dir" },
                { "/some/dir/subdir/name",  "/some/dir/" },
                { "/some/dir/subdir/name",  "/some/dir/subdir" },
                { "/some/dir/subdir/name",  "/some/dir/subdir/name" },
                { "/some/dir/subdir/name",  "/some/dir/subdir/name/" },
                { "/some/dir/subdir/name",  "/some/other" },

                // With single-char for name
                { "/some/dir/subdir/a",  "" },
                { "/some/dir/subdir/a",  "/" },
                { "/some/dir/subdir/a",  "/some" },
                { "/some/dir/subdir/a",  "/some/" },
                { "/some/dir/subdir/a",  "/some/dir" },
                { "/some/dir/subdir/a",  "/some/dir/" },
                { "/some/dir/subdir/a",  "/some/dir/subdir" },
                { "/some/dir/subdir/a",  "/some/dir/subdir/a" },
                { "/some/dir/subdir/a",  "/some/dir/subdir/a/" },
                { "/some/dir/subdir/a",  "/some/other" },

                // Bad input (trailing slash)
                { "/some/dir/subdir/a/",  "" },
                { "/some/dir/subdir/a/",  "/" },
                { "/some/dir/subdir/a/",  "/some" },
                { "/some/dir/subdir/a/",  "/some/" },
                { "/some/dir/subdir/a/",  "/some/dir" },
                { "/some/dir/subdir/a/",  "/some/dir/" },
                { "/some/dir/subdir/a/",  "/some/dir/subdir/" },
                { "/some/dir/subdir/a/",  "/some/dir/subdir/a" },
                { "/some/dir/subdir/a/",  "/some/dir/subdir/a/" },
                { "/some/dir/subdir/a/",  "/some/other" },
            }
        );
    }


    // Test some copying and deletion
    if (args.found("system"))
    {
        const fileName dirA("dirA");
        const fileName dirB("dirB");
        const fileName dirC("dirC");
        const fileName dirD("dirD");
        const fileName lnA("lnA");
        const fileName lnB("lnB");

        // Purge anything existing
        Foam::rmDir(dirA, true);
        Foam::rmDir(dirB, true);
        Foam::rmDir(dirC, true);
        Foam::rmDir(dirD, true);
        Foam::rm(lnA);
        Foam::rm(lnB);

        {
            fileName name(dirA/dirB/dirC/"abc");
            Foam::mkDir(name.path());
            touchFileContent(name);           // Create real file

            Foam::mkDir(dirB/dirB/dirB/dirB);
            Foam::ln("test", dirB/"linkB");   // Create dead link

            Foam::mkDir(dirC);
            Foam::ln("../dirD", dirC/"linkC");  // Create real link

            Foam::mkDir(dirD);

            for (const fileName& d : { dirA, dirB, dirC, dirD })
            {
                Info<< "Directory: " << d << " = "
                    << readDir(d, fileName::UNDEFINED, false, false) << nl;

                if (Foam::rmDir(d, false, true))
                {
                    Info<< "   Removed empty dir" << nl;
                }
                else
                {
                    Info<< "   Could not remove empty dir" << nl;
                }
            }

            // Force removal before continuing
            Foam::rmDir(dirA, true);
            Foam::rmDir(dirB, true);
            Foam::rmDir(dirC, true);
            Foam::rmDir(dirD, true);
        }

        Info<< nl << "=========================" << nl
            << "Test some copying and deletion" << endl;


        Info<< "Creating directory " << dirA << endl;
        Foam::mkDir(dirA);

        Info<< "Populating with various files" << endl;
        for
        (
            const std::string name
          : { "file-1", "file-2", "bad name one", "bad name 2" }
        )
        {
            // Full path, but without any stripping
            const fileName file
            (
                (static_cast<const std::string&>(dirA) + "/" + name),
                false
            );

            Info<<"    create: " << file << endl;
            touchFileContent(file);
        }

        const int oldDebug = POSIX::debug;
        POSIX::debug = 1;


        // Create link and test it
        Info<< "Creating softlink " << lnA << endl;
        Foam::ln(dirA, lnA);

        fileName::Type lnAType = lnA.type(false);

        if (lnAType != fileName::SYMLINK)
        {
            FatalErrorIn("Test-fileName") << "Type of softlink " << lnA
                << " should be " << fileName::SYMLINK
                << " but is " << lnAType << exit(FatalError);
        }

        fileName::Type dirAType = lnA.type(true);

        if (dirAType != fileName::DIRECTORY)
        {
            FatalErrorIn("Test-fileName") << "Type of what softlink " << lnA
                << " points to should be " << fileName::DIRECTORY
                << " but is " << dirAType << exit(FatalError);
        }

        // Copy link only
        {
            Info<< "Copying (non-follow) softlink " << lnA << " to " << lnB
                << endl;

            Foam::cp(lnA, lnB, false);
            if (lnB.type(false) != fileName::SYMLINK)
            {
                FatalErrorIn("Test-fileName") << "Type of softlink " << lnB
                    << " should be " << fileName::SYMLINK
                    << " but is " << lnB.type(false) << exit(FatalError);
            }
            if (lnB.type(true) != fileName::DIRECTORY)
            {
                FatalErrorIn("Test-fileName") << "Type of softlink " << lnB
                    << " should be " << fileName::DIRECTORY
                    << " but is " << lnB.type(true) << exit(FatalError);
            }

            // Delete (link)
            Foam::rm(lnB);
        }

        // Copy contents of link
        {
            Info<< "Copying (contents of) softlink " << lnA << " to " << lnB
                << endl;

            Foam::cp(lnA, lnB, true);
            if (lnB.type(false) != fileName::DIRECTORY)
            {
                FatalErrorIn("Test-fileName") << "Type of softlink " << lnB
                    << " should be " << fileName::DIRECTORY
                    << " but is " << lnB.type(false) << exit(FatalError);
            }

            // Delete (directory, not link)
            Foam::rmDir(lnB);
        }

        POSIX::debug = oldDebug;

        // Verify that rmDir works with bad names too
        Foam::rmDir(dirA);
        Foam::rm(lnA);
    }


    if (!defaultTests)
    {
        return 0;
    }

    DynamicList<word> wrdList
    {
        "hello",
        "hello1",
        "hello2",
        "hello3",
        "hello4.ext"
    };

    fileName pathName(wrdList);

    Info<< "pathName = " << pathName << nl
        << "nameOp   = " << nameOp<fileName>()(pathName) << nl
        << "pathName.name()     = >" << pathName.name() << "<\n"
        << "pathName.path()     = "  << pathName.path() << nl
        << "pathName.ext()      = >" << pathName.ext() << "<\n"
        << "pathName.stem       = >" << pathName.stem() << "<\n";

    Info<< "pathName.components() = " << pathName.components() << nl
        << "pathName.component(2) = " << pathName.component(2) << nl
        << endl;

    pathName.replace_name("newName.ext");
    Info<< "new name = " << pathName << nl;
    Info<< "has ext      = " << Switch::name(pathName.has_ext()) << nl;
    Info<< "has ext('')  = " << Switch::name(pathName.has_ext("")) << nl;
    Info<< "has ext(foo) = " << Switch::name(pathName.has_ext("foo")) << nl;
    Info<< "has ext(ext) = " << Switch::name(pathName.has_ext("ext")) << nl;

    pathName.replace_ext("png");
    Info<< "new ext = " << pathName << nl;

    pathName.replace_ext("");  // Same as remove_ext
    Info<< "new ext = " << pathName << nl;

    Info<< "has path = " << Switch::name(pathName.has_path()) << nl;
    pathName.remove_path();
    pathName.removePath();  // second type should be a no-op
    Info<< "removed path = " << pathName << nl;

    Info<< nl << nl;

    // try with different combination
    // The final one should emit warnings
    for (label start = 0; start <= wrdList.size(); ++start)
    {
        fileName instance, local;
        word name;

        fileName path(SubList<word>(wrdList, wrdList.size()-start, start));
        fileName path2 = "."/path;

        IOobject::fileNameComponents
        (
            path,
            instance,
            local,
            name
        );

        Info<< "IOobject::fileNameComponents for " << path << nl
            << "  instance = " << instance << nl
            << "  local    = " << local << nl
            << "  name     = " << name << endl;

        IOobject::fileNameComponents
        (
            path2,
            instance,
            local,
            name
        );

        Info<< "IOobject::fileNameComponents for " << path2 << nl
            << "  instance = " << instance << nl
            << "  local    = " << local << nl
            << "  name     = " << name << endl;

    }


    // test findEtcFile
    Info<< "\n\nfindEtcFile tests:" << nl
        << " controlDict => " << findEtcFile("controlDict") << nl
        << " badName => " << findEtcFile("badName") << endl;

    {

        Info<< nl << "Expect a FatalError for findEtcFile() with a bad name:"
            << nl;

        const bool oldThrowingError = FatalError.throwing(true);

        try
        {
            Info<< " badName(die) => " << flush
                << findEtcFile("<very-badName>", true) << nl
                << endl;
        }
        catch (const Foam::error& err)
        {
            Info<< nl << "findEtcFile() Caught FatalError "
                << err << nl << endl;
        }
        FatalError.throwing(oldThrowingError);
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
