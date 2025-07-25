/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Test OFstream. Primarily atomic operations

\*---------------------------------------------------------------------------*/

#include "Fstream.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "argList.H"
#include "clock.H"
#include "Switch.H"
#include "ListOps.H"

using namespace Foam;

std::string time_stamp;

void listFiles(const fileName& dir)
{
    wordList files = ListOps::create<word>
    (
        readDir(dir, fileName::FILE),
        nameOp<fileName>()
    );

    Info
        << nl
        << "files:" << nl
        << files << nl
        << "ls" << nl
        << "============" << endl;
    Foam::system("ls -al " + dir);
    Info<< "============" << endl;
}


OSstream& printInfo(OFstream& os)
{
    InfoErr
        << "open: " << os.name() << nl
        << "appending: " << Switch::name(os.is_appending())
        << " tellp: "<< os.stdStream().tellp()
        << " gz: " << Switch::name(os.compression()) << nl;

    return InfoErr.stream();
}


void withHeader(OFstream& os)
{
    const auto tellp = os.stdStream().tellp();

    if (tellp <= 0)
    {
        InfoErr
            << "Add header" << nl;
        os  << "HEADER: " << time_stamp.c_str() << nl;
    }
}


template<class OSstreamType>
void generateLines(OSstreamType& os, label count = 1)
{
    for (label line = 1; line <= count; ++line)
    {
        os  << "[" << line
            << "] =============================================" << nl;
    }
}


template<class OSstreamType>
void generateContent
(
    OSstreamType& os,
    const bool with_seekend,
    const bool test_overwrite = false,
    const int64_t seek_out = -1
)
{
    if (with_seekend)
    {
        os.stdStream().seekp(0, std::ios_base::end);
        // OR? os.seek_end();
    }

    printInfo(os);

    withHeader(os);

    if (test_overwrite && seek_out >= 0)
    {
        InfoErr<< "... seekp(" << seek_out << ")" << nl;

        auto& oss = os.stdStream();

        // Actually std::streampos, but cannot increment that

        int64_t pos(seek_out);

        const int64_t tellp_end = oss.tellp();

        if (pos >= 0 && pos < tellp_end)
        {
            InfoErr
                << "... fill from " << label(pos)
                << " to " << label(tellp_end) << nl;

            oss.seekp(pos);

            while (pos < tellp_end)
            {
                // Fill with char 'X', rely on streambuf buffering
                oss << 'X';
                ++pos;
            }

            oss.seekp(seek_out);
            os << "More content [at " << seek_out << ']' << endl;
        }
    }

    generateLines(os, 4);

    printInfo(os)
        << "... sleep" << endl;

    listFiles(os.name().path());

    sleep(2);

    os << "[new content] +++++++++++++++++++++++++++++++++++" << endl;
}


template<class OSstreamType>
void generateOverwriteContent
(
    OSstreamType& os,
    const bool with_seekend,
    const int64_t seek_out = -1
)
{
    generateContent(os, with_seekend, true, seek_out);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption("gz", "Use compression");
    argList::addBoolOption("append-app", "Use append app mode");
    argList::addBoolOption("append-ate", "Use append ate mode");
    argList::addBoolOption("seekend", "Seek to end after non-append open");
    argList::addOption("seek", "value", "Seek from start (default: 100)");
    argList::addBoolOption("atomic", "Use atomic");
    argList::addBoolOption("keep", "Do not remove test directory");
    argList::addOption("write", "file", "test writing to file");

    #include "setRootCase.H"

    // Same time-stamp for all generated files
    time_stamp = clock::dateTime();

    const fileName baseDir("Test-OFstream-directory");

    Foam::mkDir(baseDir);

    InfoErr<< "mkdir: " << baseDir << endl;

    Info<< "start:" << nl;
    listFiles(baseDir);

    const bool with_seekend = args.found("seekend");

    const int seek_out = args.getOrDefault<int>("seek", 100);

    IOstreamOption streamOpt;

    if (args.found("gz"))
    {
        streamOpt.compression(IOstreamOption::COMPRESSED);
    }

    IOstreamOption::appendType append =
    (
        args.found("append-app") ? IOstreamOption::APPEND_APP
      : args.found("append-ate") ? IOstreamOption::APPEND_ATE
      : IOstreamOption::NO_APPEND
    );

    IOstreamOption::atomicType atomic =
    (
        args.found("atomic")
      ? IOstreamOption::ATOMIC
      : IOstreamOption::NON_ATOMIC
    );

    {
        OFstream(baseDir/"dummy")() << "Some file content" << endl;

        Foam::ln("dummy", baseDir/"Test3.txt");
        Foam::ln("dummy", baseDir/"Test4.txt");
        Foam::ln("dummy", baseDir/"Test4.txt.gz");
        Foam::ln("dummy", baseDir/"Test5.txt");
        Foam::ln("dummy", baseDir/"Test5.txt.gz");
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test1.txt",
            streamOpt,
            append
        );

        generateOverwriteContent(os, with_seekend, seek_out);
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test1-app.txt",
            streamOpt,
            IOstreamOption::APPEND_APP
        );

        generateOverwriteContent(os, with_seekend, seek_out);
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test1-ate.txt",
            streamOpt,
            IOstreamOption::APPEND_ATE
        );

        generateOverwriteContent(os, with_seekend, seek_out);
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test2.txt",
            streamOpt
        );

        generateContent(os, with_seekend);
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test3.txt",
            streamOpt,
            IOstreamOption::APPEND_APP
        );

        generateContent(os, with_seekend, with_seekend);
    }
    {
        OFstream os
        (
            baseDir/"Test4.txt",
            IOstreamOption::ASCII,
            IOstreamOption::COMPRESSED
        );

        // No seekend with COMPRESSED
        generateContent(os, false);
    }
    {
        OFstream os
        (
            IOstreamOption::ATOMIC,
            baseDir/"Test5.txt"
        );

        generateContent(os, with_seekend);
    }

    Info<< nl << "done:" << endl;

    listFiles(baseDir);

    if (args.found("keep"))
    {
        InfoErr<< "keep: " << baseDir << endl;
    }
    else
    {
        InfoErr<< "rmdir: " << baseDir << endl;
        Foam::rmDir(baseDir);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
