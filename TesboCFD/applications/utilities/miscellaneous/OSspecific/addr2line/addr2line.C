/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 Alexey Matveichev
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
    addr2line

Description
    A simple, partial emulation of addr2line utility for Mac-OS.

\*---------------------------------------------------------------------------*/

#include <getopt.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

#include "regExp.H"
#include "SubStrings.H"

static void usage();
static void version();
static std::string getLine(const std::string&, const std::string&);
static std::string pipeOpen(const std::string& cmd, const int lineNum = 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    int optHelp = 0, optFunctions = 0, optVersion = 0;
    int ch;
    std::string filename = "a.out";
    std::vector<std::string> addresses;

    static struct option opts[] =
    {
        { "target", required_argument, nullptr, 'b' },
        { "demangle", required_argument, nullptr, 'C' },
        { "exe", required_argument, nullptr, 'e' },
        { "functions", no_argument, &optFunctions, 1 },
        { "version", no_argument, &optVersion, 1 },
        { "basename", no_argument, nullptr, 's' },
        { "inlines", no_argument, nullptr, 'i' },
        { "section", required_argument, nullptr, 'j' },
        { "help", no_argument, &optHelp, 1 },
        { nullptr, 0, nullptr, 0 }
    };

    while ((ch = getopt_long(argc, argv, "b:C:e:fVsij:H", opts, nullptr)) != -1)
    {
        switch (ch)
        {
            case 'e':
                filename = std::string(optarg);
                break;
            case 'C':
                // Ignoring this flag for now
                break;
            case 'f':
                // Functions are demangled in printStack
                break;
            case 0:
                if (optHelp) usage();
                if (optVersion) version();
                break;
            default:
                usage();
                break;
        }
    }

    if (optind >= argc)
    {
        usage();
    }

    argc -= optind;
    argv += optind;

    while (argc > 0)
    {
        addresses.push_back(std::string(*argv));
        ++argv;
        --argc;
    }

    for (const auto& addr : addresses)
    {
        std::cout<< '\n' << getLine(filename, addr).c_str() << '\n';
    }

    return 0;
}


void usage()
{
    std::cout
        << "usage: addr2line [-e filename|--exe=filename]"
           " address [address...]\n" << std::endl;
    std::exit(1);
}


void version()
{
    std::cout<< "OpenFOAM addr2line emulator\n" << std::endl;
    std::exit(0);
}


// Read up to and including lineNum from the piped command
// Return the final line read
std::string pipeOpen(const std::string& cmd, int lineNum)
{
    std::string str;

    FILE* handle = popen(cmd.c_str(), "r");
    if (!handle) return str;

    char* buf = nullptr;
    size_t len = 0;
    ssize_t nread;

    // Read lineNum number of lines
    for
    (
        int cnt = 0;
        cnt <= lineNum && (nread = ::getline(&buf, &len, handle)) >= 0;
        ++cnt
    )
    {
        if (cnt == lineNum)
        {
            // Retain the last line, trimming trailing newline
            str.assign(buf);

            if (str.size())
            {
                str.resize(str.size()-1);
            }
        }
    }

    free(buf);
    pclose(handle);

    return str;
}


std::string getLine(const std::string& filename, const std::string& addr)
{
    std::string line =
        pipeOpen
        (
            "echo 'image lookup -va " + addr
          + "'"
          + " | xcrun lldb "
          + "-O 'target create --no-dependents -a x86_64 "
          + filename
          + "' -o '"
          + "target module load -f "
          + filename
          + " __TEXT 0x0' 2>/dev/null"
          + " | grep LineEntry"
        );


    Foam::regExp re(".+LineEntry: .+: (.+):([0-9]+):[0-9]+");

    Foam::regExp::results_type groups;
    if (!re.match(line, groups))
    {
        line = "??:0";
    }
    else
    {
        line = groups[1].str() + ":" + groups[2].str();
    }

    return line;
}


// ************************************************************************* //
