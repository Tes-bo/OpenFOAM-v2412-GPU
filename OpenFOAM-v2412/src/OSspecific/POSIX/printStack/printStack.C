/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "error.H"
#include "OSspecific.H"

#include <cinttypes>
#include <sstream>
#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Read up to and including lineNum from the piped command
// Return the final line read
std::string pipeOpen(const std::string& cmd, const int lineNum = 0)
{
    std::string str;

    FILE *handle = popen(cmd.c_str(), "r");
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


inline std::string addressToWord(const uintptr_t addr)
{
    std::ostringstream buf;
    buf.setf(std::ios_base::hex, std::ios_base::basefield);

    buf << "0x";  // Same as setf(std::ios::showbase)

    #ifdef __APPLE__
    buf << uint64_t(addr);
    #else
    buf << addr;
    #endif

    return buf.str();  // Needs no stripping
}


// Note: demangle requires symbols only - without extra '(' etc.
inline std::string demangleSymbol(const char* sn)
{
    int st = 0;

    char* cxx_sname = abi::__cxa_demangle(sn, nullptr, nullptr, &st);

    if (st == 0 && cxx_sname)
    {
        std::string demangled(cxx_sname);
        free(cxx_sname);

        return demangled;
    }

    return sn;
}


inline Foam::string& shorterPath(Foam::string& s)
{
    s.replace(Foam::cwd() + '/', "");
    s.replace(Foam::home(), "~");
    return s;
}


void printSourceFileAndLine
(
    Foam::Ostream& os,
    const Foam::fileName& filename,
    const Dl_info& info,
    void *addr
)
{
    uintptr_t address = uintptr_t(addr);
    std::string myAddress = addressToWord(address);

    // Can use relative addresses for executables and libraries with the
    // Darwin addr2line implementation.
    // On other systems (Linux), only use relative addresses for libraries.

    #ifndef __APPLE__
    if (filename.has_ext("so"))
    #endif
    {
        // Convert address into offset into dynamic library
        uintptr_t offset = uintptr_t(info.dli_fbase);
        intptr_t relativeAddress = address - offset;
        myAddress = addressToWord(relativeAddress);
    }

    if (filename[0] == '/')
    {
        Foam::string line = pipeOpen
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + myAddress,
            1
        );

        if (line.empty())
        {
            os  << " addr2line failed";
        }
        else if (line == "??:0" || line == "??:?" )
        {
            line = filename;
            os  << " in " << shorterPath(line).c_str();
        }
        else
        {
            os  << " at " << shorterPath(line).c_str();
        }
    }
}


// Uses 'which' to find executable on PATH
// - could also iterate through PATH directly
inline Foam::fileName whichPath(const char* fn)
{
    Foam::fileName fname(fn);

    if (!fname.empty() && fname[0] != '/' && fname[0] != '~')
    {
        std::string s = pipeOpen("which " + fname);

        if (s[0] == '/' || s[0] == '~')
        {
            fname = s;
        }
    }

    return fname;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::error::safePrintStack(std::ostream& os, int size)
{
    // Get raw stack symbols
    void *callstack[100];
    size = backtrace(callstack, (size > 0 && size < 100) ? size + 1 : 100);

    char **strings = backtrace_symbols(callstack, size);
    size_t rdelim;

    // Frame 0 is 'printStack()' - report something more meaningful
    os  << "[stack trace]" << std::endl
        << "=============" << std::endl;

    for (int i = 1; i < size; ++i)
    {
        std::string str(strings[i]);

        os  << '#' << i << '\t';

        // Possibly shorten paths that appear to correspond to OpenFOAM
        // locations (platforms).
        //
        // Eg, "/path/openfoam/platforms/linux64GccDPInt32Opt/lib/libxyz.so"
        // --> "platforms/linux64GccDPInt32Opt/lib/libxyz.so"

        auto ldelim = str.find('(');
        auto beg = str.find("/platforms/");

        if (beg == std::string::npos || !beg || beg > ldelim)
        {
            beg = 0;
        }
        else
        {
            ++beg;
        }

        if
        (
            (ldelim != std::string::npos)
         && (rdelim = str.find('+', ldelim+1)) != std::string::npos
         && (rdelim > ldelim+1)
        )
        {
            // Found function between () e.g. "(__libc_start_main+0xd0)"
            // - demangle function name (before the '+' offset)
            // - preserve trailing [0xAddr]

            os  << str.substr(beg, ldelim-beg)
                << ' '
                << demangleSymbol
                   (
                       str.substr(ldelim+1, rdelim-ldelim-1).c_str()
                   );

            if ((rdelim = str.find('[', rdelim)) != std::string::npos)
            {
                os  << ' ' << str.substr(rdelim);
            }
        }
        else if (beg)
        {
            // With shortened path name
            os  << str.substr(beg);
        }
        else
        {
            // No modification to string
            os  << str;
        }
        os  << std::endl;
    }

    os  << "=============" << std::endl;

    free(strings);
}


void Foam::error::printStack(Ostream& os, int size)
{
    // Get raw stack symbols
    void *callstack[100];
    size = backtrace(callstack, (size > 0 && size < 100) ? size + 1 : 100);

    Dl_info info;
    fileName fname;

    // Frame 0 is 'printStack()' - report something more meaningful
    os  << "[stack trace]" << nl
        << "=============" << nl;

    for (int i = 1; i < size; ++i)
    {
        int st = dladdr(callstack[i], &info);

        os  << '#' << i << "  ";
        if (st != 0 && info.dli_fname != nullptr && *(info.dli_fname))
        {
            fname = whichPath(info.dli_fname);

            if (info.dli_sname)
            {
                os  << demangleSymbol(info.dli_sname).c_str();
            }
            else
            {
                os  << '?';
            }
        }
        else
        {
            fname = "???";
            os  << '?';
        }

        printSourceFileAndLine(os, fname, info, callstack[i]);
        os  << nl;
    }

    os  << "=============" << nl;
}


// ************************************************************************* //
