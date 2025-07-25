/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Test atomic copyFile as per timeActivatedFileUpdate

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "Switch.H"

using namespace Foam;

#ifdef _WIN32
#undef DebugInfo    // Windows name clash with OpenFOAM messageStream
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

unsigned maxPath = MAX_PATH;

// Prefix "\\?\" for extended-length path and widen
// The prefix is only valid with absolute paths
//
// In the future, this code will be relocated in MSwindows.C

static inline std::wstring longName(const fileName& file)
{
    const auto len = file.length();

    std::wstring out;
    out.reserve(4 + len);

    if (len > maxPath)
    {
        if (file.isAbsolute())
        {
            out.append(L"\\\\?\\");
        }
        else
        {
            Warning
                << "Relative fileName exceeds " << maxPath
                << " characters" << nl
                << "    " << file << nl << nl;
        }
    }

    // Character-wise copy to get widening
    for (const auto c : file)
    {
        out += c;
    }

    return out;
}


bool win_mv(const fileName& src, const fileName& dst)
{
    constexpr const int flags
    (
        MOVEFILE_COPY_ALLOWED
      | MOVEFILE_REPLACE_EXISTING
      | MOVEFILE_WRITE_THROUGH
    );

    if (src.length() > maxPath || dst.length() > maxPath)
    {
        const std::wstring srcName(longName(src));
        const std::wstring dstName(longName(dst));

        Info<< "Windows mv (wide)" << nl;
        return ::MoveFileExW(srcName.c_str(), dstName.c_str(), flags);
    }

    Info<< "Windows mv (ansi)" << nl;
    return ::MoveFileExA(src.c_str(), dst.c_str(), flags);
}

#else

bool win_mv(const fileName& src, const fileName& dst)
{
    Info<< "No Windows mv, using Foam::mv" << nl;
    return Foam::mv(src, dst);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addNote("Test atomic copyFile methods");

    #ifdef _WIN32
    argList::addBoolOption("win1", "Foam cp, Windows mv");
    argList::addOption("maxPath", "length", "Test with shorter MAX_PATH");
    #endif

    argList::addArgument("srcFile");
    argList::addArgument("dstFile");

    #include "setRootCase.H"

    #ifdef _WIN32
    args.readIfPresent("maxPath", maxPath);
    #endif

    const auto srcFile = args.get<fileName>(1);
    const auto dstFile = args.get<fileName>(2);
    const fileName tmpFile(dstFile + Foam::name(pid()));

    Info<< "src   : " << srcFile << nl
        << "tmp   : " << tmpFile << nl
        << "dst   : " << dstFile << nl
        #ifdef _WIN32
        << "test with max-path: " << maxPath << nl
        #endif
        << nl;

    bool cpOk = false, mvOk = false;

    if (args.found("win1"))
    {
        const auto usage = argList::optionUsage.cfind("win1");
        if (usage.good())
        {
            Info<< "    " << (*usage).c_str() << nl;
        }

        cpOk = Foam::cp(srcFile, tmpFile);
        mvOk = win_mv(tmpFile, dstFile);
    }
    else
    {
        Info<< "    Foam cp, Foam mv" << nl;

        cpOk = Foam::cp(srcFile, tmpFile);
        mvOk = Foam::mv(tmpFile, dstFile);
    }

    Info<< nl
        << "cp: " << Switch(cpOk) << nl
        << "mv: " << Switch(mvOk) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
