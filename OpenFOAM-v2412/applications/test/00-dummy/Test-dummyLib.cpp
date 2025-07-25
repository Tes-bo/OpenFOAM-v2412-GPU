/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Minimal compilation test with wmake, without OpenFOAM libraries.

    The application and library can also serve as a minimal test case for
    wmake, or to provide a minimal library/executable target for testing.

\*---------------------------------------------------------------------------*/

#include "dummyLib.H"

#include <cstring>
#include <iostream>

constexpr char nl = '\n';
constexpr const char* const bold = "\\fB";  // nroff
constexpr const char* const norm = "\\fR";  // nroff

constexpr const char* const website = "www.openfoam.com";

using std::cout;
using dummyLib = Foam::Detail::dummyLib;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static void printMan(const char* exeName)
{
    cout
        << ".TH \"" << exeName << "\" 1 "
        << "\"OpenFOAM-v" << OPENFOAM
        << "\" \"" << website << "\" \"OpenFOAM Commands Manual\""
        << nl;

    cout
        << ".SH NAME" << nl
        << exeName
        << " \\- part of " << bold << "OpenFOAM" << norm
        << " (The Open Source CFD Toolbox)." << nl
        << ".SH SYNOPSIS" << nl
        << bold << exeName << norm << " [OPTIONS]" << nl;

    cout
        << ".SH DESCRIPTION" << nl
        << ".nf" << nl
        << "Minimal compilation test with wmake, without OpenFOAM libraries."
        << nl
        << ".fi" << nl;

    cout
        << ".SH OPTIONS" << nl
        << ".TP" << nl
        << "-parallel" << nl
        << "Run parallel and provide simple report" << nl;

    if (!Foam::Detail::dummyLib::hasMPI())
    {
        cout << "[warning: no mpi]" << nl;
    }

    cout
        << ".TP" << nl
        << "-help-man" << nl
        << "Display manpage" << nl;

    cout
        << ".SH INFORMATION" << nl
        << ".nf" << nl
        << "label    = " << dummyLib::label_size << nl
        << "scalar   = " << dummyLib::scalar_size;

    if
    (
        dummyLib::solveScalar_size
     && dummyLib::solveScalar_size != dummyLib::scalar_size
    )
    {
        cout
            << " [solve=" << dummyLib::solveScalar_size << "]";
    }
    cout
        << " (" << dummyLib::precision << ')' << nl
        << "arch     = " << dummyLib::arch << nl
        << "compiler = " << dummyLib::compiler << nl;

    cout
        << nl
        << "archComp     = " << dummyLib::archComp << nl
        << "archCompBase = " << dummyLib::archCompBase << nl
        << "archCompFull = " << dummyLib::archCompFull << nl;
    cout
        << ".fi" << nl;

    cout
        << ".SH \"SEE ALSO\"" << nl
        << "Online documentation https://" << website << "/documentation/"
        << nl;
}


static void printInfo()
{
    cout
        << nl
        << "OPENFOAM  = " << OPENFOAM << nl
        << "label     = " << dummyLib::label_size << nl
        << "scalar    = " << dummyLib::scalar_size
        << " (" << dummyLib::precision << ')' << nl;

    if
    (
        dummyLib::solveScalar_size
     && dummyLib::solveScalar_size != dummyLib::scalar_size
    )
    {
        cout
            << "solve     = " << dummyLib::solveScalar_size << nl;
    }

    cout
        << "arch      = " << dummyLib::arch << nl
        << "compiler  = " << dummyLib::compiler << nl;

    cout
        << nl
        << "archComp     = " << dummyLib::archComp << nl
        << "archCompBase = " << dummyLib::archCompBase << nl
        << "archCompFull = " << dummyLib::archCompFull << nl;

    cout<< nl;
}


int main(int argc, char *argv[])
{
    bool master = true;

    if (argc > 1)
    {
        if (strcmp(argv[1], "-help-man") == 0)
        {
            printMan("Test-dummyLib");
            return 0;
        }

        if (strcmp(argv[1], "-parallel") == 0)
        {
            master = dummyLib::printMPI();
        }
    }

    if (master)
    {
        printInfo();
    }

    return 0;
}


// ************************************************************************* //
