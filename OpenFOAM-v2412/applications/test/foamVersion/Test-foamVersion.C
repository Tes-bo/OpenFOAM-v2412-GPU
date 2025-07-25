/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Test-foamVersion

Description
    Print the OpenFOAM version information.

\*---------------------------------------------------------------------------*/

#include <typeinfo>
#include "foamVersion.H"
#include "Switch.H"
#include "IOstreams.H"

using namespace Foam;


// Test extraction
void testExtraction(const std::string& str)
{
    Info<< "Extract: " << str << " =>"
        << " label: " << foamVersion::labelByteSize(str) << " bytes"
        << " scalar: " << foamVersion::scalarByteSize(str) << " bytes"
        << nl;
}


int main()
{
    Info<< "\nVersion information (function)" << nl;
    foamVersion::printBuildInfo(Info.stdStream());

    Info
        << "\nVersion information (macros)" << nl
        << "version   " << Foam::FOAMversion << nl
        << "build     " << Foam::FOAMbuild << nl
        << "buildArch " << Foam::FOAMbuildArch << nl;

    Info
        << "\nVersion information (namespace)" << nl
        << "patched?  = " << Switch(foamVersion::patched()) << nl
        << "api       " << foamVersion::api << nl
        << "patch     " << foamVersion::patch << nl
        << "version   " << foamVersion::version << nl
        << "build     " << foamVersion::build << nl
        << "buildArch " << foamVersion::buildArch << nl;

    Info
        << "\nTypes" << nl
        << "version   " << typeid(foamVersion::version).name() << nl
        << "build     " << typeid(foamVersion::build).name() << nl
        << "buildArch " << typeid(foamVersion::buildArch).name() << nl
        << "FOAMversion " << typeid(Foam::FOAMversion).name() << nl
        << "FOAMbuild   " << typeid(Foam::FOAMbuild).name() << nl;

    Info
        << "\nVerify memory addresses are identical:" << nl
        << "macro     " << name(Foam::FOAMversion) << nl
        << "namespace " << name(&(foamVersion::version[0])) << nl;


    // Test extraction
    {
        Info<< "\nTest size extraction routines" << nl;

        for
        (
            const std::string& str :
            {
                "MSB;label=32;scalar=64",
                "LSB;label=64;scalar=32",
                "LSB;label=;scalar=junk",
                "LSB;label==;scalar=128",
                "",
                "LSB;label;scalar",
                "LSB label=32 scalar=64",
            }
        )
        {
            testExtraction(str);
        }
    }

    Info
        << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
