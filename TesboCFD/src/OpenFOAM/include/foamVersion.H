/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Namespace
    Foam::foamVersion

Description
    Namespace for OpenFOAM version information

Note
    Compile-time version information is conveyed by the \b OPENFOAM define
    provided in the wmake rules "General/general"
    and compile-time configuration of some paths via FOAM_EXTRA_CXXFLAGS:
      - FOAM_CONFIGURED_PROJECT_DIR
      - FOAM_CONFIGURED_PROJECT_ETC

    For example,
    \verbatim
    FOAM_EXTRA_CXXFLAGS='-DFOAM_CONFIGURED_PROJECT_ETC=\"/etc/openfoam\"'
    \endverbatim

    The foamVersion.H file is located directly in the src/OpenFOAM/include
    directory for easier use by external packages and to allow easier
    modification during packaging.

SourceFiles
    foamConfig.Cver

\*---------------------------------------------------------------------------*/

#ifndef Foam_foamVersion_H
#define Foam_foamVersion_H

#include <iostream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- The directory name for user-resources within the HOME directory
//
//  Default: ".OpenFOAM"
//  Used by foamEtcFiles(), stringOps::expand(), Foam::JobInfo
#define FOAM_RESOURCE_USER_CONFIG_DIRNAME ".OpenFOAM"

//- The env name for site-resources to obtain a site-resources directory.
//
//  Default: "WM_PROJECT_SITE"
//  Used by foamEtcFiles() and stringOps::expand()
#define FOAM_RESOURCE_SITE_ENVNAME "WM_PROJECT_SITE"

//- The env name for determining a fallback directory name for site-resources
//- when the directory corresponding to FOAM_RESOURCE_SITE_ENVNAME is empty.
//  The fallback search appends "/site" to the directory.
//
//  Default: "WM_PROJECT_DIR"
//  Used by foamEtcFiles() and stringOps::expand()
#define FOAM_RESOURCE_SITE_FALLBACK_ENVNAME "WM_PROJECT_DIR"

// Fallback project directory name (hard-coded)
//
//  Default: undefined
//  Used by foamEtcFiles()
/* #undef FOAM_CONFIGURED_PROJECT_DIR */

// Fallback project etc/ directory name (hard-coded)
//
//  Default: undefined
//  Used by foamEtcFiles()
/* #undef FOAM_CONFIGURED_PROJECT_ETC */


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Version information
    namespace foamVersion
    {
        //- OpenFOAM api number (integer) corresponding to the value of OPENFOAM
        //- at the time of compilation.
        //  The value is 0 if OPENFOAM was not defined.
        extern const int api;

        //- OpenFOAM patch number as a std::string
        extern const std::string patch;

        //- OpenFOAM build information as a std::string
        extern const std::string build;

        //- OpenFOAM build architecture information
        //- (machine endian, label/scalar sizes) as a std::string
        extern const std::string buildArch;

        //- OpenFOAM version (name or stringified number) as a std::string
        extern const std::string version;

        //- Test if the patch string appears to be in use,
        //- which is when it is defined (non-zero).
        bool patched();

        //- Extract label size (in bytes) from "label=" tag in string
        unsigned labelByteSize(const std::string& str);

        //- Extract scalar size (in bytes) from "scalar=" tag in string
        unsigned scalarByteSize(const std::string& str);

        //- Print information about version, build, arch to output stream
        //
        //  Eg,
        //  \code
        //      Using: OpenFOAM-<VER> (API) - visit www.openfoam.com
        //      Build: <BUILD> (patch=...)
        //      Arch:  <ARCH_INFO>
        //  \endcode
        //
        //  \param os   the output stream
        //  \param full includes Arch information
        void printBuildInfo(std::ostream& os, const bool full=true);

        //- Compile-time definition of the OpenFOAM project directory
        //- or empty if not defined.
        //  Functional equivalent to WM_PROJECT_DIR.
        std::string configuredProjectDir();

        //- Compile-time definition of the OpenFOAM etc/ directory
        //- or empty if not defined.
        //  Functional equivalent to WM_PROJECT_DIR/etc
        std::string configuredEtcDir();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compatibility names (1806 and earlier).
// Historically assumed to be called within the Foam namespace, they are
// thus defined without the Foam namespace qualifier.
//
// - FOAMversion: c-string
// - FOAMbuild:   c-string
// - FOAMbuildArch: std::string

#define FOAMversion     foamVersion::version.c_str()
#define FOAMbuild       foamVersion::build.c_str()
#define FOAMbuildArch   foamVersion::buildArch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
