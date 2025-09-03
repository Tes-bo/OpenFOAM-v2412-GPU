/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "etcFiles.H"
#include "foamVersion.H"
#include "OSspecific.H"

#ifdef FULLDEBUG
#ifndef FOAM_RESOURCE_USER_CONFIG_DIRNAME
# warning FOAM_RESOURCE_USER_CONFIG_DIRNAME undefined (was this intentional?)
#endif

#ifndef FOAM_RESOURCE_SITE_ENVNAME
# warning FOAM_RESOURCE_SITE_ENVNAME undefined (was this intentional?)
#endif

#ifndef FOAM_RESOURCE_SITE_FALLBACK_ENVNAME
# warning FOAM_RESOURCE_SITE_FALLBACK_ENVNAME undefined (was this intentional?)
#endif
#endif

// Always use these names
#undef  FOAM_PROJECT_ENVNAME
#define FOAM_PROJECT_ENVNAME "WM_PROJECT_DIR"

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

//
// Some of these could be exposed too (if required),
// but are fairly special purpose.
//

namespace
{

// Return the file location mode as a string.
//
// - u : location mask 0700
// - g : location mask 0070
// - o : location mask 0007
static inline std::string locationToString(unsigned short location)
{
    std::string mode;

    if (location & 0700) { mode += 'u'; } // User
    if (location & 0070) { mode += 'g'; } // Group
    if (location & 0007) { mode += 'o'; } // Other
    if (mode.empty()) { mode = "???"; }

    return mode;
}


// Error handling when a mandatory entry is not found
static inline void errorMandatoryNotFound
(
    const std::string& name,
    unsigned short location
)
{
    // Abort when mandatory entry was not found.
    // Use a direct exit, since this could occur before anything is
    // setup at all.

    std::cerr
        << "--> FOAM FATAL ERROR :\n"
        "    Could not find mandatory etc entry (mode="
        << locationToString(location) << ")\n    '"
        << name << "'\n"
        << std::endl;
    std::exit(1);
}


// Assign 'queried' parameter to the user resource directory.
// Return true if this directory exists.
//
// Corresponds to foamEtcFile -mode=u
// Looks for
//   - ~/.OpenFOAM
static inline bool userResourceDir(Foam::fileName& queried)
{
    #ifdef FOAM_RESOURCE_USER_CONFIG_DIRNAME
    queried = Foam::home()/FOAM_RESOURCE_USER_CONFIG_DIRNAME;
    if (Foam::isDir(queried))
    {
        // If home() fails, it will have actually queried "./.OpenFOAM"
        // instead.
        // But we would have worse problems elsewhere if that were the case.
        return true;
    }
    #endif

    return false;
}


// Assign 'queried' parameter to the group resource directory.
// Return true if this directory exists.
// Otherwise clears the parameter and returns false.
//
// Corresponds to foamEtcFile -mode=g
// Looks for
//   - ${WM_PROJECT_SITE}/etc
//   - ${WM_PROJECT_DIR}/site/etc
//
// Optionally (compile-time defined):
//   - FOAM_CONFIGURED_PROJECT_DIR/site/etc
static inline bool groupResourceDir(Foam::fileName& queried)
{
    #ifdef FOAM_RESOURCE_SITE_ENVNAME
    queried = Foam::getEnv(FOAM_RESOURCE_SITE_ENVNAME)/"etc";
    if (queried.size() > 4)
    {
        return Foam::isDir(queried);
    }
    #endif

    // Fallback when WM_PROJECT_SITE is unset

    #ifdef FOAM_RESOURCE_SITE_FALLBACK_ENVNAME
    queried = Foam::getEnv(FOAM_RESOURCE_SITE_FALLBACK_ENVNAME)/"site/etc";
    if (queried.size() > 9 && Foam::isDir(queried))
    {
        return true;
    }
    #endif

    // Compile-time paths
    #ifdef FOAM_CONFIGURED_PROJECT_DIR
    queried = FOAM_CONFIGURED_PROJECT_DIR "/site/etc";
    if (queried.size() > 9 && Foam::isDir(queried))
    {
        return true;
    }
    #endif

    queried.clear();
    return false;
}


// Assign 'queried' parameter to the OpenFOAM etc/ resource directory.
// Return true if it exists.
// Otherwise clears the parameter and returns false.
//
// Corresponds to foamEtcFile -mode=o
// Looks for
//   - ${WM_PROJECT_DIR}/etc
//
// Optionally (compile-time defined):
//   - FOAM_CONFIGURED_PROJECT_ETC
//   - FOAM_CONFIGURED_PROJECT_DIR/"etc"
static inline bool projectResourceDir(Foam::fileName& queried)
{
    queried = Foam::getEnv(FOAM_PROJECT_ENVNAME)/"etc";
    if (queried.size() > 4 && Foam::isDir(queried))
    {
        return true;
    }

    // Compile-time paths

    #ifdef FOAM_CONFIGURED_PROJECT_ETC
    queried = FOAM_CONFIGURED_PROJECT_ETC;
    if (Foam::isDir(queried))
    {
        return true;
    }
    #endif

    #ifdef FOAM_CONFIGURED_PROJECT_DIR
    queried = FOAM_CONFIGURED_PROJECT_DIR "/etc";
    if (queried.size() > 4 && Foam::isDir(queried))
    {
        return true;
    }
    #endif

    queried.clear();
    return false;
}


// Check if the named file/directory matches the type required.
//
// - typeRequired (UNDEFINED)        => accept either FILE or DIRECTORY
// - typeRequired (FILE | DIRECTORY) => accept only that type
static inline bool accept
(
    const Foam::fileName& name,
    const Foam::fileName::Type typeRequired
)
{
    // followLink(true), checkGzip(true)
    // -> returns (UNDEFINED | FILE | DIRECTORY), no need to check for (SYMLINK)
    const auto t = name.type(true, true);

    return
    (
        // Found something?
        Foam::fileName::Type::UNDEFINED != t
    &&
        (
            // Any particular type required?
            Foam::fileName::Type::UNDEFINED == typeRequired
          ? (Foam::fileName::Type::UNDEFINED != t)
          : (typeRequired == t)
        )
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileNameList Foam::etcDirs(bool test)
{
    // Use foamVersion::api (instead of the OPENFOAM define) to ensure this
    // stays properly synchronized with the build information
    const Foam::fileName version(std::to_string(Foam::foamVersion::api));

    Foam::fileNameList list(5);
    Foam::fileName queried;
    label nDirs = 0;

    // User resource directories
    if (userResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried/version;
        list[nDirs++] = queried;
    }

    // Group (site) resource directories
    if (groupResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried/version;
        list[nDirs++] = queried;
    }

    // Other (project) resource directory
    if (projectResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried;
    }

    list.resize(nDirs);

    return list;
}


Foam::fileNameList Foam::findEtcEntries
(
    const Foam::fileName& name,
    unsigned short location,
    const Foam::fileName::Type typeRequired,
    const bool findFirst
)
{
    // Debug Tracing
    // std::cerr
    //     << "search ("<< locationToString(location) << "): "
    //     << name.c_str() << '\n';

    if (!(location & 0777))
    {
        // Warn about bad location (mode) ... or make it FATAL?
        std::cerr
            << "--> FOAM Error :\n    "
            "No user/group/other location specified for 'etc' file"
            " or directory\n    '"
            << name.c_str() << "'\n\n" << std::endl;
    }

    // Use foamVersion::api (instead of the OPENFOAM define) to ensure this
    // stays properly synchronized with the build information
    const Foam::fileName version(std::to_string(Foam::foamVersion::api));

    Foam::fileNameList list;
    Foam::fileName queried, candidate;

    if (fileName::Type::FILE == typeRequired && name.empty())
    {
        // FILE must have a name to be found!
        return list;
    }


    // User resource directories
    if ((location & 0700) && userResourceDir(queried))
    {
        candidate = queried/version/name;
        if (accept(candidate, typeRequired))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = queried/name;
        if (accept(candidate, typeRequired))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }
    }


    // Group (site) resource directories
    if ((location & 0070) && groupResourceDir(queried))
    {
        candidate = queried/version/name;
        if (accept(candidate, typeRequired))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = queried/name;
        if (accept(candidate, typeRequired))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }
    }


    // Other (project) resource directory
    if ((location & 0007) && projectResourceDir(queried))
    {
        candidate = queried/name;
        if (accept(candidate, typeRequired))
        {
            list.append(std::move(candidate));
        }
    }

    return list;
}


Foam::fileNameList Foam::findEtcDirs
(
    const fileName& name,
    unsigned short location,
    const bool findFirst
)
{
    return findEtcEntries(name, location, fileName::Type::DIRECTORY, findFirst);
}


Foam::fileNameList Foam::findEtcFiles
(
    const fileName& name,
    const bool mandatory,
    unsigned short location,
    const bool findFirst
)
{
    // Note: findEtcEntries checks name.size() for FILE
    fileNameList list
    (
        findEtcEntries(name, location, fileName::Type::FILE, findFirst)
    );

    if (mandatory && list.empty())
    {
        errorMandatoryNotFound(name, location);
    }

    return list;
}


Foam::fileName Foam::findEtcEntry
(
    const fileName& name,
    unsigned short location,
    const Foam::fileName::Type typeRequired
)
{
    // With findFirst(true)
    fileNameList list(findEtcEntries(name, location, typeRequired, true));

    fileName found;

    if (list.size())
    {
        found = std::move(list.first());
    }

    return found;
}


Foam::fileName Foam::findEtcDir
(
    const fileName& name,
    unsigned short location
)
{
    return findEtcEntry(name, location, fileName::Type::DIRECTORY);
}


Foam::fileName Foam::findEtcFile
(
    const fileName& name,
    const bool mandatory,
    unsigned short location
)
{
    fileName found(findEtcEntry(name, location, fileName::Type::FILE));

    if (mandatory && found.empty())
    {
        errorMandatoryNotFound(name, location);
    }

    return found;
}


// ************************************************************************* //
