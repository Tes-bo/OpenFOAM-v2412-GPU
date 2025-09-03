/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Functions used by OpenFOAM that are specific to POSIX compliant
    operating systems and need to be replaced or emulated on other systems.

SourceFiles
    POSIX.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OSspecific_H
#define Foam_OSspecific_H

#include "fileNameList.H"
#include "stringList.H"

#include <sys/types.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class CStringList;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return the PID of this process
pid_t pid();

//- Return the parent PID of this process
pid_t ppid();

//- Return the group PID of this process
pid_t pgid();

//- True if environment variable of given name is defined.
//  Using an empty name is a no-op and always returns false.
bool hasEnv(const std::string& envName);

//- Get environment value for given envName.
//  \return empty string if environment is undefined or envName is empty.
string getEnv(const std::string& envName);

//- Set an environment variable, return true on success.
//  Using an empty name is a no-op and always returns false.
bool setEnv(const word& name, const std::string& value, const bool overwrite);

//- Deprecated(2020-05) check for existence of environment variable
//  \deprecated(2020-05) - use hasEnv() function
FOAM_DEPRECATED_FOR(2020-05, "hasEnv() function")
inline bool env(const std::string& envName) { return Foam::hasEnv(envName); }

//- Return the system's host name, as per hostname(1)
string hostName();

//- Deprecated(2022-01) full hostname resolution may be unreliable
//  \deprecated(2022-01) - use hostname() function without parameter
FOAM_DEPRECATED_FOR(2022-01, "hostname() function without parameter")
string hostName(bool full);

//- Deprecated(2022-01) domain name resolution may be unreliable
//  \deprecated(2022-01) - avoid usage entirely
FOAM_DEPRECATED_FOR(2022-01, "nothing : avoid using at all")
string domainName();

//- Return the user's login name
string userName();

//- Is the current user the administrator (root)
bool isAdministrator();

//- Return home directory path name for the current user
fileName home();

//- Return home directory path name for a particular user
fileName home(const std::string& userName);

//- The physical or logical current working directory path name.
//  The behaviour is controlled by the \c cwd optimisation Switch
//  A value of 0 corresponds to the physical value, which is identical
//  to what getcwd and pwd -P would deliver.
//  A value of 1 corresponds to the logical value, which corresponds
//  to the PWD environment value and to what pwd -L would deliver.
fileName cwd();

//- The physical or logical current working directory path name.
fileName cwd(bool logical);

//- Change current directory to the one specified and return true on success.
//  Using an empty name is a no-op and always returns false.
bool chDir(const fileName& dir);

//- Make a directory and return an error if it could not be created
//  and does not already exist.
//  Using an empty pathName is a no-op and always returns false.
bool mkDir(const fileName& pathName, mode_t mode=0777);

//- Set the file/directory mode, return true on success.
//  Using an empty name is a no-op and always returns false.
bool chMod(const fileName& name, const mode_t mode);

//- Return the file mode, normally following symbolic links
//  Using an empty name is a no-op and always returns 0.
mode_t mode(const fileName& name, const bool followLink=true);

//- Return the file type: DIRECTORY or FILE, normally following symbolic links
//  Using an empty name is a no-op and always returns UNDEFINED.
fileName::Type type(const fileName& name, const bool followLink=true);

//- Does the name exist (as DIRECTORY or FILE) in the file system?
//  Optionally enable/disable check for gzip file.
//  Using an empty name is a no-op and always returns false.
bool exists
(
    const fileName& name,
    const bool checkGzip = true,
    const bool followLink = true
);

//- Does the name exist as a DIRECTORY in the file system?
//  Using an empty name is a no-op and always returns false.
bool isDir(const fileName& name, const bool followLink=true);

//- Does the name exist as a FILE in the file system?
//  Optionally enable/disable check for gzip file.
//  Using an empty name is a no-op and always returns false.
bool isFile
(
    const fileName& name,
    const bool checkGzip = true,
    const bool followLink = true
);

//- Return size of file or -1 on failure (normally follows symbolic links).
//  Using an empty name is a no-op and always returns -1.
off_t fileSize(const fileName& name, const bool followLink=true);

//- Return time of last file modification (normally follows symbolic links).
//  Using an empty name is a no-op and always returns 0.
time_t lastModified(const fileName& name, const bool followLink=true);

//- Return time of last file modification
//  Using an empty name is a no-op and always returns 0.
double highResLastModified(const fileName&, const bool followLink=true);

//- Read a directory and return the entries as a fileName List.
//  Using an empty directory name returns an empty list.
fileNameList readDir
(
    //! The directory name to scan
    const fileName& directory,
    //! The content type (eg, FILE | DIRECTORY)
    const fileName::Type type = fileName::Type::FILE,
    //! Trim '.gz' from file names
    const bool filtergz = true,
    //! Check what the symlink points to, not the symlink itself
    const bool followLink = true
);

//- Copy the source to the destination (recursively if necessary).
//  An empty source name is a no-op and always returns false.
bool cp(const fileName& src, const fileName& dst, const bool followLink=true);

//- Create a softlink. dst should not exist. Returns true if successful.
//  An empty source or destination name is a no-op that always returns false,
//  but also produces a warning.
bool ln(const fileName& src, const fileName& dst);

//- Return the contents (target) of a symlink.
//  For example, returns \c target for <code>ln -s target linkName</code>
//
//  An empty link name is a no-op that always returns an empty string.
//  Following a non-link is an error that returns an empty string.
fileName readLink(const fileName& link);

//- Rename src to dst.
//  An empty source or destination name is a no-op that always returns false.
bool mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink=false
);

//- Rename to a corresponding backup file
//  If the backup file already exists, attempt with "01" .. "99" suffix
//  An empty name or extension is a no-op that always returns false.
bool mvBak(const fileName& src, const std::string& ext = "bak");

//- Remove a file (or its gz equivalent), returning true if successful.
//  An empty name is a no-op that always returns false.
bool rm(const fileName& file);

//- Remove a directory and its contents recursively,
//  returning true if successful.
//  An empty directory name is a no-op that always returns false (silently)
bool rmDir
(
    //! The directory name to remove
    const fileName& directory,
    //! Silence warnings
    const bool silent = false,
    //! Only remove empty directories (recursive)
    const bool emptyOnly = false
);

//- Sleep for the specified number of seconds
unsigned int sleep(const unsigned int sec);

//- Close file descriptor
void fdClose(const int fd);

//- Check if machine is up by pinging given port
bool ping(const std::string& destName, const label port, const label timeOut);

//- Check if machine is up by pinging port 22 (ssh) and 222 (rsh)
bool ping(const std::string& host, const label timeOut=10);

//- Execute the specified command via the shell.
//  Uses vfork/execl internally.
//  When Foam::infoDetailLevel is zero, redirects stdout to stderr.
//
//  Where possible, use the list version instead.
//
//  \param bg return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
//      When Foam::infoDetailLevel is zero, redirects stdout to stderr.
int system(const std::string& command, const bool bg = false);

//- Execute the specified command with arguments.
//  Uses vfork/execvp internally
//  When Foam::infoDetailLevel is zero, redirects stdout to stderr.
//
//  \param bg return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
int system(const UList<string>& command, const bool bg = false);

//- Execute the specified command with arguments.
//  Uses vfork/execvp internally
//  When Foam::infoDetailLevel is zero, redirects stdout to stderr.
//
//  \param bg return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
int system(const CStringList& command, const bool bg = false);

//- Open a shared library and return handle to library.
//  A leading "lib" and ".so" suffix are added silently as required.
//  Prints warning if a library cannot be loaded (suppress with check=false)
void* dlOpen(const fileName& libName, const bool check=true);

//- Open a shared library and return handle to library.
//  A leading "lib" and ".so" suffix are added silently as required.
//  Captures any error messages in the second parameter (empty if no errors).
void* dlOpen(const fileName& libName, std::string& errorMsg);

//- Open shared libraries and return number of libraries loaded.
//  A leading "lib" and ".so" suffix are added silently as required.
//  Prints warning if a library cannot be loaded (suppress with check=false)
label dlOpen(std::initializer_list<fileName> libNames, const bool check=true);

//- Close a dlopened library using handle. Return true if successful
bool dlClose(void* handle);


//- Look for symbol in a dlopened library.
//  If the symbol is not 'required', using a null handle or an empty symbol
//  name is a no-op and returns nullptr without error.
//
//  \return the symbol or nullptr.
void* dlSymFind(void* handle, const std::string& symbol, bool required=false);

//- Lookup a symbol in a dlopened library using handle to library
//  \return the symbol or nullptr.
inline void* dlSym(void* handle, const std::string& symbol)
{
    return dlSymFind(handle, symbol, true);
}

//- Check for symbol in a dlopened library.
//  Using a null handle or an empty symbol name is a no-op and always
//  returns nullptr.
//  \return the symbol or nullptr.
inline void* dlSymFound(void* handle, const std::string& symbol)
{
    return dlSymFind(handle, symbol, false);
}


//- Return all loaded libraries
fileNameList dlLoaded();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
