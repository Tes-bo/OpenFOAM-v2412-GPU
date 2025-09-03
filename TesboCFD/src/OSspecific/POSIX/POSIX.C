/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    POSIX versions of the functions declared in OSspecific.H

\*---------------------------------------------------------------------------*/

#if defined(__sun__) && defined(__GNUC__)
    // Not certain if this is still required
    #define _SYS_VNODE_H
#endif

#include "OSspecific.H"
#include "POSIX.H"
#include "fileName.H"
#include "fileStat.H"
#include "timer.H"
#include "DynamicList.H"
#include "CStringList.H"
#include "stringOps.H"
#include "IOstreams.H"
#include "Pstream.H"

#include <fstream>
#include <cstdlib>
#include <cctype>

#include <cstdio>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <dlfcn.h>

#ifdef __APPLE__
    #define EXT_SO  "dylib"
    #include <mach-o/dyld.h>
#else
    #define EXT_SO  "so"

    // PGI does not have __int128_t
    #ifdef __PGIC__
        #define __ILP32__
    #endif

    #include <link.h>
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(POSIX, 0);
}

static bool cwdPreference_(Foam::debug::optimisationSwitch("cwd", 0));


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// After a fork in system(), before the exec() do the following
// - close stdin when executing in background (daemon-like)
// - redirect stdout to stderr when infoDetailLevel == 0
static inline void redirects(const bool bg)
{
    if (bg)
    {
        // Close stdin(0) - unchecked return value
        (void) ::close(STDIN_FILENO);
    }

    // Redirect stdout(1) to stderr(2) '1>&2'
    if (Foam::infoDetailLevel == 0)
    {
        // This is correct.  1>&2 means dup2(2, 1);
        (void) ::dup2(STDERR_FILENO, STDOUT_FILENO);
    }
}


// Library loading is normally simply via dlopen(),
// but SIP (System Integrity Protection) on Apple will generally
// clear out the DYLD_LIBRARY_PATH set from shell scripts.
// We thus have FOAM_LD_LIBRARY_PATH as a shadow parameter and use
// that to attempt loading ourselves
static inline void* loadLibrary(const Foam::fileName& libName)
{
    constexpr int ldflags = (RTLD_LAZY|RTLD_GLOBAL);

#ifdef __APPLE__
    const char* normal = nullptr;
    const char* shadow = nullptr;

    if
    (
        !libName.isAbsolute()
     && ((normal = ::getenv("DYLD_LIBRARY_PATH"))    == nullptr || !*normal)
     && ((shadow = ::getenv("FOAM_LD_LIBRARY_PATH")) != nullptr && *shadow)
    )
    {
        // SIP appears to have cleared DYLD_LIBRARY_PATH but the
        // shadow parameter is available

        const Foam::string ldPaths(shadow);
        const auto paths = Foam::stringOps::split<Foam::string>(ldPaths, ':');

        for (const auto& p : paths)
        {
            if (p.length())  // Split removes empty, but be paranoid
            {
                const Foam::fileName fullPath(p.str()/libName);
                void* handle = ::dlopen(fullPath.c_str(), ldflags);
                if (handle)
                {
                    return handle;
                }
            }
        }
    }
#endif

    // Regular loading
    return ::dlopen(libName.c_str(), ldflags);
}


// * * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * //

namespace Foam
{
namespace POSIX
{

//- A simple directory contents iterator
class directoryIterator
{
    DIR* dirptr_;

    bool exists_;

    bool hidden_;

    std::string item_;

    //- Accept file/dir name
    inline bool accept() const
    {
        return
        (
            item_.size() && item_ != "." && item_ != ".."
         && (hidden_ || item_[0] != '.')
        );
    }


public:

    // Constructors

        //- Construct for dirName, optionally allowing hidden files/dirs
        directoryIterator(const std::string& dirName, bool allowHidden = false)
        :
            dirptr_(nullptr),
            exists_(false),
            hidden_(allowHidden),
            item_()
        {
            if (!dirName.empty())
            {
                dirptr_ = ::opendir(dirName.c_str());
                exists_ = (dirptr_ != nullptr);
                next(); // Move to first element
            }
        }


    //- Destructor
    ~directoryIterator()
    {
        close();
    }


    // Member Functions

        //- Directory open succeeded
        bool exists() const noexcept
        {
            return exists_;
        }

        //- Directory pointer is valid
        bool good() const noexcept
        {
            return dirptr_;
        }

        //- Close directory
        void close()
        {
            if (dirptr_)
            {
                ::closedir(dirptr_);
                dirptr_ = nullptr;
            }
        }

        //- The current item
        const std::string& val() const noexcept
        {
            return item_;
        }

        //- Read next item, always ignoring "." and ".." entries.
        //  Normally also ignore hidden files/dirs (beginning with '.')
        //  Automatically close when there are no more items
        bool next()
        {
            struct dirent *list;

            while (dirptr_ && (list = ::readdir(dirptr_)) != nullptr)
            {
                item_ = list->d_name;

                if (accept())
                {
                    return true;
                }
            }
            close(); // No more items

            return false;
        }


    // Member Operators

        //- Same as good()
        operator bool() const noexcept
        {
            return good();
        }

        //- Same as val()
        const std::string& operator*() const noexcept
        {
            return val();
        }

        //- Same as next()
        directoryIterator& operator++()
        {
            next();
            return *this;
        }
};

} // End namespace POSIX
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pid_t Foam::pid()
{
    return ::getpid();
}


pid_t Foam::ppid()
{
    return ::getppid();
}


pid_t Foam::pgid()
{
    return ::getpgrp();
}


bool Foam::hasEnv(const std::string& envName)
{
    // An empty envName => always false
    return !envName.empty() && ::getenv(envName.c_str()) != nullptr;
}


Foam::string Foam::getEnv(const std::string& envName)
{
    // Ignore an empty envName => always ""
    char* env = envName.empty() ? nullptr : ::getenv(envName.c_str());

    if (env)
    {
        return string(env);
    }

    // Return null-constructed string rather than string::null
    // to avoid cyclic dependencies in the construction of globals
    return string();
}


bool Foam::setEnv
(
    const word& envName,
    const std::string& value,
    const bool overwrite
)
{
    // Ignore an empty envName => always false
    return
    (
        !envName.empty()
     && ::setenv(envName.c_str(), value.c_str(), overwrite) == 0
    );
}


Foam::string Foam::hostName()
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));
    return buf;
}


// DEPRECATED (2022-01)
Foam::string Foam::hostName(bool full)
{
    // implementation as per hostname from net-tools
    if (full)
    {
        char buf[128];
        ::gethostname(buf, sizeof(buf));

        struct hostent *hp = ::gethostbyname(buf);
        if (hp)
        {
            return hp->h_name;
        }
        return buf;
    }

    return Foam::hostName();
}


// DEPRECATED (2022-01)
Foam::string Foam::domainName()
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));

    // implementation as per hostname from net-tools
    struct hostent *hp = ::gethostbyname(buf);
    if (hp)
    {
        char *p = ::strchr(hp->h_name, '.');
        if (p)
        {
            ++p;
            return p;
        }
    }

    return string();
}


Foam::string Foam::userName()
{
    struct passwd* pw = ::getpwuid(::getuid());
    if (pw != nullptr)
    {
        return pw->pw_name;
    }

    return string();
}


bool Foam::isAdministrator()
{
    return (::geteuid() == 0);
}


Foam::fileName Foam::home()
{
    char* env = ::getenv("HOME");
    if (env)
    {
        return fileName(env);
    }

    struct passwd* pw = ::getpwuid(::getuid());
    if (pw)
    {
        return pw->pw_dir;
    }

    return fileName();
}


Foam::fileName Foam::home(const std::string& userName)
{
    // An empty userName => same as home()
    if (userName.empty())
    {
        return Foam::home();
    }

    struct passwd* pw = ::getpwnam(userName.c_str());
    if (pw)
    {
        return pw->pw_dir;
    }

    return fileName();
}


namespace Foam
{

//- The physical current working directory path name (pwd -P).
static Foam::fileName cwd_P()
{
    label pathLengthLimit = POSIX::pathLengthChunk;
    List<char> path(pathLengthLimit);

    // Resize path if getcwd fails with an ERANGE error
    while (pathLengthLimit == path.size())
    {
        if (::getcwd(path.data(), path.size()))
        {
            return path.data();
        }
        else if (errno == ERANGE)
        {
            // Increment path length up to the pathLengthMax limit
            if
            (
                (pathLengthLimit += POSIX::pathLengthChunk)
             >= POSIX::pathLengthMax
            )
            {
                FatalErrorInFunction
                    << "Attempt to increase path length beyond limit of "
                    << POSIX::pathLengthMax
                    << exit(FatalError);
            }

            path.resize(pathLengthLimit);
        }
        else
        {
            break;
        }
    }

    FatalErrorInFunction
        << "Couldn't get the current working directory"
        << exit(FatalError);

    return fileName();
}


//- The logical current working directory path name.
// From the PWD environment, same as pwd -L.
static Foam::fileName cwd_L()
{
    const char* env = ::getenv("PWD");

    // Basic check
    if (!env || env[0] != '/')
    {
        WarningInFunction
            << "PWD is invalid - reverting to physical description"
            << nl;

        return cwd_P();
    }

    fileName dir(env);

    // Check for "/."
    for
    (
        std::string::size_type pos = 0;
        std::string::npos != (pos = dir.find("/.", pos));
        /*nil*/
    )
    {
        pos += 2;

        if
        (
            // Ends in "/." or has "/./"
            !dir[pos] || dir[pos] == '/'

            // Ends in "/.." or has "/../"
         || (dir[pos] == '.' && (!dir[pos+1] || dir[pos+1] == '/'))
        )
        {
            WarningInFunction
                << "PWD contains /. or /.. - reverting to physical description"
                << nl;

            return cwd_P();
        }
    }

    // Finally, verify that PWD actually corresponds to the "." directory
    if (!fileStat(dir, true).sameINode(fileStat(".", true)))
    {
        WarningInFunction
            << "PWD is not the cwd() - reverting to physical description"
            << nl;

        return cwd_P();
    }


    return fileName(dir);
}

} // End namespace Foam


Foam::fileName Foam::cwd()
{
    return cwd(cwdPreference_);
}


Foam::fileName Foam::cwd(bool logical)
{
    if (logical)
    {
        return cwd_L();
    }

    return cwd_P();
}


bool Foam::chDir(const fileName& dir)
{
    // Ignore an empty dir name => always false
    return !dir.empty() && ::chdir(dir.c_str()) == 0;
}


bool Foam::mkDir(const fileName& pathName, mode_t mode)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : pathName:" << pathName << " mode:" << mode
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // empty names are meaningless
    if (pathName.empty())
    {
        return false;
    }

    // Construct path directory if does not exist
    if (::mkdir(pathName.c_str(), mode) == 0)
    {
        // Directory made OK so return true
        return true;
    }

    switch (errno)
    {
        case EPERM:
        {
            FatalErrorInFunction
                << "The filesystem containing " << pathName
                << " does not support the creation of directories."
                << exit(FatalError);
            break;
        }

        case EEXIST:
        {
            // Directory already exists so simply return true
            return true;
        }

        case EFAULT:
        {
            FatalErrorInFunction
                << "" << pathName
                << " points outside your accessible address space."
                << exit(FatalError);
            break;
        }

        case EACCES:
        {
            FatalErrorInFunction
                << "The parent directory does not allow write "
                   "permission to the process,"<< nl
                << " or one of the directories in " << pathName
                << " did not allow search (execute) permission."
                << exit(FatalError);
            break;
        }

        case ENAMETOOLONG:
        {
            FatalErrorInFunction
                << "" << pathName << " is too long."
                << exit(FatalError);
            break;
        }

        case ENOENT:
        {
            // Part of the path does not exist so try to create it
            if (pathName.path().size() && mkDir(pathName.path(), mode))
            {
                return mkDir(pathName, mode);
            }

            FatalErrorInFunction
                << "Couldn't create directory " << pathName
                << exit(FatalError);
            break;
        }

        case ENOTDIR:
        {
            FatalErrorInFunction
                << "A component used as a directory in " << pathName
                << " is not, in fact, a directory."
                << exit(FatalError);
            break;
        }

        case ENOMEM:
        {
            FatalErrorInFunction
                << "Insufficient kernel memory was available to make directory "
                << pathName << '.'
                << exit(FatalError);
            break;
        }

        case EROFS:
        {
            FatalErrorInFunction
                << "" << pathName
                << " refers to a file on a read-only filesystem."
                << exit(FatalError);
            break;
        }

        case ELOOP:
        {
            FatalErrorInFunction
                << "Too many symbolic links were encountered in resolving "
                << pathName << '.'
                << exit(FatalError);
            break;
        }

        case ENOSPC:
        {
            FatalErrorInFunction
                << "The device containing " << pathName
                << " has no room for the new directory or "
                << "the user's disk quota is exhausted."
                << exit(FatalError);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Couldn't create directory " << pathName
                << exit(FatalError);
            break;
        }
    }

    return false;
}


bool Foam::chMod(const fileName& name, const mode_t m)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return !name.empty() && ::chmod(name.c_str(), m) == 0;
}


mode_t Foam::mode(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always 0
    if (!name.empty())
    {
        fileStat fileStatus(name, followLink);
        if (fileStatus.good())
        {
            return fileStatus.status().st_mode;
        }
    }

    return 0;
}


Foam::fileName::Type Foam::type
(
    const fileName& name,
    const bool followLink
)
{
    // Ignore an empty name => always UNDEFINED
    if (name.empty())
    {
        return fileName::Type::UNDEFINED;
    }

    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
    }

    mode_t m = mode(name, followLink);

    if (S_ISREG(m))
    {
        return fileName::Type::FILE;
    }
    else if (S_ISLNK(m))
    {
        return fileName::Type::SYMLINK;
    }
    else if (S_ISDIR(m))
    {
        return fileName::Type::DIRECTORY;
    }

    return fileName::Type::UNDEFINED;
}


bool Foam::exists
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkGzip:" << checkGzip
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return
    (
        !name.empty()
     && (mode(name, followLink) || isFile(name, checkGzip, followLink))
    );
}


bool Foam::isDir(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return !name.empty() && S_ISDIR(mode(name, followLink));
}


bool Foam::isFile
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkGzip:" << checkGzip
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return
    (
        !name.empty()
     && (
            S_ISREG(mode(name, followLink))
         || (checkGzip && S_ISREG(mode(name + ".gz", followLink)))
        )
    );
}


off_t Foam::fileSize(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    if (!name.empty())
    {
        fileStat fileStatus(name, followLink);
        if (fileStatus.good())
        {
            return fileStatus.status().st_size;
        }
    }

    return -1;
}


time_t Foam::lastModified(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    return name.empty() ? 0 : fileStat(name, followLink).modTime();
}


double Foam::highResLastModified(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    return name.empty() ? 0 : fileStat(name, followLink).dmodTime();
}


Foam::fileNameList Foam::readDir
(
    const fileName& directory,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
)
{
    // Initial filename list size and the increment when resizing the list
    constexpr int maxNnames = 100;

    fileNameList dirEntries;

    // Iterate contents (ignores an empty directory name)

    POSIX::directoryIterator dirIter(directory);
    if (!dirIter.exists())
    {
        if (POSIX::debug)
        {
            InfoInFunction
                << "cannot open directory " << directory << endl;
        }

        return dirEntries;
    }

    if (POSIX::debug)
    {
        // InfoInFunction
        Pout<< FUNCTION_NAME << " : reading directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    label nFailed = 0;     // Entries with invalid characters
    label nEntries = 0;    // Number of selected entries
    dirEntries.resize(maxNnames);

    // Process the directory entries
    for (/*nil*/; dirIter; ++dirIter)
    {
        const std::string& item = *dirIter;

        // Validate filename without spaces, quotes, etc in the name.
        // No duplicate slashes to strip - dirent will not have them anyhow.

        fileName name(fileName::validate(item));
        if (name != item)
        {
            ++nFailed;
        }
        else if
        (
            (type == fileName::Type::DIRECTORY)
         || (type == fileName::Type::FILE && !fileName::isBackup(name))
        )
        {
            fileName::Type detected = (directory/name).type(followLink);

            if (detected == type)
            {
                // Only strip '.gz' from non-directory names
                if
                (
                    filtergz
                 && (detected != fileName::Type::DIRECTORY)
                 && name.has_ext("gz")
                )
                {
                    name.remove_ext();
                }

                if (nEntries >= dirEntries.size())
                {
                    dirEntries.resize(dirEntries.size() + maxNnames);
                }

                dirEntries[nEntries] = std::move(name);
                ++nEntries;
            }
        }
    }

    // Finalize the length of the entries list
    dirEntries.resize(nEntries);

    if (nFailed && POSIX::debug)
    {
        std::cerr
            << "Foam::readDir() : reading directory " << directory << nl
            << nFailed << " entries with invalid characters in their name"
            << std::endl;
    }

    return dirEntries;
}


bool Foam::cp(const fileName& src, const fileName& dest, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : src:" << src << " dest:" << dest << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Make sure source exists - this also handles an empty source name
    if (!exists(src))
    {
        return false;
    }

    const fileName::Type srcType = src.type(followLink);

    fileName destFile(dest);

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile /= src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams. Enforce binary for extra safety
        std::ifstream srcStream(src, ios_base::in | ios_base::binary);
        if (!srcStream)
        {
            return false;
        }

        std::ofstream destStream(destFile, ios_base::out | ios_base::binary);
        if (!destStream)
        {
            return false;
        }

        // Copy character data.
        char ch;
        while (srcStream.get(ch))
        {
            destStream.put(ch);
        }

        // Final check.
        if (!srcStream.eof() || !destStream)
        {
            return false;
        }
    }
    else if (srcType == fileName::SYMLINK)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile /= src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        Foam::ln(src, destFile);
    }
    else if (srcType == fileName::DIRECTORY)
    {
        if (destFile.type() == fileName::DIRECTORY)
        {
            // Both are directories. Could mean copy contents or copy
            // recursively.  Don't actually know what the user wants,
            // but assume that if names are identical == copy contents.
            //
            // So: "path1/foo" "path2/foo"  copy contents
            // So: "path1/foo" "path2/bar"  copy directory

            const word srcDirName = src.name();
            if (destFile.name() != srcDirName)
            {
                destFile /= srcDirName;
            }
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        char* realSrcPath = realpath(src.c_str(), nullptr);
        char* realDestPath = realpath(destFile.c_str(), nullptr);
        const bool samePath = strcmp(realSrcPath, realDestPath) == 0;

        if (POSIX::debug && samePath)
        {
            InfoInFunction
                << "Attempt to copy " << realSrcPath << " to itself" << endl;
        }

        if (realSrcPath)
        {
            free(realSrcPath);
        }

        if (realDestPath)
        {
            free(realDestPath);
        }

        // Do not copy over self when src is actually a link to dest
        if (samePath)
        {
            return false;
        }

        // Copy files
        fileNameList files = readDir(src, fileName::FILE, false, followLink);
        for (const fileName& item : files)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/item
                    << " to " << destFile/item << endl;
            }

            // File to file.
            Foam::cp(src/item, destFile/item, followLink);
        }

        // Copy sub directories.
        fileNameList dirs = readDir
        (
            src,
            fileName::DIRECTORY,
            false,
            followLink
        );

        for (const fileName& item : dirs)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/item
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            Foam::cp(src/item, destFile, followLink);
        }
    }
    else
    {
        return false;
    }

    return true;
}


bool Foam::ln(const fileName& src, const fileName& dst)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : Create symlink from : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (src.empty())
    {
        WarningInFunction
            << "source name is empty: not linking." << endl;
        return false;
    }

    if (dst.empty())
    {
        WarningInFunction
            << "destination name is empty: not linking." << endl;
        return false;
    }

    if (exists(dst))
    {
        WarningInFunction
            << "destination " << dst << " already exists. Not linking."
            << endl;
        return false;
    }

    if (src.isAbsolute() && !exists(src))
    {
        WarningInFunction
            << "source " << src << " does not exist." << endl;
        return false;
    }

    if (::symlink(src.c_str(), dst.c_str()) == 0)
    {
        return true;
    }

    WarningInFunction
        << "symlink from " << src << " to " << dst << " failed." << endl;
    return false;
}


Foam::fileName Foam::readLink(const fileName& link)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : Returning symlink destination for : " << link << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (link.empty())
    {
        // Treat an empty path as a no-op.
        return fileName();
    }

    fileName result;
    result.resize(1024);  // Should be large enough (mostly relative anyhow)

    ssize_t len = ::readlink(link.c_str(), &(result.front()), result.size());
    if (len > 0)
    {
        result.resize(len);
        return result;
    }

    // Failure: return empty result
    return fileName();
}


bool Foam::mv(const fileName& src, const fileName& dst, const bool followLink)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Move : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore empty names => always false
    if (src.empty() || dst.empty())
    {
        return false;
    }

    if
    (
        dst.type() == fileName::DIRECTORY
     && src.type(followLink) != fileName::DIRECTORY
    )
    {
        const fileName dstName(dst/src.name());

        return (0 == std::rename(src.c_str(), dstName.c_str()));
    }

    return (0 == std::rename(src.c_str(), dst.c_str()));
}


bool Foam::mvBak(const fileName& src, const std::string& ext)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : moving : " << src << " to extension " << ext << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name or extension => always false
    if (src.empty() || ext.empty())
    {
        return false;
    }

    if (exists(src, false))
    {
        constexpr const int maxIndex = 99;
        char index[4];

        for (int n = 0; n <= maxIndex; ++n)
        {
            fileName dstName(src + "." + ext);
            if (n)
            {
                ::snprintf(index, 4, "%02d", n);
                dstName += index;
            }

            // avoid overwriting existing files, except for the last
            // possible index where we have no choice
            if (!exists(dstName, false) || n == maxIndex)
            {
                return (0 == std::rename(src.c_str(), dstName.c_str()));
            }
        }
    }

    // fallthrough: nothing to do
    return false;
}


bool Foam::rm(const fileName& file)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Removing : " << file << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    if (file.empty())
    {
        return false;
    }

    // If removal of plain file name fails, try with .gz

    return
    (
        0 == ::remove(file.c_str())
     || 0 == ::remove((file + ".gz").c_str())
    );
}


bool Foam::rmDir
(
    const fileName& directory,
    const bool silent,
    const bool emptyOnly
)
{
    if (directory.empty())
    {
        return false;
    }

    // Iterate contents (ignores an empty directory name)
    // Also retain hidden files/dirs for removal

    POSIX::directoryIterator dirIter(directory, true);
    if (!dirIter.exists())
    {
        if (!silent && !emptyOnly)
        {
            WarningInFunction
                << "Cannot open directory " << directory << endl;
        }

        return false;
    }

    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : removing directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Process each directory entry, counting any errors encountered
    int nErrors = 0;

    for (/*nil*/; dirIter; ++dirIter)
    {
        const std::string& item = *dirIter;

        // Allow invalid characters (spaces, quotes, etc),
        // otherwise we cannot remove subdirs with these types of names.
        // -> const fileName path = directory/name; <-

        const fileName path(fileName::concat(directory, item));

        fileName::Type detected = path.type(false);  // No followLink

        if (detected == fileName::Type::DIRECTORY)
        {
            // Call silently for lower levels
            if (!rmDir(path, true, emptyOnly))
            {
                ++nErrors;
            }
        }
        else if (emptyOnly)
        {
            // Only remove empty directories (not files)
            ++nErrors;

            // Check for dead symlinks
            if (detected == fileName::Type::SYMLINK)
            {
                detected = path.type(true);  // followLink

                if (detected == fileName::Type::UNDEFINED)
                {
                    --nErrors;

                    if (!rm(path))
                    {
                        ++nErrors;
                    }
                }
            }
        }
        else
        {
            if (!rm(path))
            {
                ++nErrors;
            }
        }
    }

    if (nErrors == 0)
    {
        // No errors encountered - try to remove the top-level

        if (!rm(directory))
        {
            nErrors = -1;  // A top-level error
        }
    }

    if (nErrors && !silent && !emptyOnly)
    {
        WarningInFunction
            << "Failed to remove directory " << directory << endl;

        if (nErrors > 0)
        {
            Info<< "could not remove " << nErrors << " sub-entries" << endl;
        }
    }

    return (nErrors == 0);
}


unsigned int Foam::sleep(const unsigned int sec)
{
    return ::sleep(sec);
}


void Foam::fdClose(const int fd)
{
    if (close(fd) != 0)
    {
        FatalErrorInFunction
            << "close error on " << fd << endl
            << abort(FatalError);
    }
}


bool Foam::ping
(
    const std::string& destName,
    const label destPort,
    const label timeOut
)
{
    struct hostent *hostPtr;
    volatile int sockfd;
    struct sockaddr_in destAddr;      // will hold the destination addr
    u_int addr;

    if ((hostPtr = ::gethostbyname(destName.c_str())) == nullptr)
    {
        FatalErrorInFunction
            << "gethostbyname error " << h_errno << " for host " << destName
            << abort(FatalError);
    }

    // Get first of the SLL of addresses
    addr = (reinterpret_cast<struct in_addr*>(*(hostPtr->h_addr_list)))->s_addr;

    // Allocate socket
    sockfd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
    {
        FatalErrorInFunction
            << "socket error"
            << abort(FatalError);
    }

    // Fill sockaddr_in structure with dest address and port
    std::memset(reinterpret_cast<char *>(&destAddr), '\0', sizeof(destAddr));
    destAddr.sin_family = AF_INET;
    destAddr.sin_port = htons(ushort(destPort));
    destAddr.sin_addr.s_addr = addr;


    timer myTimer(timeOut);

    if (timedOut(myTimer))
    {
        // Setjmp from timer jumps back to here
        fdClose(sockfd);
        return false;
    }

    if
    (
        ::connect
        (
            sockfd,
            reinterpret_cast<struct sockaddr*>(&destAddr),
            sizeof(struct sockaddr)
        ) != 0
    )
    {
        // Connection refused. Check if network was actually used or not.

        int connectErr = errno;

        fdClose(sockfd);

        if (connectErr == ECONNREFUSED)
        {
            return true;
        }
        //perror("connect");

        return false;
    }

    fdClose(sockfd);

    return true;
}


bool Foam::ping(const std::string& host, const label timeOut)
{
    return ping(host, 222, timeOut) || ping(host, 22, timeOut);
}


namespace Foam
{
//! \cond fileScope
static int waitpid(const pid_t pid)
{
    // child status, return code from the exec etc.
    int status = 0;

    // in parent - blocking wait
    // modest treatment of signals (in child)
    // treat 'stopped' like exit (suspend/continue)

    while (true)
    {
        pid_t wpid = ::waitpid(pid, &status, WUNTRACED);

        if (wpid == -1)
        {
            FatalErrorInFunction
                << "some error occurred in child"
                << exit(FatalError);
            break;
        }

        if (WIFEXITED(status))
        {
            // child exited, get its return status
            return WEXITSTATUS(status);
        }

        if (WIFSIGNALED(status))
        {
            // child terminated by some signal
            return WTERMSIG(status);
        }

        if (WIFSTOPPED(status))
        {
            // child stopped by some signal
            return WSTOPSIG(status);
        }

        FatalErrorInFunction
            << "programming error, status from waitpid() not handled: "
            << status
            << exit(FatalError);
    }

    return -1;  // should not happen
}
//! \endcond
}


int Foam::system(const std::string& command, const bool bg)
{
    if (command.empty())
    {
        // Treat an empty command as a successful no-op.
        // From 'man sh' POSIX (man sh):
        //   "If the command_string operand is an empty string,
        //    sh shall exit with a zero exit status."
        return 0;
    }

    // TBD: vfork is deprecated as of macOS 12.0
    const pid_t child_pid = ::vfork();   // NB: vfork, not fork!

    if (child_pid == -1)
    {
        FatalErrorInFunction
            << "vfork() failed for system command " << command
            << exit(FatalError);

        return -1;  // fallback error value
    }
    else if (child_pid == 0)
    {
        // In child

        // Close or redirect file descriptors
        redirects(bg);

        // execl uses the current environ
        (void) ::execl
        (
            "/bin/sh",          // Path of the shell
            "sh",               // Command-name (name for the shell)
            "-c",               // Read commands from command_string operand
            command.c_str(),    // Command string
            reinterpret_cast<char*>(0)
        );

        // Obviously failed, since exec should not return
        FatalErrorInFunction
            << "exec failed: " << command
            << exit(FatalError);

        return -1;  // fallback error value
    }


    // In parent
    // - started as background process, or blocking wait for the child

    return (bg ? 0 : waitpid(child_pid));
}


int Foam::system(const CStringList& command, const bool bg)
{
    if (command.empty())
    {
        // Treat an empty command as a successful no-op.
        // For consistency with POSIX (man sh) behaviour for (sh -c command),
        // which is what is mostly being replicated here.
        return 0;
    }

    // NB: use vfork, not fork!
    // vfork behaves more like a thread and avoids copy-on-write problems
    // triggered by fork.
    // The normal system() command has a fork buried in it that causes
    // issues with infiniband and openmpi etc.

    // TBD: vfork is deprecated as of macOS 12.0
    const pid_t child_pid = ::vfork();

    if (child_pid == -1)
    {
        FatalErrorInFunction
            << "vfork() failed for system command " << command[0]
            << exit(FatalError);

        return -1;  // fallback error value
    }
    else if (child_pid == 0)
    {
        // In child

        // Close or redirect file descriptors
        redirects(bg);

        // execvp searches the path, uses the current environ
        (void) ::execvp(command[0], command.strings());

        // Obviously failed, since exec should not return
        FatalErrorInFunction
            << "exec(" << command[0] << ", ...) failed"
            << exit(FatalError);

        return -1;  // fallback error value
    }


    // In parent
    // - started as background process, or blocking wait for the child

    return (bg ? 0 : waitpid(child_pid));
}


int Foam::system(const Foam::UList<Foam::string>& command, const bool bg)
{
    if (command.empty())
    {
        // Treat an empty command as a successful no-op.
        return 0;
    }

    // Make a deep copy as C-strings
    const CStringList cmd(command);
    return Foam::system(cmd, bg);
}


void* Foam::dlOpen(const fileName& libName, const bool check)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlopen() of " << libName << std::endl;
    }

    void* handle = loadLibrary(libName);

    if (!handle)
    {
        fileName libso;

        if (!libName.has_path() && !libName.starts_with("lib"))
        {
            // Try with 'lib' prefix
            libso = "lib" + libName;
            handle = loadLibrary(libso);

            if (POSIX::debug)
            {
                std::cout
                    << "   dlopen() as " << libso << std::endl;
            }
        }
        else
        {
            libso = libName;
        }

        // With canonical library extension ("so" or "dylib"), which remaps
        // "libXX" to "libXX.so" as well as "libXX.so" -> "libXX.dylib"
        if (!handle && !libso.has_ext(EXT_SO))
        {
            libso.replace_ext(EXT_SO);
            handle = loadLibrary(libso);

            if (POSIX::debug)
            {
                std::cout
                    << "   dlopen() as " << libso << std::endl;
            }
        }
    }

    if (!handle && check)
    {
        WarningInFunction
            << "dlopen error : " << ::dlerror() << endl;
    }

    if (POSIX::debug)
    {
        std::cout
            << "dlopen() of " << libName
            << " handle " << handle << std::endl;
    }

    return handle;
}


void* Foam::dlOpen(const fileName& libName, std::string& errorMsg)
{
    // Call without emitting error message - we capture that ourselves
    void* handle = Foam::dlOpen(libName, false);

    if (!handle)
    {
        // Capture error message
        errorMsg = ::dlerror();
    }
    else
    {
        // No errors
        errorMsg.clear();
    }

    return handle;
}


Foam::label Foam::dlOpen
(
    std::initializer_list<fileName> libNames,
    const bool check
)
{
    label nLoaded = 0;

    for (const fileName& libName : libNames)
    {
        if (Foam::dlOpen(libName, check))
        {
            ++nLoaded;
        }
    }

    return nLoaded;
}


bool Foam::dlClose(void* handle)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlClose(void*)"
            << " : dlclose of handle " << handle << std::endl;
    }
    return ::dlclose(handle) == 0;
}


void* Foam::dlSymFind(void* handle, const std::string& symbol, bool required)
{
    if (!required && (!handle || symbol.empty()))
    {
        return nullptr;
    }

    if (POSIX::debug)
    {
        std::cout
            << "dlSymFind(void*, const std::string&, bool)"
            << " : dlsym of " << symbol << std::endl;
    }

    // Clear any old errors - see manpage dlopen
    (void) ::dlerror();

    // Get address of symbol
    void* fun = ::dlsym(handle, symbol.c_str());

    // Any error?
    char *err = ::dlerror();

    if (err)
    {
        if (!required)
        {
            return nullptr;
        }

        WarningInFunction
            << "Cannot lookup symbol " << symbol << " : " << err
            << endl;
    }

    return fun;
}


#ifndef __APPLE__
static int collectLibsCallback
(
    struct dl_phdr_info *info,
    size_t size,
    void *data
)
{
    Foam::DynamicList<Foam::fileName>* ptr =
        reinterpret_cast<Foam::DynamicList<Foam::fileName>*>(data);
    ptr->append(info->dlpi_name);
    return 0;
}
#endif


Foam::fileNameList Foam::dlLoaded()
{
    DynamicList<fileName> libs;
    #ifdef __APPLE__
    for (uint32_t i=0; i < _dyld_image_count(); ++i)
    {
       libs.append(_dyld_get_image_name(i));
    }
    #else
    dl_iterate_phdr(collectLibsCallback, &libs);
    #endif

    if (POSIX::debug)
    {
        std::cout
            << "dlLoaded()"
            << " : determined loaded libraries :" << libs.size() << std::endl;
    }
    return libs;
}


// ************************************************************************* //
