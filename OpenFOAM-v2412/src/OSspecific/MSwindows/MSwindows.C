/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2011 Symscape
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
    MS-Windows versions of the functions declared in OSspecific.H

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "MSwindows.H"
#include "fileName.H"
#include "fileStat.H"
#include "DynamicList.H"
#include "CStringList.H"
#include "IOstreams.H"
#include "Pstream.H"
#undef DebugInfo    // Windows name clash with OpenFOAM messageStream

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <unordered_map>

// Windows headers
#define WIN32_LEAN_AND_MEAN
#include <csignal>
#include <io.h>     // For _close
#include <windows.h>

#define EXT_SO  "dll"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MSwindows, 0);
}


namespace Foam
{
    // Don't abort under windows, causes abort dialog to popup.
    // Instead just exit with exitCode.
    static void sigAbortHandler(int exitCode)
    {
        ::exit(exitCode);
    }

    static bool installAbortHandler()
    {
        // If it didn't succeed there's not much we can do,
        // so don't check result.
        ::signal(SIGABRT, &sigAbortHandler);
        return true;
    }

    static bool const abortHandlerInstalled = installAbortHandler();


    // Move file, overwriting existing
    static bool renameFile(const std::string& src, const std::string& dst)
    {
        constexpr const int flags
        (
            MOVEFILE_COPY_ALLOWED
          | MOVEFILE_REPLACE_EXISTING
          | MOVEFILE_WRITE_THROUGH
        );

        // TODO: handle extra-long paths with ::MoveFileExW

        return ::MoveFileExA(src.c_str(), dst.c_str(), flags);
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * //

namespace Foam
{
namespace MSwindows
{

//- A simple directory contents iterator
class directoryIterator
{
    HANDLE handle_;

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

    //- Construct for dirName, optionally allowing hidden files/dirs
    directoryIterator(const std::string& dirName, bool allowHidden = false)
    :
        handle_(INVALID_HANDLE_VALUE),
        exists_(false),
        hidden_(allowHidden),
        item_()
    {
        if (!dirName.empty())
        {
            WIN32_FIND_DATA findData;

            handle_ = ::FindFirstFile((dirName/"*").c_str(), &findData);

            if (good())
            {
                exists_ = true;
                item_ = findData.cFileName;

                // If first element is not acceptable, get another one
                if (!accept())
                {
                    next();
                }
            }
        }
    }


    //- Destructor
    ~directoryIterator()
    {
        close();
    }


    // Member Functions

        //- Directory existed for opening
        bool exists() const noexcept
        {
            return exists_;
        }

        //- Directory pointer is valid
        bool good() const noexcept
        {
            return (INVALID_HANDLE_VALUE != handle_);
        }

        //- Close directory
        void close()
        {
            if (good())
            {
                ::FindClose(handle_);
                handle_ = INVALID_HANDLE_VALUE;
            }
        }

        //- The current item
        const std::string& val() const noexcept
        {
            return item_;
        }

        //- Read next item, always ignoring "." and ".." entries.
        //  Normally also ignore hidden files/dirs (beginning with '.')
        //  Automatically close when it runs out of items
        bool next()
        {
            if (good())
            {
                WIN32_FIND_DATA findData;

                while (::FindNextFile(handle_, &findData))
                {
                    item_ = findData.cFileName;

                    if (accept())
                    {
                        return true;
                    }
                }
                close();  // No more items
            }

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

} // End namespace MSwindows
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::string Foam::MSwindows::lastError()
{
    // Retrieve the system error message for the last-error code

    // Based on an example at:
    // http://msdn2.microsoft.com/en-us/library/ms680582(VS.85).aspx

    LPVOID lpMsgBuf;
    DWORD dw = GetLastError();

    FormatMessage
    (
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,  // source
        dw,    // message-id
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),  // language-id
        reinterpret_cast<LPTSTR>(&lpMsgBuf),
        0, NULL
    );

    const char* fmt = "Error %d: %s";

    // Use snprintf with zero to establish the size (without '\0') required
    std::string output;

    int n = ::snprintf(nullptr, 0, fmt, static_cast<LPCTSTR>(lpMsgBuf));

    if (n > 0)
    {
        output.resize(n+1);

        // Print directly into buffer
        n = ::snprintf(&(output[0]), n+1, fmt, static_cast<LPCTSTR>(lpMsgBuf));
        output.resize(n);
    }

    LocalFree(lpMsgBuf);

    return output;
}


std::string Foam::MSwindows::userName()
{
    const DWORD bufLen = 256;
    TCHAR buf[bufLen];
    DWORD len = bufLen;

    if (::GetUserName(buf, &len))
    {
        return buf;
    }

    std::string str;

    if
    (
        ERROR_INSUFFICIENT_BUFFER == ::GetLastError()
     && len < 2048
    )
    {
        // The len is with trailing '\0'
        str.resize(len);

        // Retrieve directly into buffer
        ::GetUserName(&(str[0]), &len);

        // Without trailing '\0'
        str.resize(len-1);
    }

    return str;
}


pid_t Foam::pid()
{
    const DWORD processId = ::GetCurrentProcessId();
    return processId;
}


pid_t Foam::ppid()
{
    // No equivalent under windows.

    if (MSwindows::debug)
    {
        Info<< "ppid not supported under MSwindows" << endl;
    }

    return 0;
}


pid_t Foam::pgid()
{
    // No equivalent under windows.

    if (MSwindows::debug)
    {
        Info<< "pgid not supported under MSwindows" << endl;
    }

    return 0;
}


bool Foam::hasEnv(const std::string& envName)
{
    // An empty envName => always false
    return !envName.empty() &&
        ::GetEnvironmentVariable(envName.c_str(), nullptr, 0);
}


Foam::string Foam::getEnv(const std::string& envName)
{
    std::string env;

    auto len = ::GetEnvironmentVariable(envName.c_str(), nullptr, 0);

    // len [return] = size with trailing nul char, or zero on failure
    if (len)
    {
        env.resize(len);

        // len [in] = size with trailing nul char
        // len [return] = size without trailing nul char
        len = ::GetEnvironmentVariable(envName.c_str(), &(env[0]), len);

        env.resize(len);
        return fileName::validate(env);
    }

    return env;
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
     && ::SetEnvironmentVariable(envName.c_str(), value.c_str())
    );
}


Foam::string Foam::hostName()
{
    const DWORD bufLen = MAX_COMPUTERNAME_LENGTH + 1;
    TCHAR buf[bufLen];
    DWORD len = bufLen;

    return ::GetComputerName(buf, &len) ? buf : string();
}


// DEPRECATED (2022-01)
Foam::string Foam::hostName(bool)
{
    return Foam::hostName();
}


// DEPRECATED (2022-01)
Foam::string Foam::domainName()
{
    // Could use ::gethostname and ::gethostbyname like POSIX.C, but would
    // then need to link against ws_32. Prefer to minimize dependencies.

    return string();
}


Foam::string Foam::userName()
{
    string name = Foam::getEnv("USERNAME");

    if (name.empty())
    {
        name = MSwindows::userName();
    }

    return name;
}


bool Foam::isAdministrator()
{
    // Assume worst case
    return true;
}


Foam::fileName Foam::home()
{
    fileName env = Foam::getEnv("HOME");

    if (env.empty())
    {
        env = Foam::getEnv("USERPROFILE");
    }

    return env;
}


Foam::fileName Foam::home(const std::string& userName)
{
    return Foam::home();
}


Foam::fileName Foam::cwd()
{
    string path;
    auto len = ::GetCurrentDirectory(0, nullptr);

    // len [return] = size with trailing nul char, or zero on failure
    if (len)
    {
        path.resize(len);

        // len [in] = size with trailing nul char
        // len [return] = size without trailing nul char
        len = ::GetCurrentDirectory(len, &(path[0]));

        path.resize(len);
        return fileName::validate(path);
    }

    FatalErrorInFunction
        << "Couldn't get the current working directory"
        << exit(FatalError);

    return fileName();
}


Foam::fileName Foam::cwd(bool logical)
{
    return Foam::cwd();
}


bool Foam::chDir(const fileName& dir)
{
    // Ignore an empty dir name => always false
    return !dir.empty() && ::SetCurrentDirectory(dir.c_str());;
}


bool Foam::mkDir(const fileName& pathName, const mode_t mode)
{
    // empty names are meaningless
    if (pathName.empty())
    {
        return false;
    }

    bool ok = ::CreateDirectory(pathName.c_str(), NULL);

    if (ok)
    {
        Foam::chMod(pathName, mode);
        return true;
    }

    const DWORD error = ::GetLastError();
    switch (error)
    {
        case ERROR_ALREADY_EXISTS:
        {
            ok = true;
            break;
        }

        case ERROR_PATH_NOT_FOUND:
        {
            // Part of the path does not exist so try to create it
            const fileName& parentName = pathName.path();

            if (parentName.size() && mkDir(parentName, mode))
            {
                ok = mkDir(pathName, mode);
            }
            break;
        }
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Couldn't create directory: " << pathName
            << " " << MSwindows::lastError()
            << exit(FatalError);
    }

    return ok;
}


bool Foam::chMod(const fileName& name, const mode_t m)
{
    // Ignore an empty name => always false
    return !name.empty() && _chmod(name.c_str(), m) == 0;
}


mode_t Foam::mode(const fileName& name, const bool followLink)
{
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


// Windows equivalent to S_ISDIR
#define ms_isdir(a) \
    ((m != INVALID_FILE_ATTRIBUTES) && (m & FILE_ATTRIBUTE_DIRECTORY))

// Windows equivalent to S_ISREG
#define ms_isreg(s) \
    ((m != INVALID_FILE_ATTRIBUTES) && !(m & FILE_ATTRIBUTE_DIRECTORY))


Foam::fileName::Type Foam::type
(
    const fileName& name,
    const bool /* followLink */
)
{
    // Ignore an empty name => always UNDEFINED
    if (name.empty())
    {
        return fileName::Type::UNDEFINED;
    }

    const DWORD m = ::GetFileAttributes(name.c_str());

    if (ms_isreg(m))
    {
        return fileName::Type::FILE;
    }
    else if (ms_isdir(m))
    {
        return fileName::Type::DIRECTORY;
    }

    return fileName::Type::UNDEFINED;
}


// Local check for gz file
static bool isGzFile(const std::string& name)
{
    const DWORD m = ::GetFileAttributes((name + ".gz").c_str());
    return ms_isreg(m);
}


bool Foam::exists
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    // Ignore an empty name => always false
    if (name.empty())
    {
        return false;
    }

    const DWORD m = ::GetFileAttributes(name.c_str());

    return (ms_isdir(m) || ms_isreg(m) || (checkGzip && isGzFile(name)));
}


bool Foam::isDir(const fileName& name, const bool followLink)
{
    // Ignore an empty name => always false
    if (name.empty())
    {
        return false;
    }

    const DWORD m = ::GetFileAttributes(name.c_str());

    return ms_isdir(m);
}


bool Foam::isFile
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    // Ignore an empty name => always false
    if (name.empty())
    {
        return false;
    }

    const DWORD m = ::GetFileAttributes(name.c_str());

    return (ms_isreg(m) || (!ms_isdir(m) && checkGzip && isGzFile(name)));
}


off_t Foam::fileSize(const fileName& name, const bool followLink)
{
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
    // Ignore an empty name
    return name.empty() ? 0 : fileStat(name, followLink).modTime();
}


double Foam::highResLastModified(const fileName& name, const bool followLink)
{
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
    // Initial filename list size
    // also used as increment if initial size found to be insufficient
    static constexpr int maxNnames = 100;

    fileNameList dirEntries;

    // Iterate contents (ignores an empty directory name)

    MSwindows::directoryIterator dirIter(directory);

    if (!dirIter.exists())
    {
        if (MSwindows::debug)
        {
            InfoInFunction
                << "cannot open directory " << directory << endl;
        }

        return dirEntries;
    }

    if (MSwindows::debug)
    {
        InfoInFunction
            << " : reading directory " << directory << endl;
    }

    label nFailed = 0;     // Entries with invalid characters
    label nEntries = 0;    // Number of selected entries
    dirEntries.resize(maxNnames);

    // Process all the directory entries
    for (/*nil*/; dirIter; ++dirIter)
    {
        const std::string& item = *dirIter;

        // Validate filename without quotes, etc in the name.
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
            fileName::Type detected = (directory/name).type();

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

    if (nFailed && MSwindows::debug)
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
    // Make sure source exists - this also handles an empty source name
    if (!exists(src))
    {
        return false;
    }

    fileName destFile(dest);

    const fileName::Type srcType = src.type(followLink);

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams.
        // - use binary mode to avoid any issues
        std::ifstream srcStream(src, ios_base::in | ios_base::binary);
        if (!srcStream)
        {
            return false;
        }

        // - use binary mode to avoid any issues
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

        // Make sure the destination directory extists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        // Copy files
        fileNameList files = readDir(src, fileName::FILE, false, followLink);
        for (const fileName& item : files)
        {
            if (MSwindows::debug)
            {
                Info<< "Copying : " << src/item
                    << " to " << destFile/item << endl;
            }

            // File to file.
            Foam::cp(src/item, destFile/item);
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
            if (MSwindows::debug)
            {
                Info<< "Copying : " << src/item
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            Foam::cp(src/item, destFile);
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
    // links are poorly supported, or need administrator privileges.
    // Skip for now.

    if (MSwindows::debug)
    {
        Info<< "MSwindows does not support softlinks" << endl;
    }

    return false;
}


Foam::fileName Foam::readLink(const fileName& link)
{
    if (MSwindows::debug)
    {
        Info<< "MSwindows does not support softlinks" << endl;
    }

    return fileName();
}


bool Foam::mv(const fileName& src, const fileName& dst, const bool followLink)
{
    if (MSwindows::debug)
    {
        Info<< "Move : " << src << " to " << dst << endl;
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

        return renameFile(src, dstName);
    }

    return renameFile(src, dst);
}


bool Foam::mvBak(const fileName& src, const std::string& ext)
{
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
                return renameFile(src, dstName);
            }
        }
    }

    // fall-through: nothing to do
    return false;
}


bool Foam::rm(const fileName& file)
{
    if (MSwindows::debug)
    {
        Info<< "Removing : " << file << endl;
    }

    // Ignore an empty name => always false
    if (file.empty())
    {
        return false;
    }


    // If removal of plain file name failed, try with .gz

    return
    (
        0 == std::remove(file.c_str())
     || 0 == std::remove((file + ".gz").c_str())
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

    MSwindows::directoryIterator dirIter(directory, true);
    if (!dirIter.exists())
    {
        if (!silent && !emptyOnly)
        {
            WarningInFunction
                << "Cannot open directory " << directory << endl;
        }

        return false;
    }

    if (MSwindows::debug)
    {
        InfoInFunction
            << " : removing directory " << directory << endl;
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

        if (detected == fileName::DIRECTORY)
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

            // POSIX has extra handling for dead symlinks...
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

        if (!::RemoveDirectory(directory.c_str()))
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
    ::Sleep(1000*sec);  // in milliseconds

    return 0;
}


void Foam::fdClose(const int fd)
{
    if (::_close(fd) != 0)
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
    // Appears that socket calls require administrator privileges.
    // Skip for now.

    if (MSwindows::debug)
    {
        Info<< "MSwindows does not support ping" << endl;
    }

    return false;
}


bool Foam::ping(const std::string& host, const label timeOut)
{
    return ping(host, 222, timeOut) || ping(host, 22, timeOut);
}


int Foam::system(const std::string& command, const bool bg)
{
    if (MSwindows::debug && bg)
    {
        InfoInFunction
            << "MSwindows does not support background (fork) tasks" << endl;
    }

    return std::system(command.c_str());
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

    const int count = command.size();

    std::string cmd;

    for (int i = 0; i < count; ++i)
    {
        if (i) cmd += ' ';
        cmd += command[i];
    }

    return system(cmd, bg);
}


int Foam::system(const UList<Foam::string>& command, const bool bg)
{
    if (command.empty())
    {
        // Treat an empty command as a successful no-op.
        return 0;
    }

    const int count = command.size();

    std::string cmd;

    for (int i = 0; i < count; ++i)
    {
        if (i) cmd += ' ';
        cmd += command[i];
    }

    return system(cmd, bg);
}


// Explicitly track loaded libraries, rather than use
// EnumerateLoadedModules64 and have to link against
// Dbghelp.dll
// Details at http://msdn.microsoft.com/en-us/library/ms679316(v=vs.85).aspx

static std::unordered_map<void*, std::string> libsLoaded;


void* Foam::dlOpen(const fileName& libName, const bool check)
{
    if (MSwindows::debug)
    {
        std::cout
            << "dlOpen(const fileName&)"
            << " : dlopen of " << libName << std::endl;
    }

    // Always remap "libXX.so" and "libXX" to "libXX.dll"
    fileName libso(libName.lessExt().ext(EXT_SO));

    void* handle = ::LoadLibrary(libso.c_str());

    if
    (
        !handle
     && !libName.has_path()
     && !libso.starts_with("lib")
    )
    {
        // Try with 'lib' prefix
        libso = "lib" + libso;
        handle = ::LoadLibrary(libso.c_str());

        if (MSwindows::debug)
        {
            std::cout
                << "dlOpen(const fileName&)"
                << " : dlopen of " << libso << std::endl;
        }
    }

    if (handle)
    {
        libsLoaded[handle] = libso.lessExt();
    }
    else if (check)
    {
        WarningInFunction
            << "dlopen error : " << MSwindows::lastError() << endl;
    }

    if (MSwindows::debug)
    {
        std::cout
            << "dlOpen(const fileName&)"
            << " : dlopen of " << libName
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
        errorMsg = MSwindows::lastError();
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


bool Foam::dlClose(void* const handle)
{
    if (MSwindows::debug)
    {
        std::cout
            << "dlClose(void*)"
            << " : dlclose of handle " << handle << std::endl;
    }

    const bool ok = ::FreeLibrary(static_cast<HMODULE>(handle));

    if (ok)
    {
        libsLoaded.erase(handle);
    }

    return ok;
}


void* Foam::dlSymFind(void* handle, const std::string& symbol, bool required)
{
    if (!required && (!handle || symbol.empty()))
    {
        return nullptr;
    }

    if (MSwindows::debug)
    {
        std::cout
            << "dlSymFind(void*, const std::string&, bool)"
            << " : dlsym of " << symbol << std::endl;
    }

    // Get address of symbol, or nullptr on failure
    void* fun =
        reinterpret_cast<void *>
        (
            ::GetProcAddress(static_cast<HMODULE>(handle), symbol.c_str())
        );

    // Any error?
    if (!fun && required)
    {
        WarningInFunction
            << "Cannot lookup symbol " << symbol << " : "
            << MSwindows::lastError() << endl;
    }

    return fun;
}


Foam::fileNameList Foam::dlLoaded()
{
    DynamicList<fileName> libs(libsLoaded.size());

    for (const auto& item : libsLoaded)
    {
        libs.append(item.second);
    }

    if (MSwindows::debug)
    {
        std::cout
            << "dlLoaded()"
            << " : determined loaded libraries :" << libs.size() << std::endl;
    }

    return libs;
}


// ************************************************************************* //
