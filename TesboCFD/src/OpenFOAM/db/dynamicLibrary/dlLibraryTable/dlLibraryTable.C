/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "dlLibraryTable.H"
#include "OSspecific.H"
#include "IOstreams.H"

// Could be constexpr in the header if required
#ifdef __APPLE__
    #define EXT_SO  "dylib"
#elif defined _WIN32
    #define EXT_SO  "dll"
#else
    #define EXT_SO  "so"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dlLibraryTable, 0);
}

int Foam::dlLibraryTable::dlcloseOnTerminate
(
    Foam::debug::optimisationSwitch("dlcloseOnTerminate", 0)
);


std::unique_ptr<Foam::dlLibraryTable> Foam::dlLibraryTable::global_(nullptr);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::dlLibraryTable::basename(const fileName& libPath)
{
    word libName(libPath.stem());
    libName.removeStart("lib");  // Remove leading 'lib' from name
    return libName;
}


Foam::word Foam::dlLibraryTable::fullname(word libName)
{
    if (libName.empty())
    {
        return libName;
    }

    // Add leading 'lib' and trailing '.so'
    return "lib" + libName.ext(EXT_SO);
}


Foam::dlLibraryTable& Foam::dlLibraryTable::libs()
{
    if (!global_)
    {
        global_.reset(new dlLibraryTable{});
    }

    return *global_;
}


bool Foam::dlLibraryTable::functionHook
(
    const bool load,
    void* handle,
    const std::string& funcName,
    const bool verbose,
    const std::string& context
)
{
    if (!handle || funcName.empty())
    {
        return false;
    }

    bool ok = false;

    void* symbol = Foam::dlSymFind(handle, funcName);

    if (symbol)
    {
        // Execute loader/unloader code
        try
        {
            loaderType fun = reinterpret_cast<loaderType>(symbol);

            if (fun)
            {
                (*fun)(load);
                ok = true;
            }
        }
        catch (...)
        {}
    }

    if (verbose && !ok)
    {
        auto& err = WarningInFunction
            << "Failed symbol lookup " << funcName.c_str() << nl;

        if (!context.empty())
        {
            err << "from " << context.c_str() << nl;
        }
    }

    return ok;
}


bool Foam::dlLibraryTable::loadHook
(
    void* handle,
    const std::string& funcName,
    const bool verbose,
    const std::string& context
)
{
    return functionHook(true, handle, funcName, verbose, context);
}


bool Foam::dlLibraryTable::unloadHook
(
    void* handle,
    const std::string& funcName,
    const bool verbose,
    const std::string& context
)
{
    return functionHook(false, handle, funcName, verbose, context);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void* Foam::dlLibraryTable::openLibrary
(
    const fileName& libName,
    bool verbose
)
{
    if (libName.empty())
    {
        return nullptr;
    }

    std::string msg;
    void* ptr = Foam::dlOpen(fileName(libName).expand(), msg);

    DebugInFunction
        << "Opened " << libName
        << " resulting in handle " << Foam::name(ptr) << nl;

    if (!ptr)
    {
        // Even with details turned off, we want some feedback about failure
        OSstream& os = (verbose > 0 ? WarningInFunction : InfoErr);
        os << "Could not load " << libName << nl << msg.c_str() << endl;
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable
(
    const UList<fileName>& libNames,
    bool verbose
)
{
    dlLibraryTable::open(libNames, verbose);
}


Foam::dlLibraryTable::dlLibraryTable
(
    std::initializer_list<fileName> libNames,
    bool verbose
)
{
    dlLibraryTable::open(libNames, verbose);
}


Foam::dlLibraryTable::dlLibraryTable
(
    const word& libsEntry,
    const dictionary& dict,
    bool verbose
)
{
    List<fileName> libNames;
    dict.readIfPresent(libsEntry, libNames, keyType::LITERAL);
    dlLibraryTable::open(libNames, verbose);
}


Foam::dlLibraryTable::dlLibraryTable
(
    const dictionary& dict,
    const word& libsEntry,
    bool verbose

)
:
    dlLibraryTable(libsEntry, dict, verbose)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    if (dlLibraryTable::dlcloseOnTerminate)
    {
        close();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::empty() const noexcept
{
    for (const void* ptr : libPtrs_)
    {
        if (ptr != nullptr)
        {
            return false;
        }
    }

    return true;
}


Foam::label Foam::dlLibraryTable::size() const noexcept
{
    label nLoaded = 0;

    for (const void* ptr : libPtrs_)
    {
        if (ptr != nullptr)
        {
            ++nLoaded;
        }
    }

    return nLoaded;
}


void Foam::dlLibraryTable::clear()
{
    libPtrs_.clear();
    libNames_.clear();
}


Foam::List<Foam::fileName> Foam::dlLibraryTable::loaded() const
{
    List<fileName> list(libNames_.size());

    label nLoaded = 0;

    forAll(libNames_, i)
    {
        void* ptr = libPtrs_[i];
        const fileName& libName = libNames_[i];

        if (ptr != nullptr && !libName.empty())
        {
            list[nLoaded] = libName;
            ++nLoaded;
        }
    }

    list.resize(nLoaded);

    return list;
}


void Foam::dlLibraryTable::close(bool verbose)
{
    label nLoaded = 0;

    forAllReverse(libPtrs_, i)
    {
        void* ptr = libPtrs_[i];

        if (ptr == nullptr)
        {
            libNames_[i].clear();
            continue;
        }

        if (Foam::dlClose(ptr))
        {
            DebugInFunction
                << "Closed [" << i << "] " << libNames_[i]
                << " with handle " << Foam::name(ptr) << nl;

            libPtrs_[i] = nullptr;
            libNames_[i].clear();
        }
        else
        {
            ++nLoaded;  // Still loaded

            if (verbose)
            {
                WarningInFunction
                    << "Failed closing " << libNames_[i]
                    << " with handle " << Foam::name(ptr) << endl;
            }
        }
    }


    // Compact the lists
    if (nLoaded && nLoaded != libPtrs_.size())
    {
        nLoaded = 0;

        forAll(libPtrs_, i)
        {
            if (libPtrs_[i] != nullptr)
            {
                if (nLoaded != i)
                {
                    libPtrs_[nLoaded] = libPtrs_[i];
                    libNames_[nLoaded] = std::move(libNames_[i]);
                }

                ++nLoaded;
            }
        }
    }

    libPtrs_.resize(nLoaded);
    libNames_.resize(nLoaded);
}


bool Foam::dlLibraryTable::push_back(const fileName& libName)
{
    if (libName.empty() || libNames_.contains(libName))
    {
        return false;
    }

    libPtrs_.push_back(nullptr);
    libNames_.push_back(libName);

    return true;
}


Foam::label Foam::dlLibraryTable::push_back(const UList<fileName>& libNames)
{
    label nAdded = 0;

    for (const fileName& libName : libNames)
    {
        if (push_back(libName))
        {
            ++nAdded;
        }
    }

    return nAdded;
}


bool Foam::dlLibraryTable::open(bool verbose)
{
    label nOpen = 0;
    label nCand = 0;  // Number of candidates (have libName but no pointer)

    forAll(libPtrs_, i)
    {
        void* ptr = libPtrs_[i];
        const fileName& libName = libNames_[i];

        if (ptr == nullptr && !libName.empty())
        {
            ++nCand;

            ptr = openLibrary(libName, verbose);

            if (ptr)
            {
                ++nOpen;
                libPtrs_[i] = ptr;
            }
            else
            {
                libNames_[i].clear();  // Avoid trying again
            }
        }
    }

    return nOpen && nOpen == nCand;
}


void* Foam::dlLibraryTable::open
(
    const fileName& libName,
    bool verbose
)
{
    // Handles empty name silently
    void* ptr = openLibrary(libName, verbose);

    if (ptr)
    {
        libPtrs_.push_back(ptr);
        libNames_.push_back(libName);
    }

    return ptr;
}


bool Foam::dlLibraryTable::open
(
    const UList<fileName>& libNames,
    bool verbose
)
{
    decltype(libNames.size()) nOpen = 0;

    for (const fileName& libName : libNames)
    {
        const label index = libNames_.find(libName);

        if (index >= 0 && libPtrs_[index] != nullptr)
        {
            // Already known and opened
            ++nOpen;
        }
        else if (dlLibraryTable::open(libName, verbose))
        {
            ++nOpen;
        }
    }

    return nOpen && nOpen == libNames.size();
}


bool Foam::dlLibraryTable::open
(
    std::initializer_list<fileName> libNames,
    bool verbose
)
{
    decltype(libNames.size()) nOpen = 0;

    for (const fileName& libName : libNames)
    {
        const label index = libNames_.find(libName);

        if (index >= 0 && libPtrs_[index] != nullptr)
        {
            // Already known and opened
            ++nOpen;
        }
        else if (dlLibraryTable::open(libName, verbose))
        {
            ++nOpen;
        }
    }

    return nOpen && nOpen == libNames.size();
}


bool Foam::dlLibraryTable::close
(
    const fileName& libName,
    bool verbose
)
{
    const label index = libNames_.rfind(libName);

    if (index < 0 || libName.empty())
    {
        return false;
    }

    void* ptr = libPtrs_[index];

    if (ptr == nullptr)
    {
        libNames_[index].clear();
        return false;
    }

    DebugInFunction
        << "Closing " << libName
        << " with handle " << Foam::name(ptr) << nl;

    const bool ok = Foam::dlClose(ptr);

    libPtrs_[index] = nullptr;
    libNames_[index].clear();

    if (ok)
    {
        // From man dlopen(3)
        // ...
        // a dynamically loaded shared object is not deallocated until
        // dlclose() has been called on it as many times as dlopen()
        // has succeeded on it.

        // Handle aliased library names
        for (label idx = 0; (idx = libPtrs_.find(ptr, idx)) >= 0; ++idx)
        {
            (void) Foam::dlClose(ptr);
            libPtrs_[idx] = nullptr;
            libNames_[idx].clear();
        }
    }
    else if (verbose)
    {
        WarningInFunction
            << "Could not close " << libName << endl;
    }

    return ok;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& libName)
{
    const label index = libNames_.rfind(libName);

    if (index < 0 || libName.empty())
    {
        return nullptr;
    }

    return libPtrs_[index];
}


bool Foam::dlLibraryTable::open
(
    const word& libsEntry,
    const dictionary& dict,
    bool verbose
)
{
    List<fileName> libNames;
    return
    (
        dict.readIfPresent(libsEntry, libNames, keyType::LITERAL)
     && dlLibraryTable::open(libNames, verbose)
    );
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    return dlLibraryTable::open(libsEntry, dict, true); // verbose = true
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<dlLibraryTable>& iproxy
)
{
    const auto& tbl = *iproxy;

    os << token::BEGIN_LIST << nl;

    // Lengths of pointers/names are guaranteed internally to be identical
    forAll(tbl.pointers(), i)
    {
        const void* ptr = tbl.pointers()[i];
        const fileName& libName = tbl.names()[i];

        // Also write out empty filenames
        // (specified with '-lib' but did not load)

        os  << Foam::name(ptr) << token::SPACE << libName << nl;
    }

    os << token::END_LIST << nl;

    return os;
}


// ************************************************************************* //
