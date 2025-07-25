/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries.

SeeAlso
    Foam::dlOpen
    Foam::dlClose

SourceFiles
    dlLibraryTable.C
    dlLibraryTableTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dlLibraryTable_H
#define Foam_dlLibraryTable_H

#include "fileName.H"
#include "DynamicList.H"
#include "InfoProxy.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dlLibraryTable;
Ostream& operator<<(Ostream&, const InfoProxy<dlLibraryTable>&);


/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Static Data

        //- Global singleton of dynamic libraries
        static std::unique_ptr<dlLibraryTable> global_;


    // Private Data

        //- Pointers to the loaded libraries
        DynamicList<void*> libPtrs_;

        //- Names of loaded libraries, or of libraries to be loaded
        DynamicList<fileName> libNames_;


    // Private Member Functions

        //- Load/unload hook.
        //  \return true if the function was found and executed
        static bool functionHook
        (
            const bool load,    //!< true = on-load, false = on-unload
            void* handle,       //!< library handle
            const std::string& funcName,
            const bool verbose,
            const std::string& context
        );


        //- Open specified library name and return pointer.
        //  Warning messages, but no additional side-effects.
        void* openLibrary(const fileName& libName, bool verbose);


public:

    // Static Data Members

        //- Use dlclose() when clearing the dlLibraryTable.
        //  This is presumably \em cleaner, but can also remove functions
        //  that are still needed (eg, to terminate the library itself)
        static int dlcloseOnTerminate;


    // Public Data Types

        //- Global loader/unloader function type (C-linkage)
        //  Called with true on load, false on unload.
        typedef void (*loaderType)(bool);


    //- Declare name of the class and its debug switch
    ClassName("dlLibraryTable");


    // Generated Methods

        //- Default construct
        dlLibraryTable() = default;

        //- No copy construct
        dlLibraryTable(const dlLibraryTable&) = delete;

        //- Move construct
        dlLibraryTable(dlLibraryTable&&) = default;

        //- No copy assignment
        void operator=(const dlLibraryTable&) = delete;

        //- Move assignment
        dlLibraryTable& operator=(dlLibraryTable&&) = default;


    // Constructors

        //- Open specified libraries, warn by default if problems occur
        //  Ignores duplicate names.
        explicit dlLibraryTable
        (
            const UList<fileName>& libNames,
            bool verbose = true
        );

        //- Open specified libraries, warn by default if problems occur
        //  Ignores duplicate names.
        explicit dlLibraryTable
        (
            std::initializer_list<fileName> libNames,
            bool verbose = true
        );

        //- Open libraries listed in 'libsEntry' entry in the dictionary,
        //- warn by default if problems occur
        dlLibraryTable
        (
            const word& libsEntry,
            const dictionary& dict,
            bool verbose = true
        );

        //- Open libraries listed in 'libsEntry' entry in the dictionary,
        //- warn by default if problems occur
        dlLibraryTable
        (
            const dictionary& dict,
            const word& libsEntry,
            bool verbose = true
        );


    //- Destructor. Closes all libraries loaded by the table.
    ~dlLibraryTable();


    // Static Member Functions

        //- Library basename without leading 'lib' or trailing '.so'
        static word basename(const fileName& libPath);

        //- Library fullname, prefix with 'lib', suffix with '.so'
        //  \note the suffix is system-dependent
        static word fullname(word libName);

        //- Table of global libraries
        static dlLibraryTable& libs();

        //- Low-level interface to execute global "void funcName(true)"
        //- from the library, typically for additional loading.
        //  If called, it should be the first step after opening a library.
        //  \return true if the function was found and executed
        static bool loadHook
        (
            void* handle,
            const std::string& funcName,
            const bool verbose = false,
            const std::string& context = "" //!< Calling context for warnings
        );

        //- Low-level interface to execute global "void funcName(false)"
        //- from the library, typically for unloading.
        //  If called, it should be the last step before closing a library.
        //  \return true if the function was found and executed
        static bool unloadHook
        (
            void* handle,
            const std::string& funcName,
            const bool verbose = false,
            const std::string& context = "" //!< Calling context for warnings
        );


    // Member Functions

        //- True if there are no libraries loaded by the table
        bool empty() const noexcept;

        //- The number of libraries loaded by the table
        label size() const noexcept;

        //- Names of the libraries in use
        List<fileName> loaded() const;

        //- Names of the libraries in use, or requested
        const UList<fileName>& names() const noexcept
        {
            return libNames_;
        }

        //- Pointers to the libraries in use. Access with caution.
        const UList<void*>& pointers() const noexcept
        {
            return libPtrs_;
        }

        //- Clears the table, without attempting to close the libraries
        void clear();


        //- Add to the list of names, but do not yet open.
        //  Ignores empty and duplicate names.
        bool push_back(const fileName& libName);

        //- Add to the list of names, but do not yet open.
        //  Ignores empty and duplicate names.
        label push_back(const UList<fileName>& libNames);


        //- Open named, but unopened libraries.
        //- These names will normally have been added with push_back().
        bool open(bool verbose = true);

        //- Open the named library, warn by default if problems occur.
        //  An empty name is a silent no-op and always returns nullptr.
        //  \return a pointer to the library opened, or nullptr on failure.
        void* open(const fileName& libName, bool verbose = true);

        //- Open the named libraries, warn by default if problems occur.
        //  Ignores duplicate names.
        bool open(const UList<fileName>& libNames, bool verbose = true);

        //- Open the named libraries, warn by default if problems occur.
        //  Ignores duplicate names.
        bool open
        (
            std::initializer_list<fileName> libNames,
            bool verbose = true
        );

        //- Close all libraries loaded by the table and remove the closed
        //- functions from the table.
        void close(bool verbose = true);

        //- Close the named library, optionally warn if problems occur
        //  Using an empty name is a no-op and always returns false.
        bool close(const fileName& libName, bool verbose = true);

        //- Find the handle of the named library
        //  Using an empty name is a no-op and always returns nullptr.
        void* findLibrary(const fileName& libName);

        //- Open libraries listed in the 'libsEntry' entry in the dictionary.
        bool open
        (
            const word& libsEntry,
            const dictionary& dict,
            bool verbose = true
        );

        //- Open libraries listed in the 'libsEntry' entry in the dictionary.
        //  Verbose = true
        bool open(const dictionary& dict, const word& libsEntry);

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary and check the additions
        //- to the given constructor table
        template<class TablePtr>
        bool open
        (
            const dictionary& dict,
            const word& libsEntry,
            const TablePtr& tablePtr,
            bool verbose = true
        );


    // Info

        //- Return info proxy,
        //- used to print library table information to a stream
        InfoProxy<dlLibraryTable> info() const noexcept { return *this; }


    // Housekeeping

        //- Add to the list of names, but do not yet open.
        //FOAM_DEPRECATED_FOR(2024-04, "push_back()")
        bool append(const fileName& libName)
        {
            return push_back(libName);
        }

        //- Add to the list of names, but do not yet open.
        //FOAM_DEPRECATED_FOR(2024-04, "push_back()")
        label append(const UList<fileName>& libNames)
        {
            return push_back(libNames);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
