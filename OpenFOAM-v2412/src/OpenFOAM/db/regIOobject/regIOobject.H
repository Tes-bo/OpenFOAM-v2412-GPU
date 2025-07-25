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

Class
    Foam::regIOobject

Description
    regIOobject is an abstract class derived from IOobject to handle
    automatic object registration with the objectRegistry.

SourceFiles
    regIOobject.C
    regIOobjectI.H
    regIOobjectRead.C
    regIOobjectWrite.C
    regIOobjectMetaData.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regIOobject_H
#define Foam_regIOobject_H

#include "IOobject.H"
#include "OSspecific.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// Forward Declarations
class dictionary;
class regIOobject;

namespace functionEntries
{
    class codeStream;
}
namespace fileOperations
{
    class uncollatedFileOperation;
}

/*---------------------------------------------------------------------------*\
                         Class regIOobject Declaration
\*---------------------------------------------------------------------------*/

class regIOobject
:
    public IOobject
{
protected:

        //- Helper: check readOpt flags and read if necessary
        bool readHeaderOk
        (
            const IOstreamOption::streamFormat fmt,
            const word& typeName
        );

        //- To flag master-only reading of objects
        static bool masterOnlyReading;


private:

    // Private Data

        //- Is this object registered with the registry
        bool registered_;

        //- Is this object owned by the registry
        bool ownedByRegistry_;

        //- The eventNo of last update
        label eventNo_;

        //- List of additional files to watch
        mutable DynamicList<fileName> watchFiles_;

        //- List of modification watch indices
        mutable DynamicList<label> watchIndices_;

        //- Dictionary for any meta-data
        autoPtr<dictionary> metaDataPtr_;

        //- Istream for reading
        autoPtr<ISstream> isPtr_;


    // Private Member Functions

        //- Construct object stream, read header if not already constructed
        void readStream(const bool readOnProc);


public:

        //- Friendship with classes needing access to masterOnlyReading
        friend class functionEntries::codeStream;
        friend class fileOperations::uncollatedFileOperation;


    // Generated Methods

        //- No copy assignment
        void operator=(const regIOobject&) = delete;


    // Static Data

        //- Runtime type information
        TypeName("regIOobject");


    // Constructors

        //- Construct from IOobject. The optional flag adds special handling
        //- if the object is the top-level regIOobject (eg, Time).
        regIOobject(const IOobject& io, const bool isTimeObject = false);

        //- Copy construct
        regIOobject(const regIOobject& rio);

        //- Copy construct, transferring registry registration to the copy
        //- if registerCopy is true
        regIOobject(const regIOobject& rio, bool registerCopy);

        //- Copy construct with new name, transferring registry registration
        //- to the copy if registerCopy is true
        regIOobject(const word& newName, const regIOobject&, bool registerCopy);

        //- Copy construct with new IO parameters
        regIOobject(const IOobject& io, const regIOobject& rio);


    //- Destructor
    virtual ~regIOobject();


    // Member Functions

        // Registration

            //- Add object to registry, if not already registered
            //  \return true if object was already registered,
            //      or was newly registered
            bool checkIn();

            //- Remove object from registry, and remove all file watches
            //  \return true if object was registered and was removed
            bool checkOut();

            //- Add file watch on object (if registered and READ_IF_MODIFIED)
            virtual void addWatch();

            //- Query the registered state (ie, has been checked in).
            //- This is not necessarily the same as registerObject(),
            //- which is just a stated preference.
            inline bool registered() const noexcept;

            //- Is this object owned by the registry?
            inline bool ownedByRegistry() const noexcept;

            //- Register object with its registry
            //- and transfer ownership to the registry.
            //  \return true if now ownedByRegistry
            inline bool store();

            //- Transfer pointer ownership to its registry.
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(Type* p);

            //- Transfer pointer ownership to its registry.
            //  Resets (clears) the parameter.
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(std::unique_ptr<Type>&& ptr);

            //- Transfer pointer ownership to its registry.
            //  Resets (clears) the parameter.
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(autoPtr<Type>& ptr);

            //- Transfer pointer ownership to its registry.
            //  Resets (clears) the parameter.
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(autoPtr<Type>&& ptr);

            //- Transfer pointer ownership to its registry.
            //  Changes parameter from PTR to CREF (do not rely on this).
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(refPtr<Type>& ptr);

            //- Transfer pointer ownership to its registry.
            //  Changes parameter from PTR to CREF (do not rely on this).
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(refPtr<Type>&& ptr);

            //- Transfer pointer ownership to its registry.
            //  Changes parameter from PTR to CREF (do not rely on this).
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(tmp<Type>& ptr);

            //- Transfer pointer ownership to its registry.
            //  Changes parameter from PTR to CREF (do not rely on this).
            //  \return reference to the stored object
            template<class Type>
            inline static Type& store(tmp<Type>&& ptr);

            //- Set object as \b not ownedByRegistry
            //  \param unregister optionally set as non-registered too
            inline void release(const bool unregister = false) noexcept;


        // Dependency Checking

            //- Event number at last update.
            inline label eventNo() const noexcept;

            //- Event number at last update.
            inline label& eventNo() noexcept;

            //- Return true if up-to-date with respect to given object
            bool upToDate(const regIOobject&) const;

            //- Return true if up-to-date with respect to given objects
            bool upToDate
            (
                const regIOobject&,
                const regIOobject&
            ) const;

            //- Return true if up-to-date with respect to given objects
            bool upToDate
            (
                const regIOobject&,
                const regIOobject&,
                const regIOobject&
            ) const;

            //- Return true if up-to-date with respect to given objects
            bool upToDate
            (
                const regIOobject&,
                const regIOobject&,
                const regIOobject&,
                const regIOobject&
            ) const;


            //- Set as up-to-date
            void setUpToDate();


        // Edit

            //- Rename
            virtual void rename(const word& newName);


        // Meta-data

            //- Return pointer to meta-data or nullptr
            const dictionary* findMetaData() const noexcept;

            //- Get or create meta-data
            dictionary& getMetaData() noexcept;

            //- Remove meta-data
            void removeMetaData();

            //- Update internal meta-data (eg, prior to writing)
            virtual void updateMetaData();


        // Reading

            //- Return complete path + object name if the file exists
            //  in the case directory otherwise null. Does not search
            //  up if parallel. Can be overridden to provide this functionality
            //  (e.g. IOdictionary)
            virtual fileName filePath() const;

            //- Read and check header info. Does not check the headerClassName.
            bool headerOk();

            //- Return Istream and check object type against that given
            Istream& readStream(const word&, const bool readOnProc = true);

            //- Close Istream
            void close();

            //- Virtual readData function.
            //  Must be defined in derived types for which
            //  re-reading is required
            virtual bool readData(Istream&);

            //- Read object
            virtual bool read();

            //- Add file watch for fileName on object if not yet watched.
            //  \return index of watch
            virtual label addWatch(const fileName&);

            //- Read access to file-monitoring handles
            inline const labelList& watchIndices() const noexcept;

            //- Write access to file-monitoring handles
            inline labelList& watchIndices() noexcept;

            //- Return true if the object's file (or files for objectRegistry)
            //- have been modified. (modified state is cached by Time)
            virtual bool modified() const;

            //- Read object if modified (as set by call to modified)
            virtual bool readIfModified();


        // Writing

            //- Pure virtual writeData function.
            //  Must be defined in derived types
            virtual bool writeData(Ostream&) const = 0;

            //- Write using stream options
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool writeOnProc
            ) const;

            //- Write using setting from DB
            virtual bool write(const bool writeOnProc = true) const;


        // Other

            //- Is object global
            virtual bool global() const
            {
                return false;
            }


    // Member Operators

        //- Copy assignment
        FOAM_DEPRECATED_STRICT(2023-12, "possibly remove in the future")
        void operator=(const IOobject& io);


    // Housekeeping

        //- Write using given format, version and compression
        FOAM_DEPRECATED_FOR(2020-02, "writeObject(IOstreamOption, bool)")
        virtual bool writeObject
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver,
            IOstreamOption::compressionType cmp,
            const bool writeOnProc
        ) const
        {
            return writeObject(IOstreamOption(fmt, ver, cmp), writeOnProc);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regIOobjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
