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

Class
    Foam::IOobject

Description
    Defines the attributes of an object for which implicit
    objectRegistry management is supported, and provides the infrastructure
    for performing stream I/O.

    An IOobject is constructed with an object name, a class name, an instance
    path, a reference to a objectRegistry, and parameters determining its
    storage status.

    \par Read options

    Define what is done on object construction and explicit reads:
      - \par NO_READ
        Do not read
      - \par MUST_READ
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or cannot be read.
        Does not check timestamp or re-read.
      - \par READ_MODIFIED (MUST_READ_IF_MODIFIED)
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or cannot be read. If object is
        registered its timestamp will be checked every timestep and possibly
        re-read.
      - \par LAZY_READ (READ_IF_PRESENT)
        Read object from Istream, but only if Istream exists. \n
        Error only if Istream exists but cannot be read.
        Does not check timestamp or re-read.

    \par Write options

    Define what is done on object destruction and explicit writes:
      - \par NO_WRITE
        No automatic writing, but can be written explicitly
      - \par AUTO_WRITE
        Object is written automatically when requested to by the
        objectRegistry.

    When serializing, the IOobject characteristics are typically written
    as a \c FoamFile header, which is a sub-dictionary with the following
    type of content:

    \table
        Property | Description                            | Type | Reqd | Deflt
        version  | The base format version                | float | no  | 2.0
        format   | The stream format (ascii/binary)       | word  | yes |
        arch     | The architecture string                | string | no |
        note     | Descriptive note about the object      | string | no |
        location | The relative location of the object    | string | no |
        class    | The type of the object                 | word | yes |
        object   | The name of the object                 | word | yes |
    \endtable

Note
    Specifying registered does not result in the IOobject itself being
    registered. It is only serves as guidance for a regIOobject using it.

See also
    Foam::objectRegistry
    Foam::regIOobject

SourceFiles
    IOobject.C
    IOobjectReadHeader.C
    IOobjectWriteHeader.C
    IOobjectPrint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOobject_H
#define Foam_IOobject_H

#include "fileName.H"
#include "typeInfo.H"
#include "refPtr.H"  // For autoPtr, refPtr, tmp, stdFoam
#include "Enum.H"
#include "InfoProxy.H"
#include "IOobjectOption.H"
#include "IOstreamOption.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class dictionary;
class objectRegistry;
class IOobject;

template<>
Ostream& operator<<(Ostream&, const InfoProxy<IOobject>&);

// Traits

//- Trait for specifying global vs. local IOobject file types
template<class T>
struct is_globalIOobject : std::false_type {};


/*---------------------------------------------------------------------------*\
                           Class IOobject Declaration
\*---------------------------------------------------------------------------*/

class IOobject
:
    public IOobjectOption
{
public:

    // Public Data Types

        //- Enumeration defining the valid states of an IOobject
        enum objectState : char
        {
            GOOD,
            BAD
        };

        //- Enumeration defining the file checking options
        enum fileCheckTypes : char
        {
            timeStamp,
            timeStampMaster,
            inotify,
            inotifyMaster
        };

        //- Names for the fileCheckTypes
        static const Enum<fileCheckTypes> fileCheckTypesNames;


private:

    // Static Data Members

        //- Use an output file banner, enabled by default
        static bool bannerEnabled_;


    // Private Data (NB: byte-flags first for better alignment)

        //- The IOobject state
        objectState objState_;

        //- The sizeof (label) in bytes, possibly read from the header
        unsigned char sizeofLabel_;

        //- The sizeof (scalar) in bytes, possibly read from the header
        unsigned char sizeofScalar_;

        //- Name
        word name_;

        //- Class name read from header
        word headerClassName_;

        //- Optional note
        string note_;

        //- Instance path component
        fileName instance_;

        //- Local path component
        fileName local_;

        //- Reference to the objectRegistry
        const objectRegistry& db_;


    // Private Member Functions

        //- Construct from registry, io options. Without name, instance, local
        IOobject(const objectRegistry& registry, IOobjectOption ioOpt);

        //- Read header and check its info.
        //  Optionally checks headerClassName against the type-name.
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        bool readAndCheckHeader
        (
            const bool isGlobal,
            const word& typeName,
            const bool checkType = true,
            const bool search = true,
            const bool verbose = true
        );


protected:

    // Protected Member Functions

        //- Helper: write content for FoamFile IOobject header
        //- with optional meta information.
        static void writeHeaderContent
        (
            Ostream& os,
            const IOobject& io,
            const word& objectType,
            const dictionary* metaDataDict = nullptr
        );

        //- Helper: write dictionary content for FoamFile header
        //- with optional meta information.
        static void writeHeaderContent
        (
            dictionary& dict,
            const IOobject& io,
            const word& objectType,
            IOstreamOption streamOpt,
            const dictionary* metaDataDict = nullptr
        );

        //- Set the object state to bad
        void setBad(const string& s);


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("IOobject");


    // Static Data

        //- Character for scoping object names (':' or '_')
        //  Change with caution.
        static char scopeSeparator;

        //- Type of file modification checking
        static fileCheckTypes fileModificationChecking;

        //- Time skew (seconds) for file modification checks
        static float fileModificationSkew;

        //- Max number of times to poll for file modification changes
        static int maxFileModificationPolls;


    // Static Functions

        //- Status of output file banner
        static bool bannerEnabled() noexcept
        {
            return bannerEnabled_;
        }

        //- Enable/disable output file banner
        //  \return the previous value
        static bool bannerEnabled(bool on) noexcept
        {
            bool old(bannerEnabled_);
            bannerEnabled_ = on;
            return old;
        }

        //- Split path into instance, local, name components
        //
        // The splitting behaviour is as follows:
        // \verbatim
        //    input       | instance   | local | name
        //    ----------- | ---------- | ----- | ----
        //    a           |            |       | a
        //    a/b         | a          |       | b
        //    a/b/c/d     | a          | b/c   | d
        //    /a/b/c      | /a/b       |       | c
        //    ./a/b/c     | PWD/a/b    |       | c
        //    ../a/b/c    | PWD/../a/b |       | c
        //    a/b/        | ERROR      |       |
        // \endverbatim
        // where PWD is the Foam::cwd() current working directory
        static bool fileNameComponents
        (
            const fileName& path,
            fileName& instance,
            fileName& local,
            word& name
        );

        //- Create dot-delimited name.group string
        //  An empty group is ignored.
        template<class StringType>
        static inline word groupName(StringType base, const word& group);

        //- Return group (extension part of name)
        static word group(const word& name);

        //- Return member (name without the extension)
        static word member(const word& name);

        //- Create scope:name or scope_name string
        //  An empty scope is ignored.
        static inline word scopedName
        (
            const std::string& scope,
            const word& name
        );

        //- Create scope:name1:name2 or scope_name1_name2 string
        //  An empty scope is ignored.
        static inline word scopedName
        (
            const std::string& scope,
            const word& name1,
            const word& name2
        );

        //- Return the IOobject, but also consider an alternative file name.
        //
        //  \param io The expected IOobject to use
        //  \param altFile Alternative fileName (ignored if empty).
        //  \param ioName  The alternative name for the IOobject when
        //      the altFile resolves to a directory.
        //
        //  \note If the alternative fileName is a non-empty string,
        //      it defines the location but uses all other properties of the
        //      expected IOobject.
        //      The location may be an absolute or a relative path.
        //      If it corresponds to a directory, the name of the
        //      expected IOobject will be used in its resolution.
        //      This expected name can provided via the ioName parameter.
        static IOobject selectIO
        (
            const IOobject& io,
            const fileName& altFile,
            const word& ioName = ""
        );


    // Generated Methods

        //- Copy construct
        IOobject(const IOobject&) = default;

        //- Destructor
        virtual ~IOobject() = default;


    // Constructors

        //- Construct from name, instance, registry, io options
        //  (default: NO_READ, NO_WRITE, REGISTER, non-global)
        IOobject
        (
            const word& name,
            const fileName& instance,
            const objectRegistry& registry,
            IOobjectOption ioOpt = IOobjectOption()
        );

        //- Construct from name, instance, local, registry, io options
        //  (default: NO_READ, NO_WRITE, REGISTER, non-global)
        IOobject
        (
            const word& name,
            const fileName& instance,
            const fileName& local,
            const objectRegistry& registry,
            IOobjectOption ioOpt = IOobjectOption()
        );

        //- Construct from path, registry, io options.
        //  (default: NO_READ, NO_WRITE, REGISTER, non-global)
        //
        //  Uses fileNameComponents() to split path into components.
        //  A path that starts with a '/' is regarded as a file system path.
        //  Paths starting with either './' or '../' are relative to
        //  current working directory (and replaced with absolute equivalents).
        //  All other paths are considered to be relative to the case.
        IOobject
        (
            const fileName& path,
            const objectRegistry& registry,
            IOobjectOption ioOpt = IOobjectOption()
        );

        //- Construct from name, instance, registry, io options
        inline IOobject
        (
            const word& name,
            const fileName& instance,
            const objectRegistry& registry,
            IOobjectOption::readOption rOpt,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            bool registerObject = true,  // == IOobjectOption::LEGACY_REGISTER
            bool globalObject = false
        );

        //- Construct from name, instance, local, registry, io options
        inline IOobject
        (
            const word& name,
            const fileName& instance,
            const fileName& local,
            const objectRegistry& registry,
            IOobjectOption::readOption rOpt,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            bool registerObject = true,  // == IOobjectOption::LEGACY_REGISTER
            bool globalObject = false
        );

        //- Construct from path, registry, io options.
        //  Uses fileNameComponents() to split path into components.
        //  A path that starts with a '/' is regarded as a file system path.
        //  Paths starting with either './' or '../' are relative to
        //  current working directory (and replaced with absolute equivalents).
        //  All other paths are considered to be relative to the case.
        inline IOobject
        (
            const fileName& path,
            const objectRegistry& registry,
            IOobjectOption::readOption rOpt,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            bool registerObject = true,  // == IOobjectOption::LEGACY_REGISTER
            bool globalObject = false
        );

        //- Copy construct, resetting registry
        IOobject(const IOobject& io, const objectRegistry& registry);

        //- Copy construct, resetting name
        IOobject(const IOobject& io, const word& name);

        //- Copy construct, resetting name and local component
        inline IOobject
        (
            const IOobject& io,
            const word& name,
            const fileName& local
        );

        //- Copy construct, resetting read/write options
        inline IOobject
        (
            const IOobject& io,
            IOobjectOption::readOption rOpt,
            IOobjectOption::writeOption wOpt
        );

        //- Copy construct, resetting register option
        inline IOobject
        (
            const IOobject& io,
            IOobjectOption::registerOption regOpt
        );


        //- Clone
        autoPtr<IOobject> clone() const
        {
            return autoPtr<IOobject>::New(*this);
        }

        //- Clone resetting registry
        autoPtr<IOobject> clone(const objectRegistry& registry) const
        {
            return autoPtr<IOobject>::New(*this, registry);
        }


    // Member Functions

    // General Access

        //- Return the local objectRegistry
        const objectRegistry& db() const noexcept;

        //- Return Time associated with the objectRegistry
        const Time& time() const noexcept;

        //- Return the object name
        inline const word& name() const noexcept;

        //- Return name of the class name read from header
        inline const word& headerClassName() const noexcept;

        //- Return non-constant access to the class name read from header
        inline word& headerClassName() noexcept;

        //- Return the optional note
        inline const string& note() const noexcept;

        //- Modifiable access to the optional note
        inline string& note() noexcept;

        //- Rename the object
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }

        //- The sizeof (label) in bytes, possibly read from the header
        inline unsigned labelByteSize() const noexcept;

        //- The sizeof (scalar) in bytes, possibly read from the header
        inline unsigned scalarByteSize() const noexcept;

        //- Clear various bits (headerClassName, note, sizeof...)
        //- that would be obtained when reading from a file.
        //  \param newName if non-null, optionally rename the IOobject
        void resetHeader(const word& newName = word::null);


    // Checks

        //- True if headerClassName() is non-empty (after reading)
        inline bool hasHeaderClass() const noexcept;

        //- Check if headerClassName() equals Type::typeName
        template<class Type>
        inline bool isHeaderClass() const;

        //- Same as isHeaderClass()
        template<class Type>
        bool isHeaderClassName() const { return isHeaderClass<Type>(); }


    // Meta-data

        //- Return pointer to meta-data (if any) or nullptr
        virtual const dictionary* findMetaData() const noexcept;


    // Path components

        //- Return group (extension part of name)
        inline word group() const;

        //- Return member (name without the extension)
        inline word member() const;

        //- Return the Time::rootPath()
        const fileName& rootPath() const noexcept;

        //- Return the Time::caseName()
        const fileName& caseName() const noexcept;

        //- Return the Time::globalCaseName()
        const fileName& globalCaseName() const noexcept;

        //- Read access to instance path component
        inline const fileName& instance() const noexcept;

        //- Write access to instance path component
        inline fileName& instance() noexcept;

        //- Read access to local path component
        inline const fileName& local() const noexcept;

        //- The complete path for the object (with instance, local,...).
        fileName path() const;

        //- The complete global path for the object (with instance, local,...)
        fileName globalPath() const;

        //- The complete path with alternative instance and local
        fileName path
        (
            const word& instance,
            const fileName& local = fileName::null
        ) const;

        //- The complete global path with alternative instance and local
        fileName globalPath
        (
            const word& instance,
            const fileName& local = fileName::null
        ) const;

        //- The complete path + object name
        inline fileName objectPath() const;

        //- The complete global path + object name
        inline fileName globalObjectPath() const;

        //- The object path relative to the root
        fileName objectRelPath() const;

        //- Redirect to fileHandler filePath, searching locally.
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        fileName localFilePath
        (
            const word& typeName,
            const bool search=true
        ) const;

        //- Redirect to fileHandler filePath, searching up if in parallel.
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        fileName globalFilePath
        (
            const word& typeName,
            const bool search=true
        ) const;


    // Reading

        //- Parse 'FoamFile' header contents and set the IOobject
        //- characteristics and return the stream characteristics.
        IOstreamOption parseHeader(const dictionary& headerDict);

        //- Read header ('FoamFile' dictionary) and set the
        //- IOobject and stream characteristics.
        bool readHeader(Istream& is);

        //- Read header (the 'FoamFile' dictionary) and set the
        //- IOobject and stream characteristics.
        //  Saves the header content in the given dictionary.
        bool readHeader(dictionary& headerDict, Istream& is);

        //- Read header (respects is_globalIOobject trait) and check its info.
        template<class Type>
        bool typeHeaderOk
        (
            //! Check headerClassName against the type-name
            const bool checkType = true,
            //! Also search previous instances if not found at current instance
            const bool search = true,
            //! Report any check-type failures
            const bool verbose = true
        );

        //- Call localFilePath or globalFilePath for given type
        //- depending on its is_globalIOobject trait.
        template<class Type>
        fileName typeFilePath(const bool search = true) const;

        //- Helper: warn that type does not support re-reading
        template<class Type>
        void warnNoRereading() const;


    // Writing

        //- Write the standard OpenFOAM file/dictionary banner
        //  Optionally without editor syntax hint (eg, for logs)
        static Ostream& writeBanner(Ostream& os, const bool noSyntaxHint=false);

        //- Write the standard file section divider
        static Ostream& writeDivider(Ostream& os);

        //- Write the standard end file divider
        static Ostream& writeEndDivider(Ostream& os);

        //- Write header with current type()
        bool writeHeader(Ostream& os) const;

        //- Write header with override of type
        bool writeHeader(Ostream& os, const word& objectType) const;

        //- Write header into a dictionary with current type()
        //- and given output format
        void writeHeader(dictionary& dict, IOstreamOption streamOpt) const;

        //- Write header into a dictionary with override of type
        //- and given output format
        void writeHeader
        (
            dictionary& dict,
            const word& objectType,
            IOstreamOption streamOpt
        ) const;


    // Error Handling

        //- Did last readHeader() succeed?
        inline bool good() const noexcept;

        //- Did last readHeader() fail?
        inline bool bad() const noexcept;


    // Info

        //- Return info proxy,
        //- for printing information to a stream
        InfoProxy<IOobject> info() const noexcept { return *this; }


    // Member Operators

        //- Copy assignment, copies all values (except the registry)
        void operator=(const IOobject& io);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for \c void always returns true (no headerClassName check).
template<>
inline bool IOobject::isHeaderClass<void>() const
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOobjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOobjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
