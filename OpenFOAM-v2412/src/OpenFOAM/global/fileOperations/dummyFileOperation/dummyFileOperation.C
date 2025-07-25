/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "dummyFileOperation.H"
#include "dummyISstream.H"
#include "Fstream.H"
#include "objectRegistry.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeName(dummyFileOperation);

// No runtime selection
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::dummyFileOperation::dummyFileOperation
(
    bool verbose
)
:
    fileOperation
    (
        // Use COMM_SELF for now, but COMM_NULL is probably more appropriate
        Tuple2<label, labelList>
        (
            UPstream::commSelf(),
            fileOperation::getGlobalIORanks()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::dummyFileOperation::~dummyFileOperation()
{}


// * * * * * * * * * * * * * Filesystem Operations * * * * * * * * * * * * * //

bool Foam::fileOperations::dummyFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    NotImplemented;
    return false;
}


mode_t Foam::fileOperations::dummyFileOperation::mode
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return mode_t(0);
}


Foam::fileName::Type Foam::fileOperations::dummyFileOperation::type
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return fileName::Type::UNDEFINED;
}


bool Foam::fileOperations::dummyFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    NotImplemented;
    return false;
}


off_t Foam::fileOperations::dummyFileOperation::fileSize
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return off_t(0);
}


time_t Foam::fileOperations::dummyFileOperation::lastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return time_t(0);
}


double Foam::fileOperations::dummyFileOperation::highResLastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    NotImplemented;
    return 0;
}



Foam::fileNameList Foam::fileOperations::dummyFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
) const
{
    NotImplemented;
    return fileNameList();
}


bool Foam::fileOperations::dummyFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::rm
(
    const fileName& fName
) const
{
    NotImplemented;
    return false;
}


bool Foam::fileOperations::dummyFileOperation::rmDir
(
    const fileName& dir,
    const bool silent,
    const bool emptyOnly
) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::dummyFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName,
    const bool search
) const
{
    NotImplemented;
    return fileName();
}


Foam::fileName Foam::fileOperations::dummyFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io,
    const bool search
) const
{
    NotImplemented;
    return fileName();
}


Foam::fileNameList Foam::fileOperations::dummyFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    NotImplemented;
    return fileNameList();
}


bool Foam::fileOperations::dummyFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    NotImplemented;
    io.resetHeader();
    return false;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::dummyFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const word& typeName,
    const bool readOnProc
) const
{
    NotImplemented;
    return autoPtr<ISstream>(new dummyISstream());
}


bool Foam::fileOperations::dummyFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstreamOption::streamFormat format,
    const word& typeName
) const
{
    NotImplemented;
    return false;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::dummyFileOperation::NewIFstream
(
    const fileName& filePath
) const
{
    NotImplemented;
    return autoPtr<ISstream>(new dummyISstream());
}


Foam::autoPtr<Foam::OSstream>
Foam::fileOperations::dummyFileOperation::NewOFstream
(
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    NotImplemented;
    return autoPtr<OSstream>(new OFstream(nullptr));  // ie, /dev/null
}


Foam::autoPtr<Foam::OSstream>
Foam::fileOperations::dummyFileOperation::NewOFstream
(
    IOstreamOption::atomicType atomic,
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    NotImplemented;
    return autoPtr<OSstream>(new OFstream(nullptr));  // ie, /dev/null
}


// ************************************************************************* //
