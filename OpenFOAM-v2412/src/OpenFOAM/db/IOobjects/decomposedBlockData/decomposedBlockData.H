/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::decomposedBlockData

Description
    The decomposedBlockData comprise a \c List\<char\> for each output
    processor, typically with IO on the master processor only.

    For decomposedBlockData, we make a distinction between the container
    description and the individual block contents.

    The \b FoamFile header specifies the container characteristics and thus
    has \c class = \c %decomposedBlockData and normally \c format = \c binary.
    This description refers to the \em entire file container, not the
    individual blocks.

    Each processor block is simply a binary chunk of characters and the
    first block also contains the header description for all of the blocks.
    For example,
\verbatim
FoamFile
{
    version     2.0;
    format      binary;
    arch        "LSB;label=32;scalar=64";
    class       decomposedBlockData;
    location    "constant/polyMesh";
    object      points;
}

// processor0
NCHARS
(FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
...content...
)

// processor1
NCHARS
(...content...)

...
\endverbatim


SourceFiles
    decomposedBlockData.C
    decomposedBlockDataHeader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_decomposedBlockData_H
#define Foam_decomposedBlockData_H

#include "regIOobject.H"
#include "ISstream.H"
#include "OSstream.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class decomposedBlockData Declaration
\*---------------------------------------------------------------------------*/

class decomposedBlockData
:
    public regIOobject
{
    // Private Functions

        //- Helper: write content for FoamFile IOobject header
        static void writeHeaderContent
        (
            Ostream& os,
            IOstreamOption streamOptContainer,
            const word& objectType,
            const string& note,
            const fileName& location,
            const word& objectName
        );


protected:

    // Protected Data

        //- Type to use for gather
        const UPstream::commsTypes commsType_;

        //- Communicator for all parallel comms
        const label comm_;

        //- The block content
        List<char> contentData_;


    // Protected Member Functions

        //- Helper: determine number of processors whose recvSizes fits
        //- into maxBufferSize
        static label calcNumProcs
        (
            const label comm,
            const off_t maxBufferSize,
            const labelUList& recvSizes,
            const label startProci
        );

        //- Read data into *this. ISstream is only valid on master.
        static bool readBlocks
        (
            const label comm,
            autoPtr<ISstream>& isPtr,
            List<char>& contentChars,
            const UPstream::commsTypes commsType
        );

        //- Helper: skip a block of (binary) character data
        static bool skipBlockEntry(Istream& is);

public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("decomposedBlockData");


    // Constructors

        //- Construct given an IOobject
        decomposedBlockData
        (
            const label comm,
            const IOobject& io,
            const UPstream::commsTypes = UPstream::commsTypes::scheduled
        );


    //- Destructor
    virtual ~decomposedBlockData() = default;


    // Member Functions

        //- Read object
        virtual bool read();

        //- Write separated content (assumes content is the serialised data)
        //  The serialised master data should also contain a FoamFile header
        virtual bool writeData(Ostream& os) const;

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


    // Helpers

        //- True if object type is a known collated type
        static bool isCollatedType(const word& objectType);

        //- True if object header class is a known collated type
        static bool isCollatedType(const IOobject& io);

        //- Extract number of decomposedBlockData block entries, optionally
        //- with an upper limit.
        //- The input stream should be in a rewound state
        //- (or only have read the header) before calling.
        static label getNumBlocks(Istream& is, const label maxNumBlocks = -1);

        //- True if the given block number (starts at 0) has a corresponding
        //- decomposedBlockData block entry.
        //- The input stream should be in a rewound state
        //- (or only have read the header) before calling.
        //  This will be faster than checking against getNumBlocks()
        //  since it can potentially exit without scanning the entire file.
        static bool hasBlock(Istream& is, const label blockNumber);

        //- Read header as per IOobject with additional handling of
        //- decomposedBlockData
        static bool readHeader(IOobject& io, Istream& is);

        //- Helper: write FoamFile IOobject header
        static void writeHeader
        (
            Ostream& os,
            IOstreamOption streamOptContainer,
            const word& objectType,
            const string& note,
            const fileName& location,
            const word& objectName,
            const dictionary& extraEntries
        );

        //- Helper: write FoamFile IOobject header
        static void writeHeader
        (
            Ostream& os,
            IOstreamOption streamOptData,
            const IOobject& io
        );

        //- Helper: generate additional entries for FoamFile header
        static void writeExtraHeaderContent
        (
            dictionary& dict,
            IOstreamOption streamOptData,
            const IOobject& io
        );

        //- Helper: read block of (binary) character data
        static bool readBlockEntry
        (
            Istream& is,
            List<char>& charData
        );

        //- Helper: write block of (binary) character data
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            const char* str,
            const size_t len
        );

        //- Helper: write block of (binary) character data
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            const UList<char>& s
        )
        {
            return writeBlockEntry(os, blocki, s.cdata(), s.size_bytes());
        }

        //- Helper: write block of (binary) character content
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            const std::string& s
        )
        {
            return writeBlockEntry(os, blocki, s.data(), s.size());
        }

        //- Helper: write block of (binary) character data
        //  \return -1 on error
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            IOstreamOption streamOptData,
            const regIOobject& io,
            const label blocki,
            const bool withLocalHeader
        );

        //- Read selected block (non-seeking) + header information
        static autoPtr<ISstream> readBlock
        (
            const label blocki,
            ISstream& is,
            IOobject& headerIO
        );

        //- Read master header information (into headerIO) and return
        //- data in stream. Note: isPtr is only valid on master.
        static autoPtr<ISstream> readBlocks
        (
            const label comm,
            const fileName& fName,
            autoPtr<ISstream>& isPtr,
            IOobject& headerIO,
            const UPstream::commsTypes commsType
        );

        //- Helper: gather single label. Note: using native Pstream.
        //  datas sized with num procs but undefined contents on
        //  slaves
        static void gather
        (
            const label comm,
            const label data,
            labelList& datas
        );

        //- Helper: gather data from (subset of) slaves.
        //
        //  Returns:
        //  - recvData : received data
        //  - recvOffsets : offset in data. recvOffsets is nProcs+1
        static void gatherSlaveData
        (
            const label comm,
            const UList<char>& data,
            const labelUList& recvSizes,

            const labelRange& fromProcs,

            List<int>& recvOffsets,
            DynamicList<char>& recvData
        );

        //- Write *this. Ostream only valid on master.
        //  Returns offsets of processor blocks in blockOffset
        static bool writeBlocks
        (
            const label comm,
            autoPtr<OSstream>& osPtr,
            List<std::streamoff>& blockOffset,

            const UList<char>& masterData,

            const labelUList& recvSizes,

            // Optional slave data (on master)
            const UPtrList<SubList<char>>& slaveData,

            const UPstream::commsTypes commsType,
            const bool syncReturnState = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
