/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::vtk::formatter

Description
    Abstract class for a VTK output stream formatter.

    Includes simple support for writing XML elements.
    By default uses single-quoting for XML attributes.

SourceFiles
    foamVtkFormatter.C
    foamVtkFormatterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_formatter_H
#define Foam_vtk_formatter_H

#include "int.H"
#include "label.H"
#include "uint64.H"
#include "direction.H"
#include "word.H"
#include "List.H"
#include "DynamicList.H"
#include "foamVtkCore.H"
#include "foamVtkPTraits.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

class outputOptions;

/*---------------------------------------------------------------------------*\
                       Class vtk::formatter Declaration
\*---------------------------------------------------------------------------*/

class formatter
{
public:

    //- Quoting for XML attributes
    enum quoteChar : char
    {
        DOUBLE_QUOTE = '\"',  //!< Double-quote XML attributes
        SINGLE_QUOTE = '\'',  //!< Single-quote XML attributes
    };


protected:

    // Private Data

        //- The output stream for the formatter
        std::ostream& os_;

        //- LIFO stack of current XML tags
        DynamicList<word> xmlTags_;

        //- Tag open/closed/ended state
        mutable bool inTag_;

        //- Quoting character for XML attributes
        char quote_;


protected:

    // Protected Member Functions

        //- Can write XML key/value attribute pair when inside a tag.
        //- Emit warning and return false if this condition is not met.
        bool canWriteAttr(const word& k) const;

        //- Can write tag-like top-level content (eg, comment, ...)
        //- when not already inside a tag.
        //- Emit warning and return false if this condition is not met.
        bool canWriteToplevel(const char* what) const;

        //- Write XML key/value attribute pair (implementation).
        template<class Type>
        inline void writeAttr(const word& k, const Type& v);

        //- No-op write XML attribute (for templating code).
        //  \return formatter for chaining
        inline formatter& xmlAttr();

        //- No-op XML comment loop (for templating code).
        inline void xmlCommentLoop();

        //- Loop/output XML comments
        template<class... Args>
        inline void xmlCommentLoop(const std::string& text, Args&&... args);

        //- Open XML tag (implementation), checking if not already inside
        //- another tag.
        //- Emit warning and return false if this condition is not met.
        bool openTagImpl(const word& tagName);


    // Constructors

        //- Construct and attach to an output stream
        inline explicit formatter(std::ostream& os);


public:

    // Public typedefs

        //- The header data is vtk UInt64
        typedef uint64_t headerType;

        //- Out of range position or size.
        static constexpr uint64_t npos = uint64_t(-1);


    //- Destructor
    virtual ~formatter() = default;


    // Member Functions

    // Access

        //- Access to the underlying output stream
        inline std::ostream& os() noexcept;

        //- The format-type / output options.
        virtual const vtk::outputOptions& opts() const = 0;

        //- Name for the XML output type or the legacy output type.
        virtual const char* name() const = 0;

        //- Name for the XML append encoding
        virtual const char* encoding() const = 0;

        //- Change quoting char for XML attributes (default: SINGLE_QUOTE)
        void quoting(const quoteChar quote);

        //- Increase the append data offset by numbytes and sizeof(uint64_t).
        //  The additional (uint64_t) size information is consistent with
        //  writeSize()
        //
        //  \return The previous data offset or formatter::npos for formats
        //     that do not support appending data.
        virtual uint64_t offset(const uint64_t numbytes);

        //- The encoded length for binary output is pass-through.
        virtual std::size_t encodedLength(std::size_t n) const;

        //- Write leading size for binary output
        //  \return True if used by the formatter.
        virtual bool writeSize(const uint64_t numbytes) = 0;

        virtual void write(const uint8_t val) = 0;
        virtual void write(const label val)   = 0;
        virtual void write(const float val)   = 0;
        virtual void write(const double val)  = 0;

        //- Flush encoding, write newline etc.
        virtual void flush() = 0;


    // General Output

        //- Add indenting according to the current XML tag depth
        //  Two spaces per depth.
        inline void indent();

        //- Add indenting of n spaces.
        inline void indent(label n);

        //- Write XML header
        //  \return formatter for chaining
        inline formatter& xmlHeader();

        //- Write XML comment (at the current indentation level)
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlComment(const std::string& text, Args&&... args);


        //- Start an XML tag, optionally with attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& openTag(const word& tagName, Args&&... args);

        //- Start an XML tag, optionally with attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& openTag(vtk::fileTag t, Args&&... args);

        //- Finish an XML tag, optional as an empty container.
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        formatter& closeTag(const bool isEmpty = false);

        //- An end XML tag, optional with sanity check
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        formatter& endTag(const word& tagName = word::null);

        //- An end XML tag with sanity check
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        inline virtual formatter& endTag(vtk::fileTag t);

        //- Write XML tag without any attributes. Combines openTag/closeTag.
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& tag(const word& t, Args&&... args);

        //- Write XML tag without any attributes. Combines openTag/closeTag.
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& tag(vtk::fileTag t, Args&&... args);


        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //- the contentType itself.
        //  \param leaveOpen Leave tag open for additional attributes.
        //  \return formatter for chaining
        formatter& beginVTKFile
        (
            const word& contentType,
            const word& contentVersion,
            const bool leaveOpen = false
        );

        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //- the contentType itself.
        //  \param leaveOpen Leave tag open for additional attributes.
        //  \return formatter for chaining
        inline formatter& beginVTKFile
        (
            vtk::fileTag contentType,
            const word& contentVersion,
            const bool leaveOpen = false
        );

        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //- the contentType itself.
        //  \param leaveOpen Leave tag open for additional attributes.
        //  \return formatter for chaining
        inline formatter& beginVTKFile
        (
            vtk::fileTag contentType,
            const bool leaveOpen = false
        );

        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //- the contentType itself.
        //  \param leaveOpen Leave tag open for additional attributes.
        //  \return formatter for chaining
        template<vtk::fileTag ContentType>
        inline formatter& beginVTKFile(bool leaveOpen = false);

        //- Add a "AppendedData" XML tag with the current encoding and output
        //- the requisite '_' prefix.
        //  \return formatter for chaining
        formatter& beginAppendedData();

        //- Begin "Block" XML section.
        //  \param index The index of the block
        //  \param name The name of the block (ignored if empty)
        //  \return formatter for chaining
        formatter& beginBlock(label index, std::string name = "");

        //- End "Block" XML section.
        //  \return formatter for chaining
        inline formatter& endBlock();

        //- Begin "Piece" XML section.
        //  \param index The index of the piece
        //  \param name The name of the piece (ignored if empty)
        //  \return formatter for chaining
        formatter& beginPiece(label index, std::string name = "");

        //- End "Piece" XML section.
        //  \return formatter for chaining
        inline virtual formatter& endPiece();

        //- Insert a single "DataSet" XML entry tag.
        //  \param index The index of the DataSet
        //  \param file The file name for the data (ignored if empty)
        //  \param autoName The name for the data extracted from the file name
        //      (without extension)
        //  \return formatter for chaining
        formatter& DataSet
        (
            label index,
            std::string file = "",
            bool autoName = true
        );

        //- Insert a single "DataSet" XML entry tag.
        //  \param index The index of the DataSet
        //  \param file The file name for the data (ignored if empty)
        //  \param name The name for the dataset
        //  \return formatter for chaining
        formatter& DataSet
        (
            label index,
            std::string file,
            std::string name
        );

        //- Begin "DataArray" XML section.
        //
        //  \param dataName The name of the DataArray
        //  \param payLoad  Additional payLoad information to increment
        //     the offset for an append formatter and add the "offset"
        //     attribute accordingly.
        //  \param leaveOpen Leave tag open for additional attributes.
        //
        //  \return formatter for chaining
        template<class Type, direction nComp=1, int nTuple=0>
        formatter& beginDataArray
        (
            const word& dataName,
            uint64_t payLoad = npos,
            bool leaveOpen = false
        );

        //- Begin "DataArray" XML section.
        //
        //  \param dataName The name of the DataArray as an enumeration
        //  \param payLoad  Additional payLoad information to increment
        //     the offset for an append formatter and add the "offset"
        //     attribute accordingly.
        //  \param leaveOpen Leave tag open for additional attributes.
        //
        //  \return formatter for chaining
        template<class Type, direction nComp=1, int nTuple=0>
        inline formatter& beginDataArray
        (
            const vtk::dataArrayAttr& dataName,
            uint64_t payLoad = npos,
            bool leaveOpen = false
        );

        //- End "DataArray" XML section
        //  \return formatter for chaining
        inline virtual formatter& endDataArray();

        //- Insert a single "PDataArray" XML entry tag.
        //  For some entries, the name is optional.
        //  \return formatter for chaining
        template<class Type, direction nComp=1, int nTuple=0>
        formatter& PDataArray(const word& dataName);


        //- Begin "FieldData" XML section.
        inline formatter& beginFieldData();

        //- Begin "CellData" XML section.
        inline formatter& beginCellData();

        //- Begin "PointData" XML section.
        inline formatter& beginPointData();

        //- End "FieldData" XML section.
        inline virtual formatter& endFieldData();

        //- End "CellData" XML section.
        inline virtual formatter& endCellData();

        //- End "PointData" XML section.
        inline virtual formatter& endPointData();


        //- End "AppendedData" XML section
        //  \return formatter for chaining
        formatter& endAppendedData();

        //- End "VTKFile" XML section.
        //  \return formatter for chaining
        inline virtual formatter& endVTKFile();


        //- Emit "TimeValue" for FieldData (name as per Catalyst output)
        formatter& writeTimeValue(scalar timeValue);


    // XML Attributes

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const word& k,
            const std::string& v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const word& k,
            const int32_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const word& k,
            const int64_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const word& k,
            const uint64_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const word& k,
            const scalar v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& k,
            const std::string& v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& k,
            const int32_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& k,
            const int64_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& k,
            const uint64_t v,
            Args&&... args
        );

        //- Pair-wise write of XML key/value attributes
        //  \return formatter for chaining
        template<class... Args>
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& k,
            const scalar v,
            Args&&... args
        );


    // Housekeeping

        //- Open "DataArray" XML tag and leave open (requires a closeTag).
        //  \deprecated Use beginDataArray instead (SEPT-2018)
        template<class Type, direction nComp=1, int nTuple=0>
        formatter& openDataArray(const word& dataName)
        {
            return beginDataArray<Type, nComp, nTuple>
            (
                dataName, formatter::npos, true
            );
        }

        //- Open "DataArray" XML tag and leave open (requires a closeTag).
        //  \deprecated Use beginDataArray instead (SEPT-2018)
        template<class Type, direction nComp=1, int nTuple=0>
        formatter& openDataArray(const vtk::dataArrayAttr& dataName)
        {
            return beginDataArray<Type, nComp, nTuple>
            (
                dataName, formatter::npos, true
            );
        }
};


// Global Functions

//- Commonly used calculation for header and payload sizes
template<class Type, direction nComp=1>
inline uint64_t sizeofData(label count)
{
    return (count * nComp * sizeof(Type));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkFormatterI.H"

#ifdef NoRepository
    #include "foamVtkFormatterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
