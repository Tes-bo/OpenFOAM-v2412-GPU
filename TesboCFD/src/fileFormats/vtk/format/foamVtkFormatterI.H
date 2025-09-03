/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr()
{
    return *this;
}


inline void Foam::vtk::formatter::xmlCommentLoop()
{}


template<class... Args>
inline void Foam::vtk::formatter::xmlCommentLoop
(
    const std::string& text,
    Args&&... args
)
{
    if (text.length())
    {
        indent(); indent(4);
        os_ << text << nl;
    }

    xmlCommentLoop(std::forward<Args>(args)...);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::vtk::formatter::formatter(std::ostream& os)
:
    os_(os),
    xmlTags_(),
    inTag_(false),
    quote_(SINGLE_QUOTE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline std::ostream& Foam::vtk::formatter::os() noexcept
{
    return os_;
}


inline void Foam::vtk::formatter::indent()
{
    indent(2*xmlTags_.size());
}


inline void Foam::vtk::formatter::indent(label n)
{
    while (n-- > 0)
    {
        os_ << ' ';
    }
}


inline Foam::vtk::formatter& Foam::vtk::formatter::xmlHeader()
{
    if (canWriteToplevel("xml header"))
    {
        os_ << "<?xml version='1.0'?>" << nl;
    }

    return *this;
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlComment
(
    const std::string& text,
    Args&&... args
)
{
    if (canWriteToplevel("xml comment"))
    {
        indent();
        os_ << "<!--";

        if (sizeof...(Args))
        {
            os_ << nl;

            xmlCommentLoop(text, std::forward<Args>(args)...);

            indent(); indent(2);
        }
        else
        {
            os_ << ' ' << text << ' ';
        }

        os_ << "-->" << nl;
    }

    return *this;
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::openTag
(
    const word& tagName,
    Args&&... args
)
{
    if (openTagImpl(tagName))
    {
        xmlAttr(std::forward<Args>(args)...);
    }

    return *this;
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::openTag
(
    vtk::fileTag t,
    Args&&... args
)
{
    return openTag(vtk::fileTagNames[t], std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::tag
(
    const word& t,
    Args&&... args
)
{
    openTagImpl(t);
    xmlAttr(std::forward<Args>(args)...);
    closeTag();

    return *this;
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::tag
(
    vtk::fileTag t,
    Args&&... args
)
{
    return tag(vtk::fileTagNames[t], std::forward<Args>(args)...);
}


// Begin tags

inline Foam::vtk::formatter& Foam::vtk::formatter::beginVTKFile
(
    vtk::fileTag contentType,
    const word& contentVersion,
    bool leaveOpen
)
{
    return beginVTKFile
    (
        vtk::fileTagNames[contentType],
        contentVersion,
        leaveOpen
    );
}


inline Foam::vtk::formatter& Foam::vtk::formatter::beginVTKFile
(
    vtk::fileTag contentType,
    bool leaveOpen
)
{
    return beginVTKFile
    (
        vtk::fileTagNames[contentType],
        vtk::fileContentVersions[contentType],
        leaveOpen
    );
}


template<Foam::vtk::fileTag ContentType>
inline Foam::vtk::formatter& Foam::vtk::formatter::beginVTKFile(bool leaveOpen)
{
    return beginVTKFile
    (
        vtk::fileTagNames[ContentType],
        vtk::fileContentVersions[ContentType],
        leaveOpen
    );
}


template<class Type, Foam::direction nComp, int nTuple>
inline Foam::vtk::formatter& Foam::vtk::formatter::beginDataArray
(
    const vtk::dataArrayAttr& dataName,
    uint64_t payLoad,
    bool leaveOpen
)
{
    return
        beginDataArray<Type, nComp, nTuple>
        (
            vtk::dataArrayAttrNames[dataName],
            payLoad,
            leaveOpen
        );
}


inline Foam::vtk::formatter& Foam::vtk::formatter::beginCellData()
{
    return tag(vtk::fileTag::CELL_DATA);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::beginPointData()
{
    return tag(vtk::fileTag::POINT_DATA);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::beginFieldData()
{
    return tag(vtk::fileTag::FIELD_DATA);
}


// End tags

inline Foam::vtk::formatter& Foam::vtk::formatter::endTag(vtk::fileTag t)
{
    return endTag(vtk::fileTagNames[t]);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endCellData()
{
    return endTag(vtk::fileTag::CELL_DATA);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endPointData()
{
    return endTag(vtk::fileTag::POINT_DATA);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endFieldData()
{
    return endTag(vtk::fileTag::FIELD_DATA);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endDataArray()
{
    return endTag(vtk::fileTag::DATA_ARRAY);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endBlock()
{
    return endTag(vtk::fileTag::BLOCK);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endPiece()
{
    return endTag(vtk::fileTag::PIECE);
}


inline Foam::vtk::formatter& Foam::vtk::formatter::endVTKFile()
{
    return endTag(vtk::fileTag::VTK_FILE);
}


// Attributes

template<class Type>
inline void Foam::vtk::formatter::writeAttr(const word& k, const Type& v)
{
    os_ << ' ' << k << '=' << quote_ << v << quote_;
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const std::string& v,
    Args&&... args
)
{
    if (!canWriteAttr(k)) return *this;

    writeAttr(k, v.c_str());
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const int32_t v,
    Args&&... args
)
{
    if (!canWriteAttr(k)) return *this;

    writeAttr(k, v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const int64_t v,
    Args&&... args
)
{
    if (!canWriteAttr(k)) return *this;

    writeAttr(k, v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const uint64_t v,
    Args&&... args
)
{
    if (!canWriteAttr(k)) return *this;

    writeAttr(k, v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const word& k,
    const scalar v,
    Args&&... args
)
{
    if (!canWriteAttr(k)) return *this;

    writeAttr(k, v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& k,
    const std::string& v,
    Args&&... args
)
{
    if (!canWriteAttr(vtk::fileAttrNames[k])) return *this;

    writeAttr(vtk::fileAttrNames[k], v.c_str());
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& k,
    const int32_t v,
    Args&&... args
)
{
    if (!canWriteAttr(vtk::fileAttrNames[k])) return *this;

    writeAttr(vtk::fileAttrNames[k], v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& k,
    const int64_t v,
    Args&&... args
)
{
    if (!canWriteAttr(vtk::fileAttrNames[k])) return *this;

    writeAttr(vtk::fileAttrNames[k], v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& k,
    const uint64_t v,
    Args&&... args
)
{
    if (!canWriteAttr(vtk::fileAttrNames[k])) return *this;

    writeAttr(vtk::fileAttrNames[k], v);
    return xmlAttr(std::forward<Args>(args)...);
}


template<class... Args>
inline Foam::vtk::formatter& Foam::vtk::formatter::xmlAttr
(
    const vtk::fileAttr& k,
    const scalar v,
    Args&&... args
)
{
    if (!canWriteAttr(vtk::fileAttrNames[k])) return *this;

    writeAttr(vtk::fileAttrNames[k], v);
    return xmlAttr(std::forward<Args>(args)...);
}


// ************************************************************************* //
