/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
\*---------------------------------------------------------------------------*/

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::vtk::formatter::canWriteAttr(const word& k) const
{
    if (!inTag_)
    {
        WarningInFunction
            << "xml attribute '" << k << "' but not inside a tag!" << endl;
    }

    return inTag_;
}


bool Foam::vtk::formatter::canWriteToplevel(const char* what) const
{
    if (inTag_)
    {
        WarningInFunction
            << "Cannot add " << what << " inside a tag!"
            << endl;
    }

    return !inTag_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::vtk::formatter::quoting(const quoteChar quote)
{
    quote_ = quote;
}


uint64_t Foam::vtk::formatter::offset(const uint64_t)
{
    return formatter::npos;
}


std::size_t Foam::vtk::formatter::encodedLength(std::size_t n) const
{
    return n;
}


bool Foam::vtk::formatter::openTagImpl(const word& tagName)
{
    if (inTag_)
    {
        WarningInFunction
            << "open xml tag '" << tagName << "', but already within a tag!"
            << endl;

        return false;
    }

    // Emit, before changing the stack or the state.
    indent();
    os_ << '<' << tagName;

    // Add to the stack and change the state.
    xmlTags_.append(tagName);
    inTag_ = true;

    return true;
}


Foam::vtk::formatter& Foam::vtk::formatter::closeTag(const bool isEmpty)
{
    if (!inTag_)
    {
        WarningInFunction
            << "attempt to close xml tag, but not within a tag!"
            << endl;
    }
    else
    {
        // Change the state
        inTag_ = false;

        if (isEmpty)
        {
            // Eg, <tag ... />
            xmlTags_.pop_back();
            os_ << " /";
        }
        os_ << '>' << nl;
    }

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::endTag(const word& tagName)
{
    const word curr(std::move(xmlTags_.back()));
    xmlTags_.pop_back();
    indent();

    if (inTag_)
    {
        WarningInFunction
            << "adding xml endTag '" << curr
            << "' but already in another tag!"
            << endl;

        // Also suppress further output, or not?
    }

    // Verify expected end tag
    if (!tagName.empty() && tagName != curr)
    {
        WarningInFunction
            << "expecting to end xml tag '" << tagName
            << "' but found '" << curr << "' instead"
            << endl;
    }

    os_  << "</" << curr << '>' << nl;

    inTag_ = false;

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::beginVTKFile
(
    const word& contentType,
    const word& contentVersion,
    const bool leaveOpen
)
{
    openTag(vtk::fileTag::VTK_FILE);
    xmlAttr("type",        contentType);
    xmlAttr("version",     contentVersion);
    xmlAttr("byte_order",  vtkPTraits<Foam::endian>::typeName);
    xmlAttr("header_type", vtkPTraits<headerType>::typeName);
    closeTag();

    openTag(contentType);
    if (!leaveOpen)
    {
        closeTag();
    }

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::beginAppendedData()
{
    openTag("AppendedData");
    xmlAttr("encoding", encoding());
    closeTag();
    os_ << '_';

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::endAppendedData()
{
    flush();     // Flush any pending encoded content
    os_ << nl;   // Ensure clear separation from content
    return endTag("AppendedData");
}


Foam::vtk::formatter& Foam::vtk::formatter::beginBlock
(
    label index,
    std::string name
)
{
    openTag(vtk::fileTag::BLOCK);
    if (index >= 0)
    {
        xmlAttr("index", index);
    }
    if (name.size())
    {
        xmlAttr("name", name);
    }
    closeTag();

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::beginPiece
(
    label index,
    std::string name
)
{
    openTag(vtk::fileTag::PIECE);
    if (index >= 0)
    {
        xmlAttr("index", index);
    }
    if (name.size())
    {
        xmlAttr("name", name);
    }
    closeTag();

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::DataSet
(
    label index,
    std::string file,
    bool autoName
)
{
    openTag(vtk::fileTag::DATA_SET);

    if (index >= 0)
    {
        xmlAttr("index", index);
    }
    if (file.size())
    {
        if (autoName)
        {
            xmlAttr("name", fileName::stem(file));
        }
        xmlAttr("file", file);
    }
    closeTag(true);  // Empty tag. ie, <DataSet ... />

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::DataSet
(
    label index,
    std::string file,
    std::string name
)
{
    openTag(vtk::fileTag::DATA_SET);

    if (index >= 0)
    {
        xmlAttr("index", index);
    }
    if (name.size())
    {
        xmlAttr("name", name);
    }
    if (file.size())
    {
        xmlAttr("file", file);
    }
    closeTag(true);  // Empty tag. ie, <DataSet ... />

    return *this;
}


Foam::vtk::formatter& Foam::vtk::formatter::writeTimeValue(scalar timeValue)
{
    // Emit "TimeValue" as FieldData
    // NumberOfTuples="1" (required!)

    uint64_t payLoad = vtk::sizeofData<float>(1);

    beginDataArray<float,1,1>("TimeValue");
    writeSize(payLoad);

    write(timeValue);
    flush();

    endDataArray();

    return *this;
}


// ************************************************************************* //
