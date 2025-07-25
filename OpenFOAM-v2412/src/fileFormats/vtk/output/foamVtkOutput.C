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

\*---------------------------------------------------------------------------*/

#include "foamVtkOutput.H"

#include "foamVtkFormatter.H"
#include "foamVtkAsciiFormatter.H"
#include "foamVtkBase64Formatter.H"
#include "foamVtkAppendBase64Formatter.H"
#include "foamVtkAppendRawFormatter.H"
#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkLegacyRawFormatter.H"
#include "foamVersion.H"
#include "typeInfo.H"
#include "globalIndex.H"
#include "instant.H"
#include "Fstream.H"
#include "Pstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::vtk::formatter>
Foam::vtk::newFormatter(std::ostream& os, unsigned prec)
{
    return autoPtr<vtk::formatter>::NewFrom<vtk::asciiFormatter>(os, prec);
}


Foam::autoPtr<Foam::vtk::formatter>
Foam::vtk::newFormatter
(
    std::ostream& os,
    const enum formatType fmtType,
    unsigned prec
)
{
    autoPtr<vtk::formatter> fmt;

    switch (fmtType)
    {
        case formatType::INLINE_ASCII:
            fmt.reset(new vtk::asciiFormatter(os, prec));
            break;

        case formatType::INLINE_BASE64:
            fmt.reset(new vtk::base64Formatter(os));
            break;

        case formatType::APPEND_BASE64:
            fmt.reset(new vtk::appendBase64Formatter(os));
            break;

        case formatType::APPEND_BINARY:
            fmt.reset(new vtk::appendRawFormatter(os));
            break;

        case formatType::LEGACY_ASCII:
            fmt.reset(new vtk::legacyAsciiFormatter(os, prec));
            break;

        case formatType::LEGACY_BINARY:
            fmt.reset(new vtk::legacyRawFormatter(os));
            break;
    }

    return fmt;
}


void Foam::vtk::writeIdentity
(
    vtk::formatter& fmt,
    const label len,
    label start
)
{
    // No nComponents for label, can use fmt.write() directly
    for (label i=0; i < len; ++i)
    {
        fmt.write(start);
        ++start;
    }
}


void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<uint8_t>& values
)
{
    // No nComponents for char, can use fmt.write() directly
    for (const uint8_t val : values)
    {
        fmt.write(val);
    }
}


void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const labelUList& values,
    const globalIndex& procOffset
)
{
    // Gather sizes (offsets irrelevant)
    const globalIndex procAddr(globalIndex::gatherOnly{}, values.size());

    if (Pstream::master())
    {
        // Write master data - with value offset
        const label offsetId = procOffset.localStart(0);
        for (const label val : values)
        {
            vtk::write(fmt, val + offsetId);
        }

        // Receive and write
        DynamicList<label> recvData(procAddr.maxNonLocalSize());

        for (const label proci : procAddr.subProcs())
        {
            const label procSize = procAddr.localSize(proci);

            if (procSize)
            {
                recvData.resize_nocopy(procSize);
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    recvData.data_bytes(),
                    recvData.size_bytes()
                );

                // With value offset
                const label offsetId = procOffset.localStart(proci);
                for (const label val : recvData)
                {
                    vtk::write(fmt, val + offsetId);
                }
            }
        }
    }
    else
    {
        if (values.size())
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                values.cdata_bytes(),
                values.size_bytes()
            );
        }
    }
}


// * * * * * * * * * * * * * * Legacy Functions  * * * * * * * * * * * * * * //

void Foam::vtk::legacy::fileHeader
(
    std::ostream& os,
    const std::string& title,
    bool binary
)
{
    // Line 1:
    os  << "# vtk DataFile Version 2.0" << nl;

    // OR
    // os  << "# vtk DataFile Version 5.1" << nl;

    // Line 2: title

    const auto truncate = title.find('\n');

    if (title.empty() || 0 == truncate)
    {
        // Avoid an empty title
        os << "File generated by OpenFOAM " << foamVersion::api << nl;
    }
    else if (std::string::npos == truncate)
    {
        os << title << nl;
    }
    else
    {
        os << title.substr(0, truncate) << nl;
    }

    // Line 3: format
    os  << (binary ? "BINARY" : "ASCII") << nl;
}


void Foam::vtk::legacy::fileHeader
(
    vtk::formatter& fmt,
    const std::string& title,
    const std::string& contentType
)
{
    std::ostream& os = fmt.os();

    legacy::fileHeader(os, title, isType<legacyRawFormatter>(fmt));
    if (contentType.size())
    {
        os << "DATASET " << contentType.c_str() << nl;
    }
}


// ************************************************************************* //
