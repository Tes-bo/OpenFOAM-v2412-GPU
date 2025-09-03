/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vtk::outputOptions&
Foam::vtk::outputOptions::ascii(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case formatType::INLINE_BASE64:
                fmtType_ = formatType::INLINE_ASCII;
                break;

            case formatType::APPEND_BINARY:
                fmtType_ = formatType::APPEND_BASE64;
                break;

            case formatType::LEGACY_BINARY:
                fmtType_ = formatType::LEGACY_ASCII;
                break;

            default:  // No change
                break;
        }
    }
    else
    {
        switch (fmtType_)
        {
            case formatType::INLINE_ASCII:
                fmtType_ = formatType::INLINE_BASE64;
                break;

            case formatType::APPEND_BASE64:
                fmtType_ = formatType::APPEND_BINARY;
                break;

            case formatType::LEGACY_ASCII:
                fmtType_ = formatType::LEGACY_BINARY;
                break;

            default:  // No change
                break;
        }
    }

    return *this;
}


Foam::vtk::outputOptions&
Foam::vtk::outputOptions::append(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case formatType::INLINE_ASCII:
            case formatType::LEGACY_ASCII:
                fmtType_ = formatType::APPEND_BASE64;
                break;

            case formatType::INLINE_BASE64:
            case formatType::LEGACY_BINARY:
                fmtType_ = formatType::APPEND_BINARY;
                break;

            default:  // No change
                break;
        }
    }
    else
    {
        switch (fmtType_)
        {
            case formatType::APPEND_BASE64:
                fmtType_ = formatType::INLINE_ASCII;
                break;

            case formatType::APPEND_BINARY:
                fmtType_ = formatType::INLINE_BASE64;
                break;

            default:  // No change
                break;
        }
    }

    return *this;
}


Foam::vtk::outputOptions&
Foam::vtk::outputOptions::legacy(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case formatType::INLINE_ASCII:
            case formatType::APPEND_BASE64:
                fmtType_ = formatType::LEGACY_ASCII;
                break;

            case formatType::INLINE_BASE64:
            case formatType::APPEND_BINARY:
                fmtType_ = formatType::LEGACY_BINARY;
                break;

            default:  // no change
                break;
        }
    }
    else
    {
        switch (fmtType_)
        {
            case formatType::LEGACY_ASCII:
                fmtType_ = formatType::INLINE_ASCII;
                break;

            case formatType::LEGACY_BINARY:
                fmtType_ = formatType::INLINE_BASE64;
                break;

            default:  // no change
                break;
        }
    }

    return *this;
}


Foam::vtk::outputOptions&
Foam::vtk::outputOptions::precision(unsigned prec)
{
    precision_ = prec;
    return *this;
}


Foam::string Foam::vtk::outputOptions::description() const
{
    switch (fmtType_)
    {
        case formatType::INLINE_ASCII:  return "xml ascii";
        case formatType::INLINE_BASE64: return "xml base64";
        case formatType::APPEND_BASE64: return "xml-append base64";
        case formatType::APPEND_BINARY: return "xml-append binary";
        case formatType::LEGACY_ASCII:  return "legacy ascii";
        case formatType::LEGACY_BINARY: return "legacy binary";
    }

    return "";
}


// ************************************************************************* //
