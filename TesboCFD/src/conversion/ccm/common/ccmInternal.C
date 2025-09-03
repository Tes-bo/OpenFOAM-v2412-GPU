/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "Pstream.H"
#include "ccmInternal.H" // include last to avoid any strange interactions


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ccm, 0);
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// \cond file_scope
// Return the string corresponding to the error
static const char* errorMsg(CCMIOError err)
{
    switch (err)
    {
        case kCCMIONoErr:
            return "";
            break;
        case kCCMIONoFileErr:
            return "No File Error";
            break;
        case kCCMIOPermissionErr:
            return "Permission Error";
            break;
        case kCCMIOCorruptFileErr:
            return "Corrupt File Error";
            break;
        case kCCMIOBadLinkErr:
            return "Bad Link Error";
            break;
        case kCCMIONoNodeErr:
            return "No Node Error";
            break;
        case kCCMIODuplicateNodeErr:
            return "Duplicate Node Error";
            break;
        case kCCMIOWrongDataTypeErr:
            return "Wrong Data Type Error";
            break;
        case kCCMIONoDataErr:
            return "NoData Error";
            break;
        case kCCMIOWrongParentErr:
            return "Wrong Parent Error";
            break;
        case kCCMIOBadParameterErr:
            return "Bad Parameter Error";
            break;
        case kCCMIONoMemoryErr:
            return "No Memory Error";
            break;
        case kCCMIOIOErr:
            return "IO Error";
            break;
        case kCCMIOTooManyFacesErr:
            return "Too Many Faces Error";
            break;
        case kCCMIOVersionErr:
            return "Version Error";
            break;
        case kCCMIOArrayDimensionToLargeErr:
            return "Array Dimension To Large Error";
            break;
        case kCCMIOInternalErr:
            return "Internal Error";
            break;
        default:
            return "unclassified error";
            break;
    }
}
//! \endcond


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::ccm::ccmGlobalState::assertNoError
(
    CCMIOError err,
    const char* msg
)
{
    if (err != kCCMIONoErr)
    {
        Perr<< endl
            << "--> FATAL ERROR:"
            << "\n    " << msg
            << "\n    libccmio reports -> " << errorMsg(err) << " <-\n"
            << endl;

        std::exit(1);
    }

    return (err == kCCMIONoErr);
}


bool Foam::ccm::ccmGlobalState::assertNoError
(
    CCMIOError err,
    const std::string& msg
)
{
    if (err != kCCMIONoErr)
    {
        assertNoError(err, msg.c_str());
    }

    return (err == kCCMIONoErr);
}


bool Foam::ccm::ccmGlobalState::assertNoError
(
    const char* msg
) const
{
    return assertNoError(error, msg);
}


bool Foam::ccm::ccmGlobalState::assertNoError
(
    const std::string& msg
) const
{
    return assertNoError(error, msg);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ccm::ccmGlobalState::ccmGlobalState()
:
    error(kCCMIONoErr)
{
    CCMIOInvalidateEntity(&root);
}


// ************************************************************************* //
