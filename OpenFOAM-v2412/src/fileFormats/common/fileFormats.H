/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Namespace
    Foam::fileFormats

Description
    Namespace to isolate specifics for file formats,
    and some common utilities.

SourceFiles
    fileFormats.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileFormats_H
#define Foam_fileFormats_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fileFormats
{

//- Find "formatOptions" in a top-level dictionary.
//- Extract and merge 'default' + formatName values.
//
//  \returns dictionary of merged formatOptions
dictionary getFormatOptions
(
    //! The top-level dictionary to search
    const dictionary& dict,

    //! The format name. Eg, \c ensight
    const word& formatName,

    //! Dictionary sub-entry to search for
    const word& entryName = "formatOptions"
);

//- Find "formatOptions" in a top-level dictionary,
//- and optional override dictionary.
//- Extract and merge 'default' + formatName values.
//
//  \returns dictionary of merged formatOptions
dictionary getFormatOptions
(
    //! The top-level dictionary to search
    const dictionary& dict,

    //! Additional dictionary to search
    const dictionary& altDict,

    //! The format name. Eg, \c ensight
    const word& formatName,

    //! Dictionary sub-entry to search for
    const word& entryName = "formatOptions"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
