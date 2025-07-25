/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::parsing

Description
    Collection of static functions and data related to parsing
    and an isolated namespace for lexers, parsers, scanners.

SourceFiles
    parsing.C
    parsingI.H

\*---------------------------------------------------------------------------*/

#ifndef parsing_H
#define parsing_H

#include "Enum.H"
#include <cerrno>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace parsing Declaration
\*---------------------------------------------------------------------------*/

namespace parsing
{
    //- Enumeration for possible parsing error
    enum class errorType
    {
        NONE = 0,         //!< No error encountered
        GENERAL = 1,      //!< General parsing error
        RANGE = 2,        //!< Range error
        TRAILING = 3,     //!< Trailing content detected
    };


    //- Strings corresponding to the errorType
    extern const Foam::Enum<errorType> errorNames;

    //- Sanity check after strtof, strtod, etc.
    //  Should set errno = 0 prior to the conversion.
    inline errorType checkConversion(const char* buf, char* endptr);


} // End namespace parsing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parsingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
