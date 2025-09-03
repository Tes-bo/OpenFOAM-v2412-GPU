/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::ensight::FileName

Description
    Specification of a valid Ensight file-name.

    Spaces must be quoted,
    no '*' wildcards, no '%' (structured block continuation).

    Overall line length within case file is limited to 1024, but this is not
    yet addressed.

Note
    To assist windows usage, a colon (':') in the name is replaced with
    by an underscore.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightFileName_H
#define Foam_ensightFileName_H

#include "fileName.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ensight
{

/*---------------------------------------------------------------------------*\
                      Class ensight::FileName Declaration
\*---------------------------------------------------------------------------*/

class FileName
:
    public fileName
{
public:

    // Generated Methods

        //- Copy construct
        FileName(const FileName&) = default;

        //- Move construct
        FileName(FileName&&) = default;


    // Constructors

        //- Construct as copy of character array
        inline explicit FileName(const char* s);

        //- Copy construct from std::string
        inline explicit FileName(const std::string& s);

        //- Move construct from std::string
        inline explicit FileName(std::string&& s);


    // Member Functions

        //- Is this character valid for an ensight file-name
        inline static bool valid(char c);

        //- Strip invalid characters
        inline void stripInvalid();


    // Member Operators

        // Assignment (disabled)

        void operator=(const fileName&) = delete;
        void operator=(const word&) = delete;
        void operator=(const string&) = delete;
        void operator=(const std::string&) = delete;
        void operator=(const char*) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensight

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing for FileName is the same as string
template<> struct Hash<ensight::FileName> : string::hasher {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightFileNameI.H"

#endif

// ************************************************************************* //
