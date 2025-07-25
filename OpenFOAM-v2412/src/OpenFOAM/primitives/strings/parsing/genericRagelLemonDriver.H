/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::parsing::genericRagelLemonDriver

Description
    Generic interface code for Ragel/Lemon combination
    Subclasses should implement one or more process() methods.

    The scanner will often be implemented as localized lexer class.
    The parser may be embedded into the scanner as file-scope, or
    use a separate interface class.

SourceFiles
    genericRagelLemonDriver.C

\*---------------------------------------------------------------------------*/

#ifndef genericRagelLemonDriver_H
#define genericRagelLemonDriver_H

#include "error.H"
#include "className.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace parsing
{

/*---------------------------------------------------------------------------*\
                   Class genericRagelLemonDriver Declaration
\*---------------------------------------------------------------------------*/

class genericRagelLemonDriver
{
protected:

    // Protected Data

        //- Reference to the input string
        std::reference_wrapper<const std::string> content_;

        //- Start position within input string
        size_t start_;

        //- Length of input (sub)string
        size_t length_;

        //- The last known parser position
        size_t position_;


public:

    // Public Typedefs

        //- Type for linear addressing within parse content
        //  Naming as per bison
        typedef size_t location_type;


    // Constructors

        //- Construct null
        genericRagelLemonDriver();

        //- Copy construct
        genericRagelLemonDriver(const genericRagelLemonDriver& rhs) = default;

        //- Move construct
        genericRagelLemonDriver(genericRagelLemonDriver&& rhs) = default;


    //- Destructor
    virtual ~genericRagelLemonDriver() = default;


    // Member Functions

        //- Reset references
        void clear();

        //- Get reference to the input buffer content
        const std::string& content() const
        {
            return content_.get();
        }

        //- Set reference to the input buffer content,
        //- which acts like a std::string_view
        void content
        (
            const std::string& s,
            size_t pos = 0,
            size_t len = std::string::npos
        );

        //- Iterator to begin of content (sub)string
        std::string::const_iterator cbegin() const;

        //- Iterator to end of content (sub)string
        std::string::const_iterator cend() const;

        //- The relative parse position with the content (sub)string
        size_t parsePosition() const
        {
            return position_;
        }

        //- The relative parse position with the content (sub)string
        size_t& parsePosition()
        {
            return position_;
        }

        //- Output the input buffer string content
        Ostream& printBuffer(Ostream& os) const;

        //- Report FatalError
        void reportFatal(const std::string& msg) const;

        //- Report FatalError at parser position
        void reportFatal(const std::string& msg, size_t pos) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace parsing
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
