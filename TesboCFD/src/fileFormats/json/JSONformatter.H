/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::JSONformatter

Description
    An wrapper for Ostream that outputs content in JSON format.

SourceFiles
    JSONformatter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_JSONformatter_H
#define Foam_JSONformatter_H

#include "Ostream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class JSONformatter Declaration
\*---------------------------------------------------------------------------*/

class JSONformatter
{
    // Private Data

        //- Reference to the output stream
        Ostream& os_;


    // Private Member Functions

        //- Write token to stream
        bool writeToken(const token& t);


public:

    //- Declare type-name
    TypeName("JSONformatter");

    //- Construct from Ostream
    explicit JSONformatter(Ostream& os);

    //- Destructor
    virtual ~JSONformatter() = default;


    // Member Functions

        // Write Functions

        virtual Ostream& write(const bool val);
        virtual Ostream& write(const int32_t val);
        virtual Ostream& write(const int64_t val);
        virtual Ostream& write(const float val);
        virtual Ostream& write(const double val);
        virtual Ostream& write(const word& str);
        virtual Ostream& write(const std::string& str);
        virtual Ostream& write(const char* str);

        //- Write JSON keyword
        virtual Ostream& writeKeyword(const keyType& keyword);

        //- Write OpenFOAM dictionary to JSON dictionary
        virtual Ostream& writeDict(const dictionary& dict);

        //- Write JSON keyword-value pair (for primitive types)
        template<class Type>
        Ostream& writeEntry(const word& keyword, const Type& val)
        {
            os_.writeQuoted(keyword) << " : ";
            write(val);
            return os_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
