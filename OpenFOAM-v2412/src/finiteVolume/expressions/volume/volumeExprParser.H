/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::expressions::volumeExpr::parser

Description
    Lemon parser interface for volume expressions grammar

\*---------------------------------------------------------------------------*/

#ifndef expressions_volumeExprParser_H
#define expressions_volumeExprParser_H

#include "exprScanToken.H"
#include "volumeExprFwd.H"

namespace Foam
{
namespace expressions
{
namespace volumeExpr
{

/*---------------------------------------------------------------------------*\
                           Class parser Declaration
\*---------------------------------------------------------------------------*/

class parser
{
    // Private Data

        //- Prompt for parser tracing
        static constexpr const char* const prompt_ = "volExpr:";

        //- The lemon parser (demand-driven)
        void* lemon_;


public:

    //- Local object debugging
    int debug;


    // Constructors

        //- Construct null
        parser() : lemon_(nullptr), debug(volumeExpr::debug) {}


    //- Destructor, deletes parser backend
    ~parser()
    {
        stop();
    }


    // Static Member Functions

        //- Return the text name corresponding to the tokenId
        static word tokenName(int tokenId);

        //- Print all token names
        static void printTokenNames(Ostream& os);

        //- Print all rules
        static void printRules(Ostream& os);


    // Member Functions

        //- Start parsing, with the given driver context
        void start(parseDriver& driver_);

        //- Stop parsing, freeing the allocated parser
        void stop();

        //- Push token type to parser with default token
        void parse(int tokenId);

        //- Push token type/value to parser
        void parse(int tokenId, scanToken tok);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace volumeExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
