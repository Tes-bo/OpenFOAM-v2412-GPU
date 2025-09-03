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
    Foam::expressions::fieldExpr::scanner

Description
    Ragel lexer/scanner interface for patch expressions.

Note
    Ragel code generated with the ./createCode script.

\*---------------------------------------------------------------------------*/

#ifndef expressions_fieldExprScanner_H
#define expressions_fieldExprScanner_H

#include "fieldExprFwd.H"

namespace Foam
{
namespace expressions
{
namespace fieldExpr
{

/*---------------------------------------------------------------------------*\
                         Class scanner Declaration
\*---------------------------------------------------------------------------*/

class scanner
{
    // Private Data

        //- Wrapped lemon parser
        parser* parser_;

        //- Ragel code state, action
        int cs, act;


    // Private Member Functions

        //- Dispatch .method to parser (if known) or Fatal
        bool dispatch_method
        (
            const parseDriver& driver_,
            word ident  // Receives a copy
        ) const;

        //- Dispatch identifier to parser (if possible) or Fatal
        bool dispatch_ident
        (
            const parseDriver& driver_,
            word ident  // Receives a copy
        ) const;


public:

    //- Local debugging
    int debug;


    // Constructors

        //- Default construct, optionally setting debugging
        explicit scanner(bool withDebug = false)
        :
            parser_(nullptr),
            debug(fieldExpr::debug)
        {
            if (withDebug)
            {
                debug |= 4;
            }
        }


    //- Destructor, deletes parser
    ~scanner();


    // Member Functions

        //- Evaluate sub-string
        bool process
        (
            const std::string& str, size_t pos, size_t len,
            parseDriver& driver_
        );

        //- Evaluate sub-string
        bool process
        (
            const std::string& str, size_t pos,
            parseDriver& driver_
        )
        {
            return process(str, pos, std::string::npos, driver_);
        }

        //- Evaluate string
        bool process(const std::string& str, parseDriver& driver_)
        {
            return process(str, 0, std::string::npos, driver_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
