/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::expressions::scanToken

Description
    A low-level input/scan token content.
    No defined constructors/destructors.
    All memory management is manual!

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_scanToken_H
#define Foam_expressions_scanToken_H

#include "scalar.H"
#include "vector.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                          Class scanToken Declaration
\*---------------------------------------------------------------------------*/

struct scanToken
{
    //- Tagged union types
    enum tokenType : unsigned char
    {
        LABEL = 0,
        SCALAR,
        VECTOR,
        WORD
    };


    // Data

        //- The data content (as a union).
        //  For memory alignment have this appear as the first member.
        union
        {
            Foam::label    labelValue;
            Foam::scalar   scalarValue;
            Foam::vector*  vectorPtr;
            Foam::word*    wordPtr;
            Foam::word*    name_;
        };

        //- The token type
        tokenType type_;


    // Member Functions

        //- Return a null token (label = 0) - in lieu of a default constructor
        static scanToken null();

        //- Assign type/value to be LABEL. Does not call destroy().
        void setLabel(label val);

        //- Assign type/value to be SCALAR. Does not call destroy().
        void setScalar(scalar val);

        //- Assign type/value to be VECTOR. Does not call destroy().
        void setVector(scalar x, scalar y, scalar z);

        //- Assign type/value to be VECTOR. Does not call destroy().
        void setVector(const vector& val);

        //- Assign type/value to be WORD (name). Does not call destroy().
        void setWord(const word& val);

        //- True if a pointer type
        bool is_pointer() const noexcept
        {
            return (type_ == tokenType::VECTOR || type_ == tokenType::WORD);
        }

        //- Manual deletion of pointer types
        void destroy();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
