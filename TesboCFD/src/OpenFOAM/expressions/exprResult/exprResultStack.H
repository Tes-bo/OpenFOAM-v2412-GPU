/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::expressions::exprResultStack

Description
    A stack of polymorphic fields.
    Can be used to build a list of results one at a time.

SourceFiles
    exprResultStack.C
    exprResultStackTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResultStack_H
#define Foam_expressions_exprResultStack_H

#include "exprResult.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                       Class exprResultStack Declaration
\*---------------------------------------------------------------------------*/

class exprResultStack
:
    public expressions::exprResult
{
    // Private Data

        //- Type-checked pop value.
        //  \return True if the type check was satisfied
        template<class T>
        bool popChecked(exprResult& result);

        //- Type-checked push value.
        //  \return True if the type check was satisfied
        template<class T>
        bool pushChecked(const exprResult& result);


public:

        //- Runtime type information
        TypeName("exprResultStack");


    // Constructors

        //- Default construct
        exprResultStack();

        //- Copy construct
        exprResultStack(const exprResultStack& rhs);

        //- Construct from a dictionary
        explicit exprResultStack(const dictionary& dict);


    // Selectors

        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>
            (
                new exprResultStack(*this)
            );
        }


    //- Destructor
    virtual ~exprResultStack() = default;


    // Member Functions

        //- Pop the last value as an expression result
        exprResult pop();

        //- Push an expression result value
        void push(const exprResult& result);


    // Member Operators

        //- Copy assignment
        void operator=(const exprResultStack& rhs);

        //- Copy assignment
        void operator=(const exprResult& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "exprResultStackTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
