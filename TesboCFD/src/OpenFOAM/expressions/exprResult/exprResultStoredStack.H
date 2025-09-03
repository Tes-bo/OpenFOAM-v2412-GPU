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
    Foam::expressions::exprResultStoredStack

Description
    An exprResultStack with persistence.

SourceFiles
    exprResultStoredStack.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResultStoredStack_H
#define Foam_expressions_exprResultStoredStack_H

#include "exprResultStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                    Class exprResultStoredStack Declaration
\*---------------------------------------------------------------------------*/

class exprResultStoredStack
:
    public expressions::exprResultStack
{
protected:

    // Protected Member Functions

        //- Reset at new timestep - disabling locally
        virtual void resetImpl() {}


public:

    //- Runtime type information
    TypeNameNoDebug("exprResultStoredStack");


    // Constructors

        //- Default construct
        exprResultStoredStack();

        //- Copy construct
        exprResultStoredStack(const exprResultStoredStack& rhs);

        //- Construct from a dictionary
        explicit exprResultStoredStack(const dictionary& dict);


    // Selectors

        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>
            (
                new exprResultStoredStack(*this)
            );
        }


    //- Destructor
    virtual ~exprResultStoredStack() = default;


    // Member Operators

        //- Copy assignment
        void operator=(const exprResultStoredStack& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
