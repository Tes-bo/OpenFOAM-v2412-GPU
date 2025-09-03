/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
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
    Foam::expressions::exprResultStored

Description
    An exprResult with persistence.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        name        | The result name                       | yes |
        initialValue | The initial value                    | yes |
        value        | The expression result (dictionary)   | no  |
    \endtable

SourceFiles
    StoredExpressionResult.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResultStored_H
#define Foam_expressions_exprResultStored_H

#include "exprResult.H"
#include "exprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprResultStored Declaration
\*---------------------------------------------------------------------------*/

class exprResultStored
:
    public expressions::exprResult
{
    // Private Data

        //- The name of the expression
        word name_;

        //- The initial value expression
        expressions::exprString startExpr_;


protected:

    // Protected Member Functions

        //- Reset at new timestep - disabling locally
        virtual void resetImpl() {}


public:

    //- Runtime type information
    TypeNameNoDebug("exprResultStored");


    // Constructors

        //- Default construct
        exprResultStored();

        //- Copy construct
        exprResultStored(const exprResultStored& rhs);

        //- Construct from a dictionary
        explicit exprResultStored(const dictionary& dict);


    // Selectors

        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>
            (
                new exprResultStored(*this)
            );
        }


    //- Destructor
    virtual ~exprResultStored() = default;


    // Member Functions

        //- The name of the expression
        const word& name() const
        {
            return name_;
        }

        //- The initial value expression
        const expressions::exprString& initialValueExpression() const
        {
            return startExpr_;
        }


        void writeDict(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const exprResultStored& rhs);

        //- Copy assignment
        void operator=(const exprResult& rhs);

        //- Compare for inequality (need by Ulist)
        bool operator!=(const exprResultStored& rhs) const noexcept
        {
            return (this != &rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IO Operators
Istream& operator>>(Istream& is, expressions::exprResultStored& data);
Ostream& operator<<(Ostream& is, const expressions::exprResultStored& data);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
