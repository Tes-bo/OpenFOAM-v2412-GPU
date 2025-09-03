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
    Foam::expressions::exprResultDelayed

Description
    An exprResult with an additional delay before evaluation

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        name        | The result name                       | yes |
        delay       | The delay before starting             | yes |
        storeInterval | The storage interval (time)         | yes |
        startupValue | The initial startup value            | yes |
        value        | The expression result (dictionary)   | no  |
        storedValues | The list of stored values            | yes |
    \endtable

SourceFiles
    exprResultDelayed.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResultDelayed_H
#define Foam_expressions_exprResultDelayed_H

#include "exprResult.H"
#include "exprString.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprResultDelayed Declaration
\*---------------------------------------------------------------------------*/

class exprResultDelayed
:
    public expressions::exprResult
{
    // Typedefs

        //- Expression at a given time value.
        typedef Tuple2<scalar, expressions::exprResult> ValueAtTime;


    // Private Data

        //- The name of the expression
        word name_;

        //- The initial value expression
        expressions::exprString startExpr_;

        //- The pending result to be set
        expressions::exprResult settingResult_;

        //- The old results
        DLList<ValueAtTime> storedValues_;

        //- How often values should be stored
        scalar storeInterval_;

        //- The size of the delay
        scalar delay_;


protected:

    // Protected Member Functions

        //- Simulate virtual templated methods
        inline expressions::exprResult& target() { return settingResult_; }


public:

    //- Runtime type information
    TypeNameNoDebug("exprResultDelayed");


    // Constructors

        //- Default construct
        exprResultDelayed();

        //- Copy construct
        exprResultDelayed(const exprResultDelayed&);

        //- Construct from a dictionary
        exprResultDelayed(const dictionary& dict);

        //- Clone
        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>
            (
                new exprResultDelayed(*this)
            );
        }


    //- Destructor
    virtual ~exprResultDelayed() = default;


    // Member Functions

        //- The expression name
        const word& name() const noexcept
        {
            return name_;
        }

        //- The initial value expression
        const expressions::exprString& startupValueExpression() const noexcept
        {
            return startExpr_;
        }

        //- Update the read-value
        //  return true if there was a valid value.
        //  this does not do the work of setReadValue because we have no
        //  access to the Parser
        bool updateReadValue(const scalar& timeVal);

        //- Set the readValue with a calculated value
        void setReadValue(const exprResult& val);

        //- Add a stored value
        void storeValue(const scalar& timeVal);


    // Write

        void writeDict(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const exprResultDelayed& rhs);

        //- Copy assignment
        void operator=(const exprResult& rhs);

        //- Move assignment
        void operator=(exprResult&& rhs);

        //- Compare for inequality (need by Ulist)
        bool operator!=(const exprResultDelayed& rhs) const noexcept
        {
            return (this != &rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// IO Operators
Istream& operator>>(Istream& is, expressions::exprResultDelayed& data);
Ostream& operator<<(Ostream& os, const expressions::exprResultDelayed& data);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
