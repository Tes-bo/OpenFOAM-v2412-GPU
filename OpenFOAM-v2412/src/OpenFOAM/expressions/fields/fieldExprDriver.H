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
    Foam::expressions::fieldExpr::parseDriver

Description
    Driver for generic primitive field expressions

    In addition to the standard mathematical functions, operations and
    logical and relational operations, the volume expression support the
    following driver-specific functions:

    Functions
    \table
        Function    | Description                      | Number of arguments |
        rand        | Random field                     | 0/1 |
    \endtable

Note
    Use namespace debug switch \c fieldExpr for scanner (2), parser (4)

SourceFiles
    fieldExprDriver.C
    fieldExprDriverFields.C
    fieldExprDriverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_fieldExprDriver_H
#define expressions_fieldExprDriver_H

#include "fieldExprFwd.H"
#include "exprDriver.H"
#include "primitiveFields.H"
#include "genericRagelLemonDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace fieldExpr
{

/*---------------------------------------------------------------------------*\
                         Class parseDriver Declaration
\*---------------------------------------------------------------------------*/

class parseDriver
:
    public parsing::genericRagelLemonDriver,
    public expressions::exprDriver
{
protected:

    // Protected Data

        //- The field size
        label size_;

public:

    ClassName("fieldExpr::driver");

    // Generated Methods

        // No copy copy construct
        parseDriver(const parseDriver&) = delete;

        // No copy assignment
        void operator=(const parseDriver&) = delete;


    // Constructors

        //- Default construct (size=1), or with specified size
        explicit parseDriver(const label len = 1);

        //- Construct for specified size with given dictionary
        parseDriver(const label len, const dictionary& dict);

        //- Construct for specified size, copy of driver context
        parseDriver
        (
            const label len,
            const parseDriver& rhs,
            const dictionary& dict
        );


    //- Destructor
    virtual ~parseDriver() = default;


    // Public Member Functions

        //- The natural field size for the expression
        virtual label size() const
        {
            return size_;
        }

        //- The point field size for the expression
        virtual label pointSize() const
        {
            return size_;
        }


    // Evaluation

        //- Perform parsing on (sub) string
        using genericRagelLemonDriver::content;

        //- Execute the parser
        //  The return value currently has no meaning.
        virtual unsigned parse
        (
            const std::string& expr,
            size_t pos = 0,
            size_t len = std::string::npos
        );


    // Fields

        //- Set result
        template<class Type>
        void setResult(Field<Type>* ptr, bool pointVal = false)
        {
            result().setResult<Type>(ptr, pointVal);
        }


    // New Fields

        //- Return named field (variable) if available
        template<class Type>
        tmp<Field<Type>> getField(const word& fieldName) const;


    // Custom Field Functions

        //- A uniform random field
        tmp<scalarField> field_rand(label seed=0, bool gaussian=false) const;

        //- A Gaussian random field
        tmp<scalarField> field_randGaussian(label seed=0) const
        {
            return field_rand(seed, true);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldExprDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
