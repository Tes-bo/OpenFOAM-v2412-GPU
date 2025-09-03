/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::polynomialFunction

Description
    Polynomial function representation

    \verbatim
        poly = logCoeff*log(x) + sum(coeffs[i]*x^i)
    \endverbatim

    where <tt> 0 <= i <= N </tt>

    - integer powers, starting at zero
    - \c value(x) to evaluate the poly for a given value
    - \c integrate(x1, x2) between two scalar values
    - \c integral() to return a new, integral coeff polynomial
      - increases the size (order)
    - \c integralMinus1() to return a new, integral coeff polynomial where
      the base poly starts at order -1

See also
    Foam::Polynomial for a templated implementation

SourceFiles
    polynomialFunction.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialFunction_H
#define polynomialFunction_H

#include "scalarList.H"
#include "Ostream.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polynomialFunction;

Istream& operator>>(Istream&, polynomialFunction& poly);
Ostream& operator<<(Ostream&, const polynomialFunction& poly);

/*---------------------------------------------------------------------------*\
                     Class polynomialFunction Declaration
\*---------------------------------------------------------------------------*/

class polynomialFunction
:
    private scalarList
{
    // Private Data

        //- Include the log term? - only activated using integralMinus1()
        bool logActive_;

        //- Log coefficient - only activated using integralMinus1()
        scalar logCoeff_;


    // Private Member Functions

        //- Return integral coefficients.
        //  Argument becomes zeroth element (constant of integration)
        static polynomialFunction cloneIntegral
        (
            const polynomialFunction&,
            const scalar intConstant = 0
        );

        //- Return integral coefficients when lowest order is -1.
        //  Argument becomes zeroth element (constant of integration)
        static polynomialFunction cloneIntegralMinus1
        (
            const polynomialFunction&,
            const scalar intConstant = 0
        );

        //- Check size is non-zero or trigger FatalErrot
        void checkSize() const;


public:

    //- Runtime type information
    TypeName("polynomialFunction");


    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Default construct as size 1 with coefficient == 0
        polynomialFunction();

        //- Construct a particular size, with all coefficients == 0
        explicit polynomialFunction(const label order);

        //- Construct from an initializer list of coefficients
        explicit polynomialFunction(std::initializer_list<scalar> coeffs);

        //- Construct from a list of coefficients
        explicit polynomialFunction(const UList<scalar>& coeffs);

        //- Construct from Istream
        explicit polynomialFunction(Istream& is);


    //- Destructor
    virtual ~polynomialFunction() = default;


    // Member Functions

        //- Non-zero number of coefficients
        using scalarList::empty;

        //- The number of coefficients
        using scalarList::size;

        //- Return coefficient at given index
        using scalarList::operator[];


        // Access

            //- True if the log term is active
            bool logActive() const;

            //- The log coefficient
            scalar logCoeff() const;


        // Evaluation

            //- Return polynomial value
            scalar value(const scalar x) const;

            //- Integrate between two values
            scalar integrate(const scalar x1, const scalar x2) const;


            //- Return integral coefficients.
            //  Argument becomes zeroth element (constant of integration)
            polynomialFunction integral
            (
                const scalar intConstant = 0
            ) const;

            //- Return integral coefficients when lowest order is -1.
            //  Argument becomes zeroth element (constant of integration)
            polynomialFunction integralMinus1
            (
                const scalar intConstant = 0
            ) const;


    // Member Operators

        //- Equality of coefficients, and logCoeff (if active)
        bool operator==(const polynomialFunction& rhs) const;

        bool operator!=(const polynomialFunction& rhs) const
        {
            return !operator==(rhs);
        }

        polynomialFunction& operator+=(const polynomialFunction&);
        polynomialFunction& operator-=(const polynomialFunction&);

        polynomialFunction& operator*=(const scalar);
        polynomialFunction& operator/=(const scalar);


    // IOstream Operators

    friend Istream& operator>>(Istream&, polynomialFunction& poly);
    friend Ostream& operator<<(Ostream&, const polynomialFunction& poly);
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

polynomialFunction operator+
(
    const polynomialFunction&,
    const polynomialFunction&
);


polynomialFunction operator-
(
    const polynomialFunction&,
    const polynomialFunction&
);


polynomialFunction operator*
(
    const scalar,
    const polynomialFunction&
);


polynomialFunction operator/
(
    const scalar,
    const polynomialFunction&
);


polynomialFunction operator*
(
    const polynomialFunction&,
    const scalar
);


polynomialFunction operator/
(
    const polynomialFunction&,
    const scalar
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
