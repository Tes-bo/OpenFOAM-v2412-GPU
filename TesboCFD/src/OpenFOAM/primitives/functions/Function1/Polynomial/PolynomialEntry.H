/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::Function1Types::PolynomialEntry

Description
    PolynomialEntry container data entry for scalars. Items are stored in a
    list of Tuple2's. Data is input in the form,
    e.g. for an entry \<entryName\> that describes y = x^2 + 2x^3

    Inline specification:
    \verbatim
        <entryName>  polynomial
        (
            (1  2)
            (2  3)
        );
    \endverbatim

    Dictionary format:
    \verbatim
        <entryName>
        {
            type   polynomial;
            coeffs
            (
                (1 2)
                (2 3)
            );
        }
    \endverbatim

SourceFiles
    PolynomialEntry.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Polynomial_H
#define Foam_Function1Types_Polynomial_H

#include "Function1.H"
#include "Tuple2.H"
#include "Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Polynomial
:
    public Function1<Type>
{
    // Private Data

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<Type, Type>> coeffs_;

        //- Flag to indicate whether polynomial can be integrated
        bool canIntegrate_;


    // Private Member Functions

        //- Check coefficients and if polynomial can be integrated
        void checkCoefficients();

        //- No copy assignment
        void operator=(const Polynomial<Type>&) = delete;


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Polynomial
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from components
        Polynomial
        (
            const word& entryName,
            const List<Tuple2<Type, Type>>& coeffs,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy constructor
        explicit Polynomial(const Polynomial& poly);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Polynomial() = default;


    // Member Functions

        //- Convert time
        virtual void userTimeToTime(const Time& t);

        //- Return Polynomial value
        virtual Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PolynomialEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
