/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Namespace
    Foam::FlatOutput

Description
    Various output adaptors, principally to output a list of items
    on a single line.

\*---------------------------------------------------------------------------*/

#ifndef FlatOutput_H
#define FlatOutput_H

#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FlatOutput
{

// Forward Declarations
template<class Container, class Delimiters> class OutputAdaptor;

} // End namespace FlatOutput


// Forward Declarations
template<class Container, class Delimiters>
inline Ostream& operator<<
(
    Ostream& os,
    const FlatOutput::OutputAdaptor<Container, Delimiters>& adaptor
);


namespace FlatOutput
{

/*---------------------------------------------------------------------------*\
                         Class Decorators Declaration
\*---------------------------------------------------------------------------*/

//- List decorators with \c open, \c close and \c separator characters
template<char OpenChar, char CloseChar, char SepChar>
struct Decorators
{
    static constexpr char open  = OpenChar;
    static constexpr char close = CloseChar;
    static constexpr char separator = SepChar;
};


#undef  makeDecorator
#define makeDecorator(Name, Open, Close, Sep)                                  \
    /*! \brief Surround with \c Open and \c Close separate with \c Sep */      \
    struct Name : public Decorators<Open, Close, Sep> {};

makeDecorator(BareComma, '\0','\0', ',');
makeDecorator(BareSpace, '\0','\0', ' ');

makeDecorator(BraceComma, '{','}', ',');
makeDecorator(BraceSpace, '{','}', ' ');

makeDecorator(ParenComma, '(',')', ',');
makeDecorator(ParenSpace, '(',')', ' ');  // Normal default

makeDecorator(PointyComma, '<','>', ',');
makeDecorator(PointySpace, '<','>', ' ');

makeDecorator(SquareComma,'[',']', ',');
makeDecorator(SquareSpace,'[',']', ' ');

#undef makeDecorator


/*---------------------------------------------------------------------------*\
                        Class OutputAdaptor Declaration
\*---------------------------------------------------------------------------*/

//- An output adaptor with a write method and an Ostream operator.
//
//  Generate single line (flat) output using the characters specified by
//  the templated Delimiters.
//  Normally called with the global flatOutput() function.
//  For example,
//  \code
//
//  /* With default parenthesis/space delimiters */
//  Info<< "Names: " << flatOutput(names) << nl;
//
//  /* Other delimiters */
//  Info<< flatOutput(names, FlatOutput::SquareComma{}) << nl;
//
//  /* User-specified delimiters */
//  Info<< flatOutput(names, FlatOutput::Decorators<'[',')',':'>{}) << nl;
//
//  \endcode
//
template<class Container, class Delimiters>
class OutputAdaptor
{
    // Private Data

        //- The container of values for output
        const Container& values;

public:

    // Constructors

        //- Construct from component
        explicit OutputAdaptor(const Container& obj)
        :
            values(obj)
        {}


    // Member Functions

        //- Write list using \c open, \c close and \c separator characters
        //- specified by Delimiters template, which generally results in
        //- a single line without line breaks.
        //
        //  \note Suppresses nul char output.
        //        No special handling for newline separators.
        inline Ostream& write(Ostream& os) const
        {
            bool started = false;

            // In c++17, can use constexpr if

            if (Delimiters::open)
            {
                os  << Delimiters::open;
            }
            for (const auto& item : values)
            {
                if (started)
                {
                    if (Delimiters::separator)
                    {
                        os  << Delimiters::separator;
                    }
                }
                else
                {
                    started = true;
                }
                os  << item;
            }
            if (Delimiters::close)
            {
                os << Delimiters::close;
            }

            return os;
        }


    // Operators

        //- Ostream Operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const OutputAdaptor<Container, Delimiters>& adaptor
        )
        {
            return adaptor.write(os);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FlatOutput
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Global flatOutput() function with specified output delimiters
template<class Container, class Delimiters>
inline FlatOutput::OutputAdaptor<Container, Delimiters>
flatOutput
(
    const Container& obj,
    Delimiters delim
)
{
    return FlatOutput::OutputAdaptor<Container, Delimiters>(obj);
}


//- Global flatOutput() function with default (parenthesis/space) delimiters
template<class Container>
inline FlatOutput::OutputAdaptor<Container, FlatOutput::ParenSpace>
flatOutput
(
    const Container& obj
)
{
    return FlatOutput::OutputAdaptor<Container, FlatOutput::ParenSpace>(obj);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
