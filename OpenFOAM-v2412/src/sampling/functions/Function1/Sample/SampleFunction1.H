/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::Function1Types::sample

Description
    Minimal example by using \c system/controlDict.functions:
    \verbatim

    \verbatim
        <entryName> sample;
        <entryName>Coeffs
        {
            field           <field name>;
            position        (0 0 0);
            interpolationScheme cell;
        }
    \endverbatim

SourceFiles
    sampleFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_Sample_H
#define Foam_Function1Types_Sample_H

#include "Function1.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Sample Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Sample
:
    public Function1<Type>
{
    // Private Data

        //- Name of field to sample
        const word fieldName_;

        //- Sample position
        const point position_;

        //- Interpolation scheme name; default = "cell"
        const word interpolationScheme_;

        //- Sample cell
        mutable label celli_;

        //- Points event number used to determine if celli should be updated
        mutable label pointEventNo_;


    // Private Member Functions

        //- Set the sample celli; error on not found
        void setSampleCell() const;

        //- No copy assignment
        void operator=(const Sample<Type>&) = delete;


public:

    //- Runtime type information
    TypeName("sample");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Sample
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from components
        Sample
        (
            const word& entryName,
            const List<Tuple2<Type, Type>>& coeffs,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy constructor
        explicit Sample(const Sample& poly);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~Sample() = default;


    // Member Functions

        //- Return Sample value
        virtual Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Write as primitive (inline) format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SampleFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
