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
    Foam::Function1Types::FunctionObjectTrigger

Description
    Returns a 0/1 value corresponding to function object trigger levels.

    Usage:
    \verbatim
        <entryName> functionObjectTrigger;
        <entryName>Coeffs
        {
            triggers        (1 3 5);
            defaultValue    false;  // Optional
        }
    \endverbatim

    Where:
    \table
        Property | Description                                  | Required
        triggers | List of active trigger states to check for   | yes |
        defaultValue | Treatment for unactivated trigger state  | no  | false
    \endtable

    In some circumstances, it can be useful to treat an unactivated trigger
    as being true. This is the role of the "defaultValue" keyword.

Note
- does not implement integrate()

SourceFiles
    FunctionObjectTrigger.C
    FunctionObjectTriggerI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_FunctionObjectTrigger_H
#define Foam_Function1Types_FunctionObjectTrigger_H

#include "Function1.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                     Class FunctionObjectTrigger Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FunctionObjectTrigger
:
    public Function1<Type>
{
    // Private Data

        //- Trigger indices when it is considered active
        labelList triggers_;

        //- Treatment for unactivated trigger state (true/false)
        bool defaultValue_;


    // Private Member Functions

        //- Is the trigger considered active?
        inline bool active() const;

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    //- Runtime type information
    TypeName("functionObjectTrigger");


    // Generated Methods

        //- No copy assignment
        void operator=(const FunctionObjectTrigger<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        FunctionObjectTrigger
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit FunctionObjectTrigger(const FunctionObjectTrigger<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~FunctionObjectTrigger() = default;


    // Member Functions

        //- Return the trigger indices
        inline const labelList& triggers() const noexcept;

        //- Change the trigger indices
        inline void resetTriggers(const labelUList& indices);


        //- Return 0/1 value at current time
        virtual inline Type value(const scalar /*unused*/) const;

        //- Integrate between two (scalar) values. Not implemented!
        virtual inline Type integrate(const scalar, const scalar) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FunctionObjectTriggerI.H"

#ifdef NoRepository
    #include "FunctionObjectTrigger.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
