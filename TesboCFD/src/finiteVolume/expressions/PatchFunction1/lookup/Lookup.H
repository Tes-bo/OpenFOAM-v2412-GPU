/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::Function1Types::Lookup

Description
    Function1 to lookup UniformDimensionedField from an objectregistry.

    The dictionary specification would typically resemble this:
    \verbatim
    entry
    {
        type        lookup;
        name        myScalar;
    }
    \endverbatim

    where the entries mean:
    \table
        Property  | Description             | Type | Reqd | Default
        type      | Function type: lookup   | word | yes |
        name      | name of variable        | word | yes |
    \endtable

SourceFiles
    Lookup.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Types_Lookup_H
#define Function1Types_Lookup_H

#include "Function1.H"
#include "lookupBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Lookup Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Lookup
:
    public lookupBase,
    public Function1<Type>
{
public:

    // Runtime type information
    TypeName("lookup");


    // Generated Methods

        //- No copy assignment
        void operator=(const Lookup<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        Lookup
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit Lookup(const Lookup<Type>& rhs);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Lookup<Type>(*this));
        }


    //- Destructor
    virtual ~Lookup() = default;


    // Member Functions

        //- Return value for time t
        virtual Type value(const scalar t) const;

        //- Integrate between two (scalar) values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        using lookupBase::writeEntries;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Lookup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
