/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::namedDictionary

Description
    A tuple of keyType and dictionary, which can be used when reading
    named or unnamed dictionary entries or simply a name.

    For example,
    \verbatim
    fields
    (
        U
        T { relax  false; }
    );
    \endverbatim

    The namedDictionary can also be used in situations where an individual
    dictionary entry should be read. The keyword() and dict() methods
    correspond to Foam::entry naming conventions.
    \verbatim
    actions
    (
        testing { action new; ... }   // An action with a name
        { action subset;  }           // Unnamed action
    );
    \endverbatim
    Normal dictionary reading would fail for this type of input since the
    leading 'testing' keyword would cause the entire content to be considered
    a single dictionary.

Note
    No distinction currently made between a missing and an empty dictionary.

SourceFiles
    namedDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_namedDictionary_H
#define Foam_namedDictionary_H

#include "dictionary.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class namedDictionary;

Istream& operator>>(Istream&, namedDictionary&);
Ostream& operator<<(Ostream&, const namedDictionary&);

/*---------------------------------------------------------------------------*\
                       Class namedDictionary Declaration
\*---------------------------------------------------------------------------*/

class namedDictionary
:
    public Tuple2<keyType, dictionary>
{
public:

    // Constructors

        //- Inherit constructors
        using Tuple2<keyType, dictionary>::Tuple2;

        //- Default construct
        namedDictionary() = default;

        //- Construct from Istream
        explicit namedDictionary(Istream& is);


    //- Destructor
    ~namedDictionary() = default;


    // Member Functions

        //- Clear keyword and dictionary
        void clear();

        //- Empty if both keyword and dictionary are empty
        bool empty() const noexcept;

        //- Return keyword
        const keyType& keyword() const noexcept { return first(); }

        //- Return non-const access to keyword
        keyType& keyword() noexcept { return first(); }

        //- Read-access to the dictionay
        const dictionary& dict() const noexcept { return second(); }

        //- Write access to the dictionay
        dictionary& dict() noexcept { return second(); }


    // IOstream Operators

        friend Istream& operator>>(Istream&, namedDictionary&);
        friend Ostream& operator<<(Ostream&, const namedDictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
