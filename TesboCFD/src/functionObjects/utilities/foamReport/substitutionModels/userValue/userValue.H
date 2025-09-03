/*---------------------------------------------------------------------------*\
  =========                |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
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
    Foam::substitutionModels::userValue

Description
    The \c userValue substitution model. Dictionaries can be retrieved from
    an object registry, e.g. time, mesh, or from a file.

    \verbatim
    userValues1
    {
        // Mandatory entries
        type        userValue;

        entries
        {
            my_keyword1 "My local string value";
            my_keyword2 "My local string value";
        }

        // Inherited entries
        ...
    }
    \endverbatim

    The entries mean:
    \table
      Property     | Description                    | Type | Reqd  | Deflt
      type         | Type name: userValue           | word |  yes  | -
      entries      | Keyword lookup pairs           | dictionary | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link substitutionModel.H \endlink

SourceFiles
    userValue.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModels_userValue_H
#define Foam_substitutionModels_userValue_H

#include "substitutionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace substitutionModels
{

/*---------------------------------------------------------------------------*\
                          Class userValue Declaration
\*---------------------------------------------------------------------------*/

class userValue
:
    public substitutionModel
{
    // Private Data

        //- Hash table for key and environment variable pairs
        HashTable<string> entries_;


    // Private Functions

        //- No copy construct
        userValue(const userValue&) = delete;

        //- No copy assignment
        void operator=(const userValue&) = delete;


public:

    //- Runtime type information
    TypeName("userValue");


    //- Constructor
    userValue(const dictionary& dict, const Time& time);


    //- Destructor
    virtual ~userValue() = default;


    // Member Functions

        //- Return true of model applies to this keyName
        virtual bool valid(const word& keyName) const;

        //- Apply substitutions to this string buffer
        virtual bool apply(const word& key, string& buffer) const;

        //- Return a word list of the keys
        virtual wordList keys() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace substitutionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //