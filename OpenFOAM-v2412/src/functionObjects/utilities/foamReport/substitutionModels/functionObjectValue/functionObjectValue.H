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
    Foam::substitutionModels::functionObjectValue

Description
    functionObjectValue substitution model.

Usage
    The \c functionObjectValue substitution model.

    The example below shows how the keywords \c cont_error_* are set by
    retrieving the values \c local, \c global, \c cumulative from the function
    object \c functionObjectValue.

    \verbatim
    functionObjectValue1
    {
        // Mandatory entries
        type        functionObjectValue;
        functionObject continuityError1;

        entries
        {
            cont_error_local    local;
            cont_error_global   global;
            cont_error_cumulative cumulative;
        }

        // Optional entries
        debugValues   <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    The entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: functionObjectValue     | word |  yes  | -
      functionObject | Name of function object to query | word |  yes  | -
      entries      | Keyword-lookup pairs               | dictionary | yes | -
      debugValues  | Flag to show available function values | bool | no | false
    \endtable

    The inherited entries are elaborated in:
      - \link substitutionModel.H \endlink

SourceFiles
    functionObjectValue.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModels_functionObjectValue_H
#define Foam_substitutionModels_functionObjectValue_H

#include "substitutionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace substitutionModels
{

/*---------------------------------------------------------------------------*\
                       Class functionObjectValue Declaration
\*---------------------------------------------------------------------------*/

class functionObjectValue
:
    public substitutionModel
{
    // Private Data

        //- Name of function object
        const word functionObject_;

        //- Hash table for key and entry-lookup pairs
        HashTable<word> entries_;

        //- Debug - shows available function values
        bool debugValues_;


    // Private Functions

        //- Get the result value from the function object
        template<class Type>
        bool getValue(OStringStream& oss, const word& lookup) const;

        //- No copy construct
        functionObjectValue(const functionObjectValue&) = delete;

        //- No copy assignment
        void operator=(const functionObjectValue&) = delete;


public:

    //- Runtime type information
    TypeName("functionObjectValue");


    //- Constructor
    functionObjectValue
    (
        const dictionary& dict,
        const Time& time
    );


    //- Destructor
    virtual ~functionObjectValue() = default;


    // Member Functions

        //- Update model local data
        virtual bool update();

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