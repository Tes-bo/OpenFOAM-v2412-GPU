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
    Foam::substitutionModels::environmentVariable

Description
    The \c environmentVariable substitution model.

    \verbatim
    environmentVariables1
    {
        // Mandatory entries
        type        environmentVariable;

        entries
        {
            home    "HOME";
            ldpath  "LD_LIBRARY_PATH";
        }

        // Inherited entries
        ...
    }
    \endverbatim

    The entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: environmentVariable     | word |  yes  | -
      entries      | Keyword lookup pairs               | dictionary | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link substitutionModel.H \endlink

SourceFiles
    environmentVariable.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModels_environmentVariable_H
#define Foam_substitutionModels_environmentVariable_H

#include "substitutionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace substitutionModels
{

/*---------------------------------------------------------------------------*\
                  Class environmentVariable Declaration
\*---------------------------------------------------------------------------*/

class environmentVariable
:
    public substitutionModel
{
    // Private Data

        //- Hash table for key and environment variable pairs
        HashTable<string> entries_;


    // Private Functions

        //- No copy construct
        environmentVariable(const environmentVariable&) = delete;

        //- No copy assignment
        void operator=(const environmentVariable&) = delete;


public:

    //- Runtime type information
    TypeName("environmentVariable");


    //- Constructor
    environmentVariable(const dictionary& dict, const Time& time);


    //- Destructor
    virtual ~environmentVariable() = default;


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