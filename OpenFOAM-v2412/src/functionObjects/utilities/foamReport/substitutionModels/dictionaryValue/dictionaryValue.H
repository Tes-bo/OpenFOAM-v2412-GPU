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
    Foam::substitutionModels::dictionaryValue

Description
    The \c dictionaryValue substitution model. Dictionaries can be retrieved
    from an object registry, e.g. time, mesh, or from file.

    The example below shows how the keywords \c p_solver and \c u_solver are set
    by retrieving values from the \c fvSolution dictionary.

    \verbatim
    dictionaryValues1
    {
        // Mandatory entries
        type        dictionaryValue;

        entries
        {
            p_solver    "solvers/p/solver";
            u_solver    "solvers/u/solver";
        }

        // Conditional entries

            // Option-1
            object      "fvSolution";  // registry-based retrieval
            // region      "fluidMesh";

            // Option-2
            // path        "<system>/fvSolution"; // file-based retrieval


        // Optional entries
        separator       <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    The entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: dictionaryValue         | word |  yes  | -
      entries      | keyword lookup pairs               | dictionary | yes | -
      object       | Name of registered dictionary      | string | no  | -
      region       | Name of mesh region                | word | no  | region0
      path         | Path to dictionary file            | string | no  | -
      separator | Sep. when lookup value has multiple tokens | word | no | -
    \endtable

    The inherited entries are elaborated in:
      - \link substitutionModel.H \endlink

SourceFiles
    dictionaryValue.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModels_dictionaryValue_H
#define Foam_substitutionModels_dictionaryValue_H

#include "substitutionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace substitutionModels
{

/*---------------------------------------------------------------------------*\
                       Class dictionaryValue Declaration
\*---------------------------------------------------------------------------*/

class dictionaryValue
:
    public substitutionModel
{
    // Private Data

        //- Dictionary name for registry-based lookup
        word object_;

        //- Region name for registry-based lookup
        word region_;

        //- Path to dictionary for file-based lookup
        fileName path_;

        //- Separator when lookup value has multiple tokens
        const word separator_;

        //- Hash table for key and entry-lookup pairs
        HashTable<string> entries_;


    // Private Functions

        //- No copy construct
        dictionaryValue(const dictionaryValue&) = delete;

        //- No copy assignment
        void operator=(const dictionaryValue&) = delete;


protected:

    // Protected Member Functions

        //- Main function to process the dictionary
        bool processDict
        (
            const dictionary& dict,
            const word& key,
            string& buffer
        ) const;


public:

    //- Runtime type information
    TypeName("dictionaryValue");


    //- Constructor
    dictionaryValue(const dictionary& dict, const Time& time);


    //- Destructor
    virtual ~dictionaryValue() = default;


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