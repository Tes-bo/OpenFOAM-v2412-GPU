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
    Foam::substitutionModel

Description
    Base class for substitution models.

    Provides a static hash table for builtin keyword-value pairs and functions
    to manipulate/interact.

SourceFiles
    substitutionModel.C
    substitutionModelNew.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModel_H
#define Foam_substitutionModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class substitutionModel Declaration
\*---------------------------------------------------------------------------*/

class substitutionModel
{
public:

    // Static Data Members

        //- Keyword starting characters
        static const word KEY_BEGIN;

        //- Keyword ending characters
        static const word KEY_END;

        //- Built-in substitutions
        static HashTable<string> builtin_;


    // Static Member Functions

        //- Return a key representation from a word
        static string keyify(const word& w);

        //- Clean the key text
        static word cleanKey(const string& str);

        //- Return all keys from a string buffer
        //  Also cleans the key strings in the buffer
        static wordList getKeys(string& buffer);

        //- Add a builtin to the hash table - does not overwrite
        static void addBuiltinStr(const word& key, const string& value);

        //- Add a builtin to the hash table - does not overwrite
        template<class Type>
        static void addBuiltin(const word& key, const Type& value);

        //- Return true if key is builtin
        static bool containsBuiltin(const word& key);

        //- Set a builtin to the hash table
        static void setBuiltinStr(const word& key, const string& value);

        //- Set a builtin to the hash table
        template<class Type>
        static void setBuiltin(const word& key, const Type& value);

        //- Replace key in string
        static bool replaceBuiltin(const word& key, string& str);

        //- Replace all occurrences of key in string
        static bool replaceBuiltin(string& str);

        //- Write all builtins to stream
        static void writeBuiltins(Ostream& os);


private:

    // Private Functions

        //- No copy construct
        substitutionModel(const substitutionModel&) = delete;

        //- No copy assignment
        void operator=(const substitutionModel&) = delete;


protected:

    // Protected Data

        //- Construction dictionary
        const dictionary dict_;

        //- Reference to the time database
        const Time& time_;


public:

    //- Runtime type information
    TypeName("substitutionModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            substitutionModel,
            dictionary,
            (
                const dictionary& dict,
                const Time& time
            ),
            (dict, time)
        );


    // Selectors

        //- Return a reference to the selected substitution model
        static autoPtr<substitutionModel> New
        (
            const dictionary& dict,
            const Time& time
        );


    //- Constructor
    substitutionModel
    (
        const dictionary& dict,
        const Time& time
    );


    //- Destructor
    virtual ~substitutionModel() = default;


    // Member Functions

        //- Update model local data
        virtual bool update() { return true; }

        //- Return true of model applies to this keyName
        virtual bool valid(const word& keyName) const = 0;

        //- Apply substitutions to this string buffer
        virtual bool apply(const word& key, string& buffer) const = 0;

        //- Return a word list of the keys
        virtual wordList keys() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "substitutionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //