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
    Foam::substitutionModels::fileRegEx

Description
    The \c fileRegEx substitution model.

    The example below shows how the keyword \c executionTime is set by
    applying a regular expression (string) to a log file.

    \verbatim
    fileRegEx1
    {
        // Mandatory entries
        type        fileRegEx;
        path        "log.simpleFoam";

        entries
        {
            executionTime    "ExecutionTime = (.*) s  Clock.*";
        }

        // Optional entries
        sectionSeparator    <string>;
        matchSeparator      <string>;
        lastMatch           <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    The entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: functionObjectValue     | word |  yes  | -
      path         | Path to file                       | string |  yes  | -
      entries      | Keyword regular-expression pairs   | dictionary | yes | -
      sectionSeparator | Marker used to separate files into sections  <!--
            --!>   | string | no | "Time ="
      matchSeparator | Separator used to join multiple values <!--
            --!>   | string | no | " "
      lastMatch    | Flag to use last file section      | bool | no    | yes
    \endtable

    The inherited entries are elaborated in:
      - \link substitutionModel.H \endlink

SourceFiles
    fileRegEx.C

---------------------------------------------------------------------------*/

#ifndef Foam_substitutionModels_fileRegEx_H
#define Foam_substitutionModels_fileRegEx_H

#include "substitutionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace substitutionModels
{

/*---------------------------------------------------------------------------*\
                       Class fileRegEx Declaration
\*---------------------------------------------------------------------------*/

class fileRegEx
:
    public substitutionModel
{
    // Private Data

        //- Path to dictionary
        const fileName path_;

        //- Hash table for key and regular expression pairs
        HashTable<string> entries_;

        //- Section separator to dive log files, e.g. into time step info
        const string sectionSeparator_;

        //- Separator to apply between (multiple) matches
        const string matchSeparator_;

        //- Last match wins flag
        bool lastMatch_;


    // Private Functions

        //- No copy construct
        fileRegEx(const fileRegEx&) = delete;

        //- No copy assignment
        void operator=(const fileRegEx&) = delete;


public:

    //- Runtime type information
    TypeName("fileRegEx");


    //- Constructor
    fileRegEx(const dictionary& dict, const Time& time);


    //- Destructor
    virtual ~fileRegEx() = default;


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