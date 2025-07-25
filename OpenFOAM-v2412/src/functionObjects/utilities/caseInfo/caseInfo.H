/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::functionObjects::caseInfo

Description
    Collects and writes case information to file.

    Example of function object specification:
    \verbatim
    caseInfo
    {
        type            caseInfo;
        libs            (utilityFunctionObjects);

        // Warn when entries are not found
        lookupMode  warn; // none | warn | error;

        // Write format
        writeFormat json; //  dictionary | json;

        dictionaries
        {
            USolver // User-specified names
            {
                // Look up using registered name
                name        "fvSolution";

                // Optionally limit to specific entries
                include
                (
                    "solvers/U/solver"
                );
            }
            fvSchemes
            {
                name        "fvSchemes";

                // include all entries by default
            }
            timeScheme
            {
                name        "fvSchemes";

                include
                (
                    "/ddtSchemes/default"
                );
            }

            turbulence
            {
                name        "turbulenceProperties";

                // include all entries by default
            }
            controlDict
            {
                // Look up using file path
                path        "<case>/system/controlDict";

                include
                (
                    "application"
                    "deltaT"
                    "startTime"
                    "endTime"
                );
            }
        }

        functionObjects (minMax1); // v2306 only
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description                         | Required | Default
        type          | Type name: caseInfo                 | yes |
        lookupMode    | Lookup mode                         | no  | warn
        writeFormat   | Write format                        | yes |
        dictionaries  | Dictionaries to process             | no  | \<none\>
        functionObjects | Function objects to process       | no  | \<none\>
    \endtable

See also
    Foam::functionObject
    Foam::timeFunctionObject

SourceFiles
    caseInfo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_caseInfo_H
#define functionObjects_caseInfo_H

#include "IOdictionary.H"
#include "writeFile.H"
#include "stateFunctionObject.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class caseInfo Declaration
\*---------------------------------------------------------------------------*/

class caseInfo
:
    public IOdictionary,
    public stateFunctionObject,
    public writeFile
{
public:

    // Public enumerations

        //- Write format enumeration
        enum class writeFormat
        {
            dict,
            json
        };

        //- Lookup mode enumeration
        enum class lookupMode
        {
            none,
            warn,
            error
        };


private:

    // Private Member Data

        //- Write format names
        static const Enum<writeFormat> writeFormatNames_;

        //- Lookup mode names
        static const Enum<lookupMode> lookupModeNames_;

        //- Write/output format, e.g. dictionary, JSON
        writeFormat writeFormat_;

        //- Lookup mode when reading entries
        lookupMode lookupMode_;

        //- Dictionaries
        dictionary dictionaries_;

        //- List of function objects to process
        wordList functionObjectNames_;


    // Private Member Functions

        //- Report
        void report(const string& str) const;

        //- Process dictionary
        void processDict
        (
            dictionary& dict,
            const dictionary& targetDict,
            const entry* includePtr,
            const entry* excludePtr
        ) const;


protected:

    // Protected Member Functions

        //- No copy construct
        caseInfo(const caseInfo&) = delete;

        //- No copy assignment
        void operator=(const caseInfo&) = delete;


    // Write data

        //- Write case meta data
        void writeMeta(dictionary& dict) const;

        //- Write registered dictionaries
        void writeRegisteredDicts
        (
            const objectRegistry& obr,
            dictionary& dict,
            dictionary& dictionaries
        ) const;

        //- Write file-based dictionaries
        void writeFileDicts
        (
            dictionary& dict,
            dictionary& dictionaries
        ) const;

        //- Write function object results
        void writeFunctionObjects(dictionary& dict) const;

        //- Write mesh statistics
        void writeMeshStats(const polyMesh& mesh, dictionary& dict) const;

        //- Write mesh patches
        void writePatches(const fvMesh& mesh, dictionary& dict) const;


public:

    //- Runtime type information
    TypeName("caseInfo");


    // Constructors

        //- Construct from Time and dictionary
        caseInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~caseInfo() = default;


    // Member Functions

        using regIOobject::read;
        using regIOobject::write;

        //- Read the controls
        virtual bool read(const dictionary& dict);

        //- Execute, does nothing
        virtual bool execute();

        //- Write the caseInfo
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //