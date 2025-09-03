/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Foam::functionObjectList

Description
    List of function objects with start(), execute() and end() functions
    that is called for each object.

    \verbatim
    functions   // sub-dictionary name under the system/controlDict file
    {
        ..optional entries..

        <userDict1>
        {
            // Mandatory entries
            type    <typeName>;
            libs    (<libName> .. <libName>);
            ...
        }

        <userDict2>
        {
            ...
        }

        ...
    }
    \endverbatim

    with optional entries:
    \table
        Property | Description                            | Type | Reqd | Deflt
        libs     | Preloaded library names                | words | no  | -
        errors   | Error handling (default/warn/ignore/strict) | word | no | inherits
        useNamePrefix | Default enable/disable scoping prefix | bool | no | no-op
    \endtable

    The optional \c errors entry controls how FatalError is caught
    during construction and execute/write. FatalIOError is unaffected.
    - \c default : warn on construction errors, fatal on runtime errors
    - \c warn : warn on construction and runtime errors
    - \c ignore : ignore construction and runtime errors
    - \c strict : fatal on construction and runtime errors
    .

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjectList_H
#define Foam_functionObjectList_H

#include "PtrList.H"
#include "functionObject.H"
#include "SHA1Digest.H"
#include "HashTable.H"
#include "IOdictionary.H"
#include "HashSet.H"
#include "functionObjectProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class argList;
class mapPolyMesh;
class wordRe;

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    private PtrList<functionObject>
{
    // Private Data

        //- A list of error/warning handling
        List<error::handlerTypes> errorHandling_;

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests/errorHandling
        HashTable<label> indices_;

        //- Track the number of warnings per function object and limit
        //  to a predefined number to avoid flooding the display.
        //  Clear on re-read of functions.
        HashTable<unsigned> warnings_;

        //- Reference to Time
        const Time& time_;

        //- The parent dictionary containing a "functions" sub-dictionary
        //- of functionObject specifications
        const dictionary& parentDict_;

        //- Function object properties - stores state information
        mutable autoPtr<functionObjects::properties> propsDictPtr_;

        //- Function objects output registry
        mutable autoPtr<objectRegistry> objectsRegistryPtr_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- List of functions
        const PtrList<functionObject>& functions() const { return *this; }

        //- List of functions
        PtrList<functionObject>& functions() { return *this; }

        //- Create properties dictionary - attached to Time.
        void createPropertiesDict() const;

        //- Create registry for output objects - attached to Time.
        void createOutputRegistry() const;

        //- Remove and return the function object pointer by name,
        //- and returns the old index (into digest) via the parameter.
        //  Returns nullptr (and index -1) if it didn't exist
        autoPtr<functionObject> remove(const word& key, label& oldIndex);

        //- Search the specified directory for functionObject
        //- configuration files, add to the given map and recurse
        static void listDir(const fileName& dir, wordHashSet& available);

        //- Like Enum::getOrDefault, but with additional code to warn if
        //- the 'key' is not a primitive entry.
        //
        //  This additional treatment is to ensure that potentially existing
        //  code with an "errors" functionObject will continue to run.
        error::handlerTypes getOrDefaultErrorHandling
        (
            const word& key,
            const dictionary& dict,
            const error::handlerTypes deflt
        ) const;


        //- No copy construct
        functionObjectList(const functionObjectList&) = delete;

        //- No copy assignment
        void operator=(const functionObjectList&) = delete;


public:

    // Static data members

        //- Default relative path ("caseDicts/postProcessing") to the
        //- directory structure containing functionObject dictionary files.
        static fileName functionObjectDictPath;


    // Constructors

        //- Construct from Time and the execution setting.
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& runTime,
            const bool execution=true
        );

        //- Construct from Time, a dictionary with a "functions" entry
        //- and the execution setting.
        //  \param[in]  runTime - the other Time instance to construct from
        //  \param[in]  parentDict - the parent dictionary containing
        //    a "functions" sub-dictionary of functionObject specifications.
        //  \param[in]  execution - whether the function objects should execute
        //    or not. Default: true.
        functionObjectList
        (
            const Time& runTime,
            const dictionary& parentDict,
            const bool execution=true
        );

        //- Construct and return a functionObjectList for an application.
        //  If the "dict" argument is specified the functionObjectList is
        //  constructed from that dictionary which is returned as
        //  controlDict otherwise the functionObjectList is constructed
        //  from the "functions" sub-dictionary of "system/controlDict"
        static autoPtr<functionObjectList> New
        (
            const argList& args,
            const Time& runTime,
            dictionary& controlDict,
            HashSet<wordRe>& requiredFields
        );


    //- Destructor
    ~functionObjectList() = default;


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Return the current trigger index (read from the propsDict)
        label triggerIndex() const;

        //- Reset/read properties dictionary for current time
        void resetPropertiesDict();

        //- Write access to the properties dictionary
        //- ("functionObjectProperties") registered on Time
        functionObjects::properties& propsDict();

        //- Const access to the properties dictionary
        //- ("functionObjectProperties") registered on Time
        const Foam::functionObjects::properties& propsDict() const;

        //- Write access to the output objects ("functionObjectObjects")
        //- registered on Time
        objectRegistry& storedObjects();

        //- Const access to the output objects ("functionObjectObjects")
        //- registered on Time
        const objectRegistry& storedObjects() const;

        //- Clear the list of function objects
        void clear();

        //- Find the ID of a given function object by name, -1 if not found.
        label findObjectID(const word& objName) const;

        //- Print a list of functionObject configuration files in the
        //- directories located using
        //- Foam::findEtcDirs("caseDicts/postProcessing")
        //
        //  -# \b user settings
        //    - ~/.OpenFOAM/{PROJECT_API}/"caseDicts/postProcessing"
        //    - ~/.OpenFOAM/"caseDicts/postProcessing"
        //  -# \b group settings
        //    - $WM_PROJECT_SITE/{PROJECT_API}/"etc/caseDicts/postProcessing"
        //    - $WM_PROJECT_SITE/"etc/caseDicts/postProcessing"
        //  -# \b other (shipped) settings
        //    - $WM_PROJECT_DIR/etc/"caseDicts/postProcessing"
        //
        // Where {PROJECT_API} is the value of the OPENFOAM define.
        // See further notes in Foam::findEtcEntries()
        static void list();

        //- Find a functionObject dictionary file in the case
        //- \<system\> directory or any directory located using
        //- Foam::findEtcDirs("caseDicts/postProcessing")
        //
        //  \return The path of the functionObject dictionary file found
        //     or an empty path
        static fileName findDict(const word& funcName);

        //- Read the specified functionObject configuration dictionary parsing
        //- the optional arguments included in the name 'funcNameArgs0',
        //- inserting 'field' or 'fields' entries as required and merging the
        //- resulting functionObject dictionary into 'functionsDict'.  Any
        //- fields required to execute the functionObject are added to
        //- 'requiredFields'
        //
        //  Uses functionObjectList::findDict() for searching
        static bool readFunctionObject
        (
            const string& funcNameArgs0,
            dictionary& functionsDict,
            HashSet<wordRe>& requiredFields,
            const word& region = word::null
        );

        //- Read and set the function objects if their data have changed
        bool read();

        //- Switch the function objects on
        void on();

        //- Switch the function objects off
        void off();

        //- Return the execution status (on/off) of the function objects
        bool status() const;

        //- Called at the start of the time-loop
        bool start();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        bool execute(bool writeProperties = true);

        //- Execute function objects using the specified subIndex.
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const label subIndex);

        //- Execute a subset of function objects using the specified subIndex.
        //  \param functionNames names or regex of existing functions to
        //      execute
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar
        bool execute(const UList<wordRe>& functionNames, const label subIndex);

        //- Called when Time::run() determines that the time-loop exits
        bool end();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        bool adjustTimeStep();

        //- Did any file get changed during execution?
        bool filesModified() const;

        //- Update for changes of mesh
        void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
