/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

Namespace
    Foam::functionObjects

Description
    Function objects are OpenFOAM utilities to ease workflow configurations and
    enhance workflows by producing additional user-requested data both during
    runtime and postprocessing calculations, typically in the form of
    additional logging to the screen, or generating text, image and field files.

    Function objects eliminate the need to store all runtime generated data,
    hence saving considerable resources. Furthermore, function objects are
    readily applied to batch-driven processes, improving reliability by
    standardising the sequence of operations and reducing the amount of manual
    interaction.

    In addition, the output of most function objects, e.g. output fields, are
    stored on the mesh database so that it can be retrieved and used for other
    applications (e.g. directly using \c wallShearStress function object output
    in \c fieldAverage function object to produce \c wallShearStressMean field).

    \section secFunctionObjects Using function objects

    Function objects can be executed by using two methods:

      - \c functions sub-dictionary in the \c system/controlDict file
      - \c postProcess command-line utility

    For the first method, each selected function object should be listed inside
    \c functions sub-dictionary of the \c system/controlDict file as a nested
    sub-dictionary, typically as in the following example:

    \verbatim
    functions   // sub-dictionary name under the system/controlDict file
    {
        ..optional entries..

        <dictName1>
        {
            // Mandatory entries
            type    <functionObjectTypeName>;
            libs    (<libType>FunctionObjects);

            // Mandatory entries defined in <functionObjectType>
            ...

            // Optional entries defined in <functionObjectType>
            ...

            // Optional (inherited) entries
            region              region0;
            enabled             true;
            log                 true;
            timeStart           0;
            timeEnd             1000;
            executeControl      timeStep;
            executeInterval     1;
            writeControl        timeStep;
            writeInterval       1;
        }

        <dictName2>
        {
            ...
        }

        ...

        <dictNameN>
        {
            ...
        }
    }
    \endverbatim

    where the entries mean:
    \table
        Property | Description                            | Type | Reqd | Deflt
        type     | Type name of function object           | word |  yes  | -
        libs     | Library name(s) for implementation     | words | no  | -
        enabled  | Switch to turn function object on/off  | bool | no | true
        errors   | Error handling (default/warn/ignore/strict) | word | no | inherits
        log     | Switch to write log info to standard output | bool | no | true
        useNamePrefix | Add scoping prefix to names       | bool | no | inherits
        region   | Name of region for multi-region cases  | word | no | region0
        timeStart | Start time for function object execution | scalar | no | 0
        timeEnd   | End time for function object execution   | scalar | no | inf
        executeControl  | See time controls below     | word | no    | timeStep
        executeInterval | Steps/time between execute phases  | label | no  | 1
        writeControl    | See time controls below     | word | no    | timeStep
        writeInterval   | Steps/time between write phases    | label | no  | 1
    \endtable

    If the \c errors entry is missing, it uses the value (if any)
    specified within the top-level functionObjectList.

    If the \c useNamePrefix entry is missing, it uses the value (if any)
    specified within the top-level functionObjectList or otherwise
    uses the current value of functionObject::defaultUseNamePrefix

    Time controls:
    \table
        Option            | Description
        none              | Trigger is disabled
        timeStep          | Trigger every 'Interval' time-steps
        writeTime         | Trigger every 'Interval' output times
        runTime           | Trigger every 'Interval' run time period
        adjustableRunTime | Currently identical to "runTime"
        clockTime         | Trigger every 'Interval' clock time period
        cpuTime           | Trigger every 'Interval' CPU time period
        onEnd             | Trigger on end of simulation run
    \endtable

    The sub-dictionary name \c \<userDefinedSubDictName\> is chosen by the user,
    and is typically used as the name of the output directory for any data
    written by the function object.

    As the base mandatory entries, the \c type entry defines the type of
    function object properties that follow.  Function objects are packaged into
    separate libraries for flexibility and the \c libs entry is used to specify
    which library should be loaded.

    Each function object has two separate run phases:

      - The \c execute phase is meant to be used for updating calculations
        or for management tasks.
      - The \c write phase is meant for writing the calculated data to disk.

    For each phase the respective time controls are provided, as listed above.


    The second method of executing function objects is to use \c postProcess
    utility.

    When specified without additional arguments, the \c postProcess utility
    executes all function objects defined in the \c system/controlDict file
    for all time directories:

    \verbatim
        postProcess
    \endverbatim

    Most function objects can be invoked directly without the need to specify
    the input dictionary using the \c -func option, e.g. to execute the Courant
    number function object:

    \verbatim
        postProcess -func CourantNo
    \endverbatim

    In addition, the \c -postProcess option is available to all solvers,
    and operates similarly to the stand-alone \c postProcess utility.
    For example, having completed a \c simpleFoam calculation, the following
    will execute all function objects defined in the \c system/controlDict file
    for all time directories:

    \verbatim
        simpleFoam -postProcess
    \endverbatim

Class
    Foam::functionObject

Description
    Abstract base-class for Time/database function objects.

See also
    - Foam::functionObjectList
    - Foam::functionObjects::timeControl

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObject_H
#define Foam_functionObject_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private Data

        //- Function object name
        const word name_;

        //- Flag to indicate that names should be prefixed
        bool useNamePrefix_;


protected:

    // Protected Member Functions

        //- Return a scoped (prefixed) name
        //  Used to construct local field names, controlled by useNamePrefix_
        word scopedName(const word& name) const;


public:

    // Forward Declarations
    class unavailableFunctionObject;


    // Public Data

        //- Flag to write log into Info
        bool log;


    // Static Data Members

        //- Flag to execute debug content
        static int debug;

        //- Global post-processing mode switch
        static bool postProcess;

        //- Global default for useNamePrefix
        static bool defaultUseNamePrefix;

        //- Directory prefix
        static word outputPrefix;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionObject,
            dictionary,
            (const word& name, const Time& runTime, const dictionary& dict),
            (name, runTime, dict)
        );


    // Constructors

        //- Construct from components
        explicit functionObject
        (
            const word& name,
            const bool withNamePrefix = defaultUseNamePrefix
        );

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Select from dictionary, based on its "type" entry
        static autoPtr<functionObject> New
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~functionObject() = default;


    // Member Functions

        //- Runtime type information
        virtual const word& type() const = 0;

        //- Return the name of this functionObject
        const word& name() const noexcept;

        //- Return the flag for adding a scoping name prefix
        bool useNamePrefix() const noexcept;

        //- Modify the flag for adding a scoping name prefix
        //  \return previous value.
        bool useNamePrefix(bool on) noexcept;

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute() = 0;

        //- Execute using the specified subIndex.
        //  The base implementation is a no-op.
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar.
        virtual bool execute(const label subIndex);

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write() = 0;

        //- Called when Time::run() determines that the time-loop exits.
        //  The base implementation is a no-op.
        virtual bool end();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Did any file get changed during execution?
        virtual bool filesModified() const;

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void movePoints(const polyMesh& mesh);
};


/*---------------------------------------------------------------------------*\
          Class functionObject::unavailableFunctionObject Declaration
\*---------------------------------------------------------------------------*/

//- Abstract functionObject to report when a real version is unavailable.
class functionObject::unavailableFunctionObject
:
    public functionObject
{
protected:

    //- Construct with name
    explicit unavailableFunctionObject(const word& name);

    //- Report it is unavailable, emitting a FatalError for try/catch
    //- in the caller
    void carp(std::string message = "") const;


public:

    // Member Functions

        //- No nothing
        virtual bool execute();

        //- No nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
