/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::Time

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    TimeNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Time_H
#define Foam_Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "unwatchedIOdictionary.H"
#include "FIFOStack.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "Enum.H"
#include "typeInfo.H"
#include "dlLibraryTable.H"
#include "functionObjectList.H"
#include "sigWriteNow.H"
#include "sigStopAtWriteNow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class argList;
class profilingTrigger;
class OSstream;

/*---------------------------------------------------------------------------*\
                             Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
public:

        //- Write control options
        enum writeControls
        {
            wcNone,                 //!< "none"
            wcTimeStep,             //!< "timeStep"
            wcRunTime,              //!< "runTime"
            wcAdjustableRunTime,    //!< "adjustable" / "adjustableRunTime"
            wcClockTime,            //!< "clockTime"
            wcCpuTime,              //!< "cpuTime"
            wcUnknown               //!< Dummy no-op
        };

        //- Stop-run control options, which are primarily used when
        //- altering the stopAt condition.
        enum stopAtControls
        {
            saEndTime,    //!< Stop when Time reaches prescribed endTime
            saNoWriteNow, //!< Adjust endTime to stop immediately w/o writing
            saWriteNow,   //!< adjust endTime to stop immediately w/ writing
            saNextWrite,  //!< stop at the next data write interval
            saUnknown     //!< Dummy no-op. Do not change current value.
        };

        //- Names for writeControls
        static const Enum<writeControls> writeControlNames;

        //- Names for stopAtControls
        static const Enum<stopAtControls> stopAtControlNames;


        //- Style for "ExecutionTime = " output
        //  0 = seconds (with trailing 's')
        //  1 = day-hh:mm:ss
        //
        // \note this is public so registered info switches can modify it.
        static int printExecutionFormat_;


private:

    // Private Data

        //- Profiling trigger for time-loop (for run, loop)
        mutable std::unique_ptr<profilingTrigger> loopProfiling_;

        //- Any loaded dynamic libraries
        //  Construct before reading controlDict
        mutable dlLibraryTable libs_;

        //- The controlDict
        unwatchedIOdictionary controlDict_;


protected:

    // Protected Data

        label startTimeIndex_;

        scalar startTime_;

        mutable scalar endTime_;

        mutable stopAtControls stopAt_;

        writeControls writeControl_;

        scalar writeInterval_;

        label purgeWrite_;

        mutable FIFOStack<word> previousWriteTimes_;

        //- The total number of sub-cycles, the current sub-cycle index,
        //- or 0 if time is not being sub-cycled
        label subCycling_;

        // One-shot writing
        bool writeOnce_;

        //- If time is being sub-cycled this is the previous TimeState
        autoPtr<TimeState> prevTimeState_;

        //- Signal handler for one-shot writing upon signal
        sigWriteNow sigWriteNow_;

        //- Signal handler for write and clean exit upon signal
        sigStopAtWriteNow sigStopAtWriteNow_;

        //- Format for time directory names (general | fixed | scientific)
        static IOstreamOption::floatFormat format_;

        //- Time directory name precision
        static int precision_;

        //- Maximum time directory name precision
        static const int maxPrecision_;


        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Set the controls from the current controlDict
        void setControls();

        //- Set file monitoring, profiling, etc
        //  Optionally force profiling without inspecting the controlDict
        void setMonitoring(const bool forceProfiling=false);

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- The write stream option (format, compression, version)
        IOstreamOption writeStreamOption_;

        //- Default graph format
        word graphFormat_;

        //- Is runtime modification of dictionaries allowed?
        Switch runTimeModifiable_;

        //- Is temporary object cache enabled?
        mutable bool cacheTemporaryObjects_;

        //- Function objects executed at start and on ++, +=
        mutable functionObjectList functionObjects_;


public:

    // Static Data Members

        //- The default control dictionary name (normally "controlDict")
        static word controlDictName;


    //- Runtime type information
    TypeName("time");


    // Constructors

    // Construct from argList

        //- Construct from argument list,
        //- reading from specified control dictionary name
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const word& ctrlDictName,
            const argList& args,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );

        //- Construct from argument list,
        //- reading from specified control dictionary name.
        //  Use specified system/constant directory names.
        Time
        (
            const word& ctrlDictName,
            const argList& args,
            const word& systemDirName,      //!< "system" directory name
            const word& constantDirName,    //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );


    // Construct from root and case paths

        //- Construct given name of dictionary to read and root/case paths.
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const word& ctrlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );

        //- Construct given name of dictionary to read and root/case paths.
        //  Use specified system/constant directory names.
        Time
        (
            const word& ctrlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemDirName,      //!< "system" directory name
            const word& constantDirName,    //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );

        //- Construct from root/case paths,
        //- with initial dictionary content
        //- \em without reading from "controlDict".
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );

        //- Construct from root/case paths,
        //- with initial dictionary content
        //- \em without reading from "controlDict".
        //  Use specified system/constant directory names.
        Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemDirName,      //!< "system" directory name
            const word& constantDirName,    //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true,
            //! Read option for controlDict and objectRegistry
            IOobjectOption::readOption rOpt = IOobjectOption::READ_MODIFIED
        );

        //- Construct from root/case paths
        //- \em without reading or requiring a controlDict or other files.
        //  Uses "system" and "constant" for its directories.
        //  Uses default controls values
        //  (write: timeStep, stop: endTime, format: ASCII, ...),
        inline Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct from root/case paths
        //- \em without reading or requiring a controlDict or other files.
        //  Use specified system/constant directory names.
        //  Uses default controls values
        //  (write: timeStep, stop: endTime, format: ASCII, ...)
        Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemDirName,      //!< "system" directory name
            const word& constantDirName,    //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );


    // Factory Methods

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New();

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New(const fileName& caseDir);

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New(const argList& args);

        //- Construct (dummy) global Time - no functionObjects or libraries,
        //- using the global path information stored in the FOAM_CASE
        //- environment. See argList::envGlobalPath()
        static autoPtr<Time> NewGlobalTime();

        //- Construct (dummy) global Time - no functionObjects or libraries,
        //- using the global path information from the referenced Time.
        static autoPtr<Time> NewGlobalTime(const Time& runTime);


    //- Destructor
    virtual ~Time();


    // Member Functions

    // TimeState Functions

        //- Return the current time name
        using TimeState::timeName;

        //- Return a time name for the given scalar time value
        //- formatted with the given precision
        static word timeName(const scalar t, const int precision = precision_);


    // Database names and paths

        //- Use name from objectRegistry, not TimePaths
        using objectRegistry::name;

        //- Return the rootPath
        using TimePaths::rootPath;

        //- Return global case name
        using TimePaths::globalCaseName;

        //- Return case name
        using TimePaths::caseName;

        //- Return path = rootPath/caseName. Same as TimePaths::path()
        fileName path() const
        {
            return TimePaths::rootPath()/TimePaths::caseName();
        }

        //- Return global path for the case = rootPath/globalCaseName.
        //- Same as TimePaths::globalPath()
        fileName globalPath() const
        {
            return TimePaths::rootPath()/TimePaths::globalCaseName();
        }

        //- Return current time path = path/timeName
        fileName timePath() const
        {
            return TimePaths::path()/TimeState::timeName();
        }

        //- Return current time global path = globalPath/timeName
        fileName globalTimePath() const
        {
            return TimePaths::globalPath()/TimeState::timeName();
        }


    // Database Functions

        //- Return read access to the controlDict dictionary
        const dictionary& controlDict() const noexcept { return controlDict_; }

        virtual const fileName& dbDir() const { return fileName::null; }

        //- Get write stream option (format, compression, version)
        inline IOstreamOption writeStreamOption() const noexcept;

        //- Get write stream format
        inline IOstreamOption::streamFormat writeFormat() const noexcept;

        //- Set the write stream format and return the previous value.
        //  This change will only effective until the next readModified.
        //  As a side-effect (eg, changing from ASCII to BINARY)
        //  it may also disable compression
        inline IOstreamOption::streamFormat
        writeFormat(IOstreamOption::streamFormat fmt) noexcept;

        //- Get the write stream compression
        inline IOstreamOption::compressionType writeCompression()
        const noexcept;

        //- Get the write stream version
        inline IOstreamOption::versionNumber writeVersion() const noexcept;

        //- Default graph format
        const word& graphFormat() const noexcept { return graphFormat_; }


    // Reading

        //- Supports re-reading
        Switch runTimeModifiable() const noexcept
        {
            return runTimeModifiable_;
        }

        //- Set re-reading support on/off (use with caution).
        //  \return the previous value
        Switch runTimeModifiable(Switch sw) noexcept
        {
            Switch old(runTimeModifiable_);
            runTimeModifiable_ = sw;
            return old;
        }

        //- Read control dictionary, update controls and time
        virtual bool read();

        //- Read the objects that have been modified
        void readModifiedObjects();


    // Searching

        //- Return time instance (location) of \c directory containing
        //- the file \c name (eg, used in reading mesh data).
        //- When \c name is empty, searches for \c directory only.
        //- Does not search beyond \c stopInstance (if set) or \c constant.
        //
        //  If the instance cannot be found:
        //  - FatalError when readOpt is MUST_READ or READ_MODIFIED
        //  - return \c stopInstance (if set and reached)
        //  - return \c constant if constant_fallback is true
        //  - return an empty word if constant_fallback is false
        //  .
        word findInstance
        (
            //! The subdirectory (local) for the search
            const fileName& directory,
            //! The filename for the search. If empty, only search for directory
            const word& name = word::null,
            //! The search type : generally MUST_READ or READ_IF_PRESENT
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            //! The search stop instance
            const word& stopInstance = word::null,
            //! Return \c "constant" instead of \c "" if the search failed
            const bool constant_fallback = true
        ) const;


    // Member Functions

        // Writing

            //- Write time dictionary to the \<time\>/uniform directory
            virtual bool writeTimeDict() const;

            //- Write using stream options
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool writeOnProc
            ) const;

            //- Write the objects immediately (not at end of iteration)
            //- and continue the run
            bool writeNow();

            //- Write the objects now (not at end of iteration) and end the run
            bool writeAndEnd();

            //- Write the objects once (one shot) and continue the run
            void writeOnce();

            //- Print the elapsed ExecutionTime (cpu-time), ClockTime
            Ostream& printExecutionTime(OSstream& os) const;


        // Access

            //- Return start time index
            virtual label startTimeIndex() const;

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;

            //- Return the stop control information
            virtual stopAtControls stopAt() const;

            //- Return true if adjustTimeStep is true
            virtual bool isAdjustTimeStep() const;

            //- Return the list of function objects
            const functionObjectList& functionObjects() const noexcept
            {
                return functionObjects_;
            }

            //- Mutable access to the loaded dynamic libraries
            dlLibraryTable& libs() const noexcept
            {
                return libs_;
            }

            //- Zero (tests as false) if time is not being sub-cycled,
            //- otherwise the current sub-cycle index or the total number of
            //- sub-cycles.
            //  The interpretation of non-zero values is dependent on the
            //  routine.
            label subCycling() const noexcept
            {
                return subCycling_;
            }

            //- Return previous TimeState if time is being sub-cycled
            const TimeState& prevTimeState() const
            {
                return *prevTimeState_;
            }


        // Check

            //- Return true if run should continue,
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.run())
            //      {
            //          ++runTime;
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool run() const;

            //- Return true if run should continue and if so increment time
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.loop())
            //      {
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool loop();

            //- Return true if end of run,
            //  does not invoke any functionObject methods
            //  \note
            //      The rounding heuristics near endTime mean that
            //      \code run() \endcode and \code !end() \endcode may
            //      not yield the same result
            virtual bool end() const;


        // Edit

            //- Adjust the current stopAtControl.
            //  \param stopCtrl the new stop control, whereby
            //      stopAtControls::saUnknown is treated as a no-op.
            //  \note this value only persists until the next time the
            //      dictionary is read.
            //  \return true if the stopAt() value was changed.
            virtual bool stopAt(const stopAtControls stopCtrl) const;

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time& t);

            //- Reset the time and time-index
            virtual void setTime(const instant& inst, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime
            (
                const dimensionedScalar& newTime,
                const label newIndex
            );

            //- Reset the time and time-index
            virtual void setTime(const scalar newTime, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar& endTime);

            //- Reset end time
            virtual void setEndTime(const scalar endTime);

            //- Reset time step, normally also calling adjustDeltaT()
            virtual void setDeltaT
            (
                const dimensionedScalar& deltaT,
                const bool adjust = true
            );

            //- Reset time step, normally also calling adjustDeltaT()
            virtual void setDeltaT
            (
                const scalar deltaT,
                const bool adjust = true
            );

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Adjust the reported sub-cycle index.
            //  \param index is the sub-cycle index.
            //      This index is ignored sub-cycling was
            //      not already registered, or if the index is zero or
            //      negative.
            virtual void subCycleIndex(const label index);

            //- Reset time after sub-cycling back to previous TimeState
            virtual void endSubCycle();

            //- Return non-const access to the list of function objects
            functionObjectList& functionObjects() noexcept
            {
                return functionObjects_;
            }


    // Member Operators

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const dimensionedScalar& deltaT);

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const scalar deltaT);

        //- Prefix increment,
        //  also invokes the functionObjectList::start() or
        //  functionObjectList::execute() method, depending on the time-index
        virtual Time& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
