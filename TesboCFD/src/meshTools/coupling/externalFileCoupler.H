/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::externalFileCoupler

Description
    Encapsulates the logic for coordinating between OpenFOAM and an
    external application.

    This class provides a simple interface for explicit coupling with an
    external application using plain text files located in the user-specified
    communications directory.
    These files are to be read/written on the master processor only.

    The explicit coupling follows a master/slave model in which OpenFOAM
    is the 'master' and the external application is the 'slave'.
    The readiness to exchange information in either direction is handled
    by a lock file (ie, OpenFOAM.lock).
    If the lock file is present, the slave (external application) should wait
    for the master (OpenFOAM) to complete.

    When the master is finished its tasks and has prepared data for
    the slave, the lock file is removed, instructing the external
    source to take control of the program execution.

    When the slave has completed its tasks, it will reinstate the lock file.

    Example of the communication specification:
    \verbatim
    communication
    {
        commsDir        "<case>/comms";
        waitInterval    1;
        timeOut         100;
        initByExternal  no;
        statusDone      done;
    }
    \endverbatim

    A typical coupling loop would look like this (on the master-side):
    \verbatim
        initialize - master takes control
        writeDataMaster() - write data for slave
        useSlave()
        waitForSlave()
        removeDataMaster() - cleanup old data from master [optional?]
        readDataMaster() - read data from slave
        useMaster()
    \endverbatim

    On the slave-side:
    \verbatim
        waitForMaster()
        readDataSlave()  - read data from master
        writeDataSlave() - write data for master
        useMaster()
    \endverbatim

    Note that since the waitForSlave() method not only waits for the lock file
    to be reinstated but also does a simple check of its contents, it can
    also serve to communicate some control from the slave to the master.

SourceFiles
    externalFileCoupler.C
    externalFileCouplerI.H

\*---------------------------------------------------------------------------*/

#ifndef externalFileCoupler_H
#define externalFileCoupler_H

#include "fileName.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class externalFileCoupler Declaration
\*---------------------------------------------------------------------------*/

class externalFileCoupler
{
public:

    // Public Data Types

        //- The run state (ie, who is currently in charge)
        enum runState
        {
            NONE,       //!< Not initialized
            MASTER,     //!< The master (OpenFOAM) is in charge
            SLAVE,      //!< The slave (external program) is in charge
            DONE        //!< Finished
        };

private:

    // Private Data

        //- The current run (and initialization) state
        mutable runState runState_;

        //- Local path to communications directory
        fileName commsDir_;

        //- Value for "status=..." on termination
        word statusDone_;

        //- Interval time between checking for return data [s]
        unsigned waitInterval_;

        //- Timeout [s] while waiting for the external application
        unsigned timeOut_;

        //- Flag to indicate values are initialized by external program
        bool slaveFirst_;

        //- Local logging/verbosity flag
        bool log;


    // Private Member Functions

        //- No copy construct
        externalFileCoupler(const externalFileCoupler&) = delete;

        //- No copy assignment
        void operator=(const externalFileCoupler&) = delete;


public:

    //- Runtime type information
    TypeName("externalFileCoupler");


    // Static data members

        //- Name of the lock file
        static word lockName;


    // Constructors

        //- Construct using standard defaults.
        //  Does not create communications directory.
        externalFileCoupler();

        //- Construct with specified communications directory.
        //  Creates the communications directory upon construction.
        externalFileCoupler(const fileName& commsDir);

        //- Construct from dictionary
        //  Creates the communications directory upon construction.
        externalFileCoupler(const dictionary& dict);


    //- Destructor
    virtual ~externalFileCoupler();


    // Member Functions

    // Initialization

        //- True if state has been initialized
        inline bool initialized() const;

        //- External application provides initial values
        inline bool slaveFirst() const;


    // File locations

        //- Return the file path to the base communications directory
        inline const fileName& commDirectory() const;

        //- Return the file path in the communications directory
        inline fileName resolveFile(const word& file) const;

        //- Return the file path to the lock file
        inline fileName lockFile() const;


    // Settings

        //- Read communication settings from dictionary
        bool readDict(const dictionary& dict);


    // Handshaking

        //- Create lock file to indicate that OpenFOAM is in charge
        //  \param wait - wait for master to complete.
        //
        //  \return Time::stopAtControls::saUnknown if not waiting, otherwise
        //      as per the waitForMaster() description.
        enum Time::stopAtControls useMaster(const bool wait=false) const;

        //- Remove lock file to indicate that the external program is in charge
        //  \param wait - wait for slave to complete.
        //
        //  \return Time::stopAtControls::saUnknown if not waiting, otherwise
        //      as per the waitForSlave() description.
        enum Time::stopAtControls useSlave(const bool wait=false) const;


        //- Wait for master to complete.
        //  This is when the lock file disappears, or exists but has
        //  \c status=done content.
        //
        //  \return Time::stopAtControls::saUnknown or if lock file
        //      contained \c status=done it returns
        //      Time::stopAtControls::saEndTime
        enum Time::stopAtControls waitForMaster() const;

        //- Wait for slave to complete.
        //  This is when the lock file appears.
        //
        //  When the lock file appears, it is checked for the following
        //  content which corresponds to particular return values:
        //  - \c status=done
        //  \return Foam::Time::saEndTime
        //  - \c action=noWriteNow
        //  \return Foam::Time::saNoWriteNow
        //  - \c action=writeNow
        //  \return Foam::Time::saWriteNow
        //  - \c action=nextWrite
        //  \return Foam::Time::saNextWrite
        //  - Anything else (empty file, no action= or status=, etc)
        //  \return Foam::Time::saUnknown
        enum Time::stopAtControls waitForSlave() const;


    // File creation, removal

        //- Read data files on master (OpenFOAM).
        //  These data files are normally created by the slave.
        virtual void readDataMaster();

        //- Read data files on slave (external program).
        //  These data files are normally created by the master.
        virtual void readDataSlave();


        //- Write data files from master (OpenFOAM)
        virtual void writeDataMaster() const;

        //- Write data files from slave (external program)
        virtual void writeDataSlave() const;


        //- Remove data files written by master (OpenFOAM)
        virtual void removeDataMaster() const;

        //- Remove data files written by slave (external program)
        virtual void removeDataSlave() const;


        //- Generate \c status=done in lock (only when run-state = master)
        //  The exact text can be specified via the statusDone keyword
        void shutdown() const;

        //- Remove files written by OpenFOAM
        void removeDirectory() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "externalFileCouplerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
