/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::adiosWrite

Description
    Writes fields and parcel clouds to an ADIOS file.

    restartFrom     (none | restartTime | latestTime);
    restartTime     float;

    stopAt          (none | now | stopTime);
    stopTime        float;

Note
    Original code concept (adios1) 2015 Norbert Podhorszki

SourceFiles
    adiosRestart.C
    adiosWrite.C
    adiosWriteField.C
    adiosWriteFunctionObject.C
    adiosRegion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adiosWrite_H
#define Foam_adiosWrite_H

#include "adiosCoreWrite.H"
#include "adiosReader.H"
#include "adiosRegionControl.H"
#include "adiosTime.H"

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"
#include "timeFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class adiosWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosWrite
:
    public functionObjects::timeFunctionObject,
    public adiosFoam::adiosCoreWrite
{
    // Private Data

        //- Requested or current restart state
        enum restartState
        {
            restartNone,        //!< No restart, or not initialized
            restartTime,        //!< Restart based on specified time
            restartLatest,      //!< Restart from latest file available
            restartDisabled     //!< Restart disabled or already executed
        };

        //- Requested or current restart state
        enum stopState
        {
            stopNone,           //!< No stopping control
            stopNow,            //!< Manual stop immediately (after ADIOS write)
            stopTime,           //!< Stop if the time at or exceeds stopTime
            stopDisabled        //!< Stopping disabled or already triggered
        };


        //- The ADIOS data directory
        fileName dataDir_;

        //- The restart state (requested or current)
        restartState restartType_;

        //- The stop state (requested or current)
        stopState stopAt_;

        //- The restart time requested
        scalar restartTime_;

        //- Stop time requested
        scalar stopTime_;

        //- The time index for restart
        //  Also avoids write immediately after restart.
        int restartIndex_;

        //- List of regions holding all data needed (mesh, names)
        SLList<adiosFoam::regionControl> regions_;


protected:

    // Protected Member Functions

    // General Functions

        //- Create and open dataset for all data at given time-step
        bool open(const Time& runTime);

        //- Restart
        bool restart();

        //- Write ADIOS data for current time-step
        //  Meshes are written if regions have been updated
        void writeData();


    // Functions for writing (fields, clouds)

        //- Define and write fields for the region
        //  Return the number of fields defined.
        label writeFields
        (
            const adiosFoam::regionControl& regCtrl,
            bool verbose = false
        );

        //- Define/write ADIOS cloud variables/attributes
        //  for the region.
        //  Return the number of clouds written.
        label writeClouds(const adiosFoam::regionControl& regCtrl);

        //- Print cloud (debug only)
        template<class CloudType>
        label printCloud(const CloudType& c);



    // Functions for restart, reading (fields, clouds)

        //- Read data from current time
        //  Attempt to read data for all variables existing in memory
        //  Returns time information on success, invalid on failure
        adiosFoam::adiosTime readData();

        adiosFoam::adiosTime readData(const fileName& bpFile);
        adiosFoam::adiosTime readData(const instant& when);


public:

    //- No copy construct
    adiosWrite(const adiosWrite&) = delete;

    //- No copy assignment
    void operator=(const adiosWrite&) = delete;


    //- Runtime type information
    TypeName("adiosWrite");


    // Constructors

        //- Construct from Time and dictionary
        adiosWrite
        (
            const word& groupName,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~adiosWrite();


    // Member Functions

        //- Read the adiosWrite specification (including restart request)
        virtual bool read(const dictionary& dict);

        //- Execute, handles restart at first call or write initial conditions.
        virtual bool execute();

        //- Write a timestep to file
        virtual bool write();


    // From polyMesh

    //- Set io to be up-to-date with points
    // virtual void setUpToDatePoints(regIOobject& io) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
