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
    Foam::adiosFoam::adiosRegionControl

Description
    Controls for managing region and cloud writing

SourceFiles
    adiosRegionControl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adiosFoamRegionControl_H
#define Foam_adiosFoamRegionControl_H

#include "adiosTime.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cloud;
class fvMesh;

namespace adiosFoam
{

/*---------------------------------------------------------------------------*\
                        Class regionControl Declaration
\*---------------------------------------------------------------------------*/

//- Per-region info variables grouped together to be able
//- to create a list of them
class regionControl
{
    // Private Data

        //- Name of the region
        word name_;

        //- Normally don't want an explicit-write (want auto-write)
        bool explicitWrite_;

        //- Names of volume fields to write
        wordRes requestedFields_;

        //- Names of clouds to write
        wordRes requestedClouds_;

        //- Names of volume fields to ignore
        wordRes ignoredFields_;

        //- Names of clouds to ignored
        wordRes ignoredClouds_;


    // Track changes to mesh status (points/faces)

        //- Time of the topology status at the last mesh change
        mutable scalar topoTime_;

        //- Time of the point status at the last mesh change
        mutable scalar pointTime_;

        //- Name of the most recent faces instance
        mutable word facesName_;

        //- Name of the most recent points instance
        mutable word pointsName_;


    // Private Member Functions

        //- Read info for one region from its own dictionary
        void read(const fvMesh& mesh, const dictionary& dict);

        //- Reset all information
        void resetAll();


public:

    // Constructors

        //- Default construct
        regionControl();

        //- Default construct with specified region name
        explicit regionControl(const word& regionName);

        //- Construct from dictionary information
        regionControl(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    ~regionControl() = default;


    // Private Member Functions

        //- Name of the region
        const word& name() const noexcept { return name_; }

        //- Time of the topology status at the last mesh change
        scalar topoTime() const noexcept { return topoTime_; }

        //- Time of the point status at the last mesh change
        scalar pointTime() const noexcept { return pointTime_; }

        //- If auto-write is in effect
        bool autoWrite() const noexcept { return !explicitWrite_; }


        //- Accept cloud based on its name
        bool acceptCloudObject(const cloud& c) const;

        //- Accept field io-object based on its name
        bool acceptFieldObject(const regIOobject& io) const;


        //- Update topology changing information (as required)
        //  Return the state (unchanged, points-moved, topo-change)
        enum polyMesh::readUpdateState updateTimes(const fvMesh&) const;


        //- Path name for the region
        fileName regionPath() const { return adiosFoam::regionPath(name_); }

        //- Simple text representation for output purposes
        word info() const
        {
            return name_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adiosFoam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
