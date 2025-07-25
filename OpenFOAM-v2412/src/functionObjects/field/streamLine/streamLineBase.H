/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::streamLineBase

SeeAlso
    Foam::functionObjects::streamLine
    Foam::functionObjects::wallBoundedStreamLine

SourceFiles
    streamLineBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_streamLineBase_H
#define Foam_functionObjects_streamLineBase_H

#include "fvMeshFunctionObject.H"
#include "DynamicList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "coordSetWriter.H"
#include "indirectPrimitivePatch.H"
#include "interpolation.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class meshSearch;
class sampledSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class streamLineBase Declaration
\*---------------------------------------------------------------------------*/

class streamLineBase
:
    public functionObjects::fvMeshFunctionObject
{
public:

    // Public data types

        //- Enumeration defining the track direction
        enum trackDirType : char
        {
            FORWARD,            //!< Use "forward" tracking
            BACKWARD,           //!< Use "backward" tracking
            BIDIRECTIONAL,      //!< Use "bidirectional" tracking
        };

        //- Names for the trackDir
        static const Enum<trackDirType> trackDirTypeNames;


protected:

        //- Seed set engine
        mutable autoPtr<sampledSet> sampledSetPtr_;

        //- Axis of the sampled points to output
        mutable word sampledSetAxis_;

        //- Input dictionary
        dictionary dict_;

        //- List of fields to sample
        wordList fields_;

        //- Field to transport particle with
        word UName_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- Whether to use +U, -U or both
        trackDirType trackDir_;

        //- Maximum lifetime (= number of cells) of particle
        label lifeTime_;

        //- Track length
        scalar trackLength_;

        //- Optional trimming of tracks
        boundBox bounds_;

        //- Optional specified name of particles
        word cloudName_;

        //- Type of seed
        word seedSet_;

        //- Names of scalar fields
        wordList scalarNames_;

        //- Names of vector fields
        wordList vectorNames_;


    // Demand Driven

        //- File writer for tracks data
        mutable autoPtr<coordSetWriter> trackWriterPtr_;


    // Generated Data

        //- All tracks. Per track the points it passed through
        DynamicList<List<point>> allTracks_;

        //- Per scalarField, per track, the sampled values
        List<DynamicList<scalarList>> allScalars_;

        //- Per vectorField, per track, the sampled values
        List<DynamicList<vectorList>> allVectors_;


    // Protected Member Functions

        //- The axis of the sampledSet. Creates sampledSet if required.
        const word& sampledSetAxis() const;

        //- Demand driven construction of the sampledSet.
        //  Also updates sampledSetAxis_
        const sampledSet& sampledSetPoints() const;

        //- Construct patch out of all wall patch faces
        autoPtr<indirectPrimitivePatch> wallPatch() const;

        //- Initialise interpolators and track storage
        //  Return velocity interpolator: standalone or part of vector
        //  interpolators
        refPtr<interpolation<vector>> initInterpolations
        (
            const label nSeeds,
            PtrList<interpolation<scalar>>& vsInterp,
            PtrList<interpolation<vector>>& vvInterp
        );

        //- Generate point and values by interpolating from existing values
        void storePoint
        (
            const label tracki,

            const scalar w,
            const label lefti,
            const label righti,

            DynamicList<point>& newTrack,
            DynamicList<List<scalar>>& newScalars,
            DynamicList<List<vector>>& newVectors
        ) const;

        //- Trim and possibly split a track
        void trimToBox
        (
            const treeBoundBox& bb,
            const label tracki,
            PtrList<DynamicList<point>>& newTracks,
            PtrList<DynamicList<scalarList>>& newScalars,
            PtrList<DynamicList<vectorList>>& newVectors
        ) const;

        //- Trim tracks to bounding box
        void trimToBox(const treeBoundBox& bb);

        //- Do the actual tracking to fill the track data
        virtual void track() = 0;

        //- Write tracks to file
        virtual bool writeToFile();

        //- Reset the field names
        virtual void resetFieldNames
        (
            const word& newUName,
            const wordList& newFieldNames
        );


public:

    //- Runtime type information
    TypeName("streamLineBase");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        streamLineBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from Time and dictionary and list of fields to sample
        streamLineBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const wordList& fieldNames
        );


    //- Destructor
    virtual ~streamLineBase();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Execute the averaging
        virtual bool execute();

        //- Track and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
