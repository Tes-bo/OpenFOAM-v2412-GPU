/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2015,2023 OpenCFD Ltd.
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
    Foam::surfaceZonesInfo

Description

SourceFiles
    surfaceZonesInfo.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceZonesInfo_H
#define surfaceZonesInfo_H

#include "Enum.H"
#include "pointField.H"
#include "word.H"
#include "PtrList.H"
#include "labelList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;
class searchableSurfaces;
class polyMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                      Class surfaceZonesInfo Declaration
\*---------------------------------------------------------------------------*/

class surfaceZonesInfo
{
public:

    //- Types of selection of area
    enum areaSelectionAlgo
    {
        INSIDE,
        OUTSIDE,
        INSIDEPOINT,
        NONE
    };

    static const Enum<areaSelectionAlgo> areaSelectionAlgoNames;

    //- How to generate faceZone name
    enum faceZoneNaming
    {
        NOZONE,
        SINGLE,
        REGION
    };

    static const Enum<faceZoneNaming> faceZoneNamingNames;

    //- What to do with faceZone faces
    enum faceZoneType
    {
        INTERNAL,
        BAFFLE,
        BOUNDARY
    };

    static const Enum<faceZoneType> faceZoneTypeNames;


private:

    // Private data

        //- Per 'interface' surface : names of faceZones to put faces into
        wordList faceZoneNames_;

        //- Per 'interface' surface : name of cellZone to put cells into
        word cellZoneName_;

        //- Per 'interface' surface : (only used if surface is closed)
        //  How to select zone cells : surface inside or outside or given
        //  inside location.
        areaSelectionAlgo zoneInside_;

        //- If zoneInside=location gives the corresponding inside point
        pointField zoneInsidePoints_;

        //- Per 'interface' surface :
        //  What to do with outside
        faceZoneType faceType_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const surfaceZonesInfo&) = delete;


public:

    // Constructors

        //- Construct from surfaces and dictionary and fully resolved
        //  region names (for optional automatic faceZone naming)
        surfaceZonesInfo
        (
            const searchableSurface& surface,
            const dictionary& surfacesDict,
            const wordList& regionNames
        );

        //- Construct from components
        surfaceZonesInfo
        (
            const wordList& faceZoneNames,
            const word& cellZoneNames,
            const areaSelectionAlgo& zoneInside,
            const pointField& zoneInsidePoints,
            const faceZoneType& faceType
        );

        //- Copy constructor
        surfaceZonesInfo(const surfaceZonesInfo&);

        //- Return clone
        autoPtr<surfaceZonesInfo> clone() const
        {
            return autoPtr<surfaceZonesInfo>::New(*this);
        }


    // Member Functions

        // Access

            //- Per 'interface' surface : empty or names of faceZones to put
            //  faces into (according to region)
            const wordList& faceZoneNames() const
            {
                return faceZoneNames_;
            }

            //- Per 'interface' surface : empty or name of cellZone to put
            //  cells into
            const word& cellZoneName() const
            {
                return cellZoneName_;
            }

            const areaSelectionAlgo& zoneInside() const
            {
                return zoneInside_;
            }

            //- Get specified inside location for surfaces with a cellZone
            const point zoneInsidePoint() const
            {
                return zoneInsidePoints_[0];
            }

            //- Get specified inside locations for surfaces with a cellZone
            const pointField& zoneInsidePoints() const
            {
                return zoneInsidePoints_;
            }

            //- How to handle face of surfaces with a faceZone
            const faceZoneType& faceType() const
            {
                return faceType_;
            }


        // Query

            //- Get indices of unnamed surfaces (surfaces without faceZoneName)
            static labelList getUnnamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named surfaces (surfaces with faceZoneName)
            static labelList getNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named surfaces without a cellZone
            static labelList getStandaloneNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of surfaces with a cellZone that are closed and
            //  have 'inside' or 'outside' selection.
            static labelList getClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are unclosed
            static labelList getUnclosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are closed.
            static labelList getAllClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that have 'insidePoint'
            //  section.
            static labelList getInsidePointNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            static label addCellZone
            (
                const word& name,
                const labelList& addressing,
                polyMesh& mesh
            );

            static labelList addCellZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );

            static label addFaceZone
            (
                const word& name,
                const labelList& addressing,
                const boolList& flipMap,
                polyMesh& mesh
            );

            static labelListList addFaceZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
