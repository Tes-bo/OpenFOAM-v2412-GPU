/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::holeToFace

Description
    A topoSetFaceSource to select a set of faces that closes a hole i.e.
    disconnects zones (specified by locations) from one another.

    Algorithm roughly according to
    "A 3D-Hole Closing Algorithm", Zouina Aktouf et al

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        points      | Per zone the list of points       | yes   |
        faceSet     | Optional blocked faces            | no    | <empty>
        cellSet     | Optional subset of cells to operate in | no    | <empty>
        erode       | Perform some cleanup on set       | no    | no
    \endtable

    Limited to max 31 zones.

SourceFiles
    holeToFace.C

\*---------------------------------------------------------------------------*/

#ifndef holeToFace_H
#define holeToFace_H

#include "topoSetFaceSource.H"
#include "bitSet.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class holeToFace Declaration
\*---------------------------------------------------------------------------*/

class holeToFace
:
    public topoSetFaceSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Per 'zone' the set of points
        List<pointField> zonePoints_;

        //- Names of faceSets specifying blockage
        wordList blockedFaceNames_;

        //- Names of cellSets specifying blockage
        wordList blockedCellNames_;

        //- Erode set
        bool erode_;


    // Private Member Functions

        // Helper functions

            //- Expand pointField into list of pointFields
            static List<pointField> expand(const pointField& pts);

            //- Write .obj file with selected faces
            void writeFaces
            (
                const word& name,
                const bitSet& faceFld
            ) const;


        //- Calculate topological distance from seedFaces. Do not cross
        //  blocked cells/faces
        void calculateDistance
        (
            const labelList& seedFaces,
            const bitSet& isBlockedCell,
            const bitSet& isBlockedFace,
            labelList& cellDist,
            labelList& faceDist
        ) const;

        //- Calculate faces on outside of hole cells. Ignore blocked faces
        //  (isSurfaceFace) and faces reachable from more than one location
        bitSet frontFaces
        (
            const bitSet& isSurfaceFace,
            const List<unsigned int>& locationFaces,
            const bitSet& isHoleCell
        ) const;

        //- Overall driver to find minimum set of faces to close the path
        //  between zones
        bitSet findClosure
        (
            const bitSet& isSurfaceFace,                // intersected faces
            const bitSet& isCandidateHoleCell,          // starting blockage
            const labelListList& locationCells          // cells per zone
        ) const;

        //- WIP. Remove excess faces. Not parallel consistent.
        bitSet erodeSet
        (
            const bitSet& isSurfaceFace,
            const bitSet& isSetFace
        ) const;


public:

    //- Runtime type information
    TypeName("holeToFace");

    // Constructors

        //- Construct from components
        holeToFace
        (
            const polyMesh& mesh,
            const List<pointField>& zonePoints,
            const wordList& blockedFaceNames,
            const wordList& blockedCellNames,
            const bool erode
        );

        //- Construct from dictionary
        holeToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        holeToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~holeToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

        //- Optional direct use to generate a faceSet
        void combine
        (
            topoSet& set,
            const bitSet& isBlockedFace,
            const bitSet& isActiveCell,
            const bool add
        ) const;

        //- Optional direct use to generate the set of faces and the method to
        //  get data from nearby blocked faces. Gets provided with the
        //  - set of points per 'zone'
        //  - set of blocked faces
        //  - global numbering for these blocked faces
        //  Returns
        //  - the set of faces to disconnect zones from one
        //    another (closureFaces)
        //  - a nullptr or the maps from closureFaces to nearest blocked face
        //    (mapDistribute and closureToBlocked). Note: closureToBlocked
        //     (index into the map) can contain -1 if no nearby valid blocked
        //    face could be found (from an edge-face-edge walk)
        static autoPtr<mapDistribute> calcClosure
        (
            const polyMesh& mesh,
            const List<pointField>& zonePoints,
            const labelList& blockedFaces,
            const globalIndex& globalBlockedFaces,
            const bool erode,
            labelList& closureFaces,
            labelList& closureToBlocked
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
