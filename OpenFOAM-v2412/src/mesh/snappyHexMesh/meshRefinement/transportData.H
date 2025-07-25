/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::transportData

Description
    Holds information (coordinate and distance). Walks out 0.5*distance.

SourceFiles
    transportDataI.H
    transportData.C

\*---------------------------------------------------------------------------*/

#ifndef transportData_H
#define transportData_H

#include "wallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class transportData Declaration
\*---------------------------------------------------------------------------*/

class transportData
:
    public wallPointData<scalar>
{
public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Per face the index of the surface hit
        const labelList& surfaceIndex_;

        trackData(const labelList& surfaceIndex)
        :
            surfaceIndex_(surfaceIndex)
        {}
    };


private:

    // Private Member Functions

        //- Evaluate distance to point.
        //  Update distSqr, origin from whomever is nearer pt.
        //  \return true if w2 is closer to point, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const transportData& w2,
            const scalar tol,
            TrackingData& td
        );

public:

    // Constructors

        //- Construct null
        inline transportData();

        //- Construct from origin, gapSize, distance
        inline transportData
        (
            const point& origin,
            const scalar gapSize,
            const scalar distSqr
        );


    // Member Functions

        // Needed by FaceCellWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of cellI
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCellI,
                const label neighbourFaceI,
                const transportData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const label neighbourCellI,
                const transportData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of faceI
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFaceI,
                const transportData& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for transportData
template<>
struct is_contiguous<transportData> : is_contiguous<wallPointData<scalar>> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transportDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
