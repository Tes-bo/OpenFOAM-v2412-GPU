/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::MeshWave

Description
    FaceCellWave plus data

SourceFiles
    MeshWave.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_MeshWave_H
#define Foam_MeshWave_H

#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class MeshWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(MeshWave);


/*---------------------------------------------------------------------------*\
                           Class MeshWave Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class TrackingData = int>
class MeshWave
:
    public MeshWaveName
{
    // Private Data

        //- Wall information for all faces
        List<Type> allFaceInfo_;

        //- Wall information for all cells
        List<Type> allCellInfo_;

        //- Wave calculation engine.
        FaceCellWave<Type, TrackingData> calc_;


    // Private Member Functions

        //- No copy construct
        MeshWave(const MeshWave&) = delete;

        //- No copy assignment
        void operator=(const MeshWave&) = delete;


public:

    // Constructors

        //- Construct from mesh and list of changed faces with the Type
        //- for these faces.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        MeshWave
        (
            const polyMesh& mesh,
            const labelList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            const label maxIter,
            TrackingData& td = FaceCellWaveBase::dummyTrackData_
        );

        //- Construct from mesh, list of changed faces with the Type
        //- for these faces and initial field.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        MeshWave
        (
            const polyMesh& mesh,
            const labelList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            const List<Type>& allCellInfo,
            const label maxIter,
            TrackingData& td = FaceCellWaveBase::dummyTrackData_
        );


    // Member Functions

        //- Get allFaceInfo
        const List<Type>& allFaceInfo() const noexcept
        {
            return allFaceInfo_;
        }

        //- Get allCellInfo
        const List<Type>& allCellInfo() const noexcept
        {
            return allCellInfo_;
        }

        //- Additional data to be passed into container
        const TrackingData& data() const noexcept
        {
            return calc_.data();
        }

        //- Iterate until no changes or maxIter reached.
        //  \return actual number of iterations.
        label iterate(const label maxIter)
        {
            return calc_.iterate(maxIter);
        }

        //- Current number of changed cells
        label nChangedCells() const noexcept { return calc_.nChangedCells(); }

        //- Current number of changed faces
        label nChangedFaces() const noexcept { return calc_.nChangedFaces(); }

        //- Number of unvisited cells
        label nUnvisitedCells() const noexcept
        {
            return calc_.nUnvisitedCells();
        }

        //- Number of unvisited faces
        label nUnvisitedFaces() const noexcept
        {
            return calc_.nUnvisitedFaces();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
