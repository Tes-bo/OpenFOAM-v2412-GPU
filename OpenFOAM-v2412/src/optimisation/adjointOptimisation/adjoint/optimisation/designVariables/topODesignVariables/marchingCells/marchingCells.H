/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    marchingCells

Description
    Marches cells from given seed patches, faceZones or the boundary faces of
    cellZones, towards the interior of the mesh, with a given step at a time.
    Used in topology optimisation, to gradually activate the design variables,
    in an attempt to prevent the optimiser from acting too agressively in the
    first optimisation cycles. Use meshWave to do the marching.

SourceFiles
    marchingCells.C

\*---------------------------------------------------------------------------*/

#ifndef marchingCells_H
#define marchingCells_H

#include "fvMesh.H"
#include "wallPointData.H"
#include "FaceCellWave.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class marchingCells Declaration
\*---------------------------------------------------------------------------*/

class marchingCells
{

protected:

    // Protected data

        const fvMesh& mesh_;

        //- Patches used as seeds in the marching algorithm
        labelHashSet seedPatches_;

        //- Cell zones, the boundary faces of which are  used as seeds in the
        //- marching algorithm
        labelList seedCellZoneIDs_;

        //- Face zones used as seeds in the marching algorithm.
        labelList seedFaceZoneIDs_;

        //- Marching step
        label marchingStep_;

        //- Whether each cell is curently active or not
        boolList isActiveCell_;

        //- Should this cell remain incative
        boolList isFixedCell_;

        //- Which are the active cells
        DynamicList<label> activeCells_;

        //- Which are the added cells
        labelList addedCells_;

        //- Has the initial seeding been conducted
        bool initialised_;

        //- Iterations conducted thus far
        label nIters_;

        //- Information for all faces
        List<wallPointData<bool>> allFaceInfo_;

        //- Information for all cells
        List<wallPointData<bool>> allCellInfo_;

        //- Engine propagating the active cells.
        //  Uses infrastructure from FaceCellWave, to easily handle parallel
        //  interfaces
        FaceCellWave<wallPointData<bool>> meshWave_;


    // Protected Member Functions

        //- Initialise the active cells from the seeding patches
        void initialise();

        //- Append cell to seed cells
        void appendSeedCell(const label cellID);

        // Updates the isActiveCells_ list and gathers newlyAddedCells
        void march
        (
            label nVisited,
            const label cI,
            labelList& newlyAddedCells
        );


private:

    // Private Member Functions

        //- No copy construct
        marchingCells(const marchingCells&) = delete;

        //- No copy assignment
        void operator=(const marchingCells&) = delete;


public:

    //- Runtime type information
    TypeName("marchingCells");


    // Constructors

        //- Construct from components
        marchingCells
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~marchingCells() = default;


    // Member Functions

        // Access

            //- Get the active cells
            inline const labelList& getActiveCells() const
            {
                return activeCells_;
            }

            //- Number of newly added cells in the last iteration
            inline label numberOfAddedCells() const
            {
                return addedCells_.size();
            }


        // Edit

            //- Update active cells
            void update(const label iters = 1);

            //- Set marching step
            inline void setMarchingStep(const label step)
            {
                marchingStep_ = step;
            }

            //- Add fixed cells through cellZone IDs
            void addFixedCells
            (
                const cellZoneMesh& cellZoneMesh,
                const labelList& fixedCellZoneIDs
            );

            //- Add fixed cells
            void addFixedCells(const labelList& fixedCellIDs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
