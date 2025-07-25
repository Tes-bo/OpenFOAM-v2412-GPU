/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::shellSurfaces

Description
    Encapsulates queries for volume refinement ('refine all cells within
    shell').

SourceFiles
    shellSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef shellSurfaces_H
#define shellSurfaces_H

#include "searchableSurface.H"
#include "Enum.H"
#include "Tuple2.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                           Class shellSurfaces Declaration
\*---------------------------------------------------------------------------*/

class shellSurfaces
{
public:

    // Public data types

        //- Volume refinement controls
        enum refineMode
        {
            INSIDE,         // Refine all inside shell
            OUTSIDE,        // ,,         outside
            DISTANCE        // Refine based on distance to shell
        };


private:

    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Dry-run mode?
        const bool dryRun_;

        //- Indices of surfaces that are shells
        labelList shells_;

        //- Per shell whether to refine inside or outside
        List<refineMode> modes_;

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;


        // Directional

            //- Per shell any additional directional refinement
            List<Tuple2<labelPair,labelVector>> dirLevels_;

            //- Per shell the smoothing direction
            vectorField smoothDirection_;

            //- Per shell the directional smoothing iterations
            labelList nSmoothExpansion_;

            //- Per shell the positional smoothing iterations
            labelList nSmoothPosition_;


        // Gap level refinement

            //- Per shell, per region the small-gap level specification
            List<List<FixedList<label, 3>>> extendedGapLevel_;

            //- Per shell, per region the small-gap level specification
            List<List<volumeType>> extendedGapMode_;

            //- Per shell, per region whether to allow selfProximity refinement
            boolListList selfProximity_;


    // Private data

        //- refineMode names
        static const Enum<refineMode> refineModeNames_;


    // Private Member Functions

        //- Helper function for initialisation of levels
        void setAndCheckLevels
        (
            const label shellI,
            const List<Tuple2<scalar, label>>&
        );

        void orient();

        //- Find first shell with a level higher than maxLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& maxLevel
        ) const;

        //- Update highest min gap level
        void findHigherGapLevel
        (
            const pointField& pt,
            const labelList& ptLevel,
            const label shellI,
            labelList& gapShell,
            List<FixedList<label, 3>>& gapInfo,
            List<volumeType>& gapMode
        ) const;

        //- Find first shell with a level lower or equal to minLevel. Update
        //  minLevel and shell.
        void findLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& minLevel,
            labelList& shell
        ) const;

public:

    // Constructors

        //- Construct from geometry and dictionary
        shellSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary& shellsDict,
            const bool dryRun = false
        );


    // Member Functions

        // Access

            //- Indices of surfaces that are shells
            const labelList& shells() const
            {
                return shells_;
            }

            //- Per shell, per region the small-gap level specification
            const List<List<FixedList<label, 3>>>& extendedGapLevel() const
            {
                return extendedGapLevel_;
            }

            //- Per shell, per region the small-gap level specification
            const List<List<volumeType>>& extendedGapMode() const
            {
                return extendedGapMode_;
            }

            //- Per shell, per region whether to test for gap with same surface
            const boolListList& gapSelf() const
            {
                return selfProximity_;
            }


        // Query

            //- Highest shell level
            label maxLevel() const;

            //- Highest shell gap level
            labelList maxGapLevel() const;

            //- Min and max cell level for directional refinement
            labelPairList directionalSelectLevel() const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& gapShell,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel. gapInfo
            //  is (0 0 0) if no shell found
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find first shell (or -1) with level equal or lower than ptLevel.
            void findLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& shell
            ) const;

            //- Find any shell (or -1) with higher wanted directional level
            void findDirectionalLevel
            (
                const pointField& pt,
                const labelList& ptLevel,   // omnidirectional level
                const labelList& dirLevel,  // directional level
                const direction dir,
                labelList& shell
            ) const;

            //- Per shell the smoothing direction
            const vectorField& smoothDirection() const;

            //- Per shell the directional smoothing iterations
            const labelList& nSmoothExpansion() const;

            //- Per shell the positional smoothing iterations
            const labelList& nSmoothPosition() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
