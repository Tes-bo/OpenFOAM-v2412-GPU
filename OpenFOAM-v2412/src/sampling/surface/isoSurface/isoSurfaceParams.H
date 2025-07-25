/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::isoSurfaceParams

Description
    Preferences for controlling iso-surface algorithms.

    Some common dictionary properties:
    \table
        Property    | Description                           | Required | Default
        isoMethod   | Algorithm (cell/topo/point/default)   | no  | default
        regularise  | Face simplification (enum or bool)    | no  | true
        mergeTol    | Point merge tolerance (cell/point)    | no  | 1e-6
        snap        | Point snapping (topo)                 | no  | true
        bounds      | Optional clip bounds                  | no  | inverted
    \endtable

    The default algorithm denotes the use of the current \em standard
    algorithm.

Filtering types (for topological iso-surface)
    - \c none : leave tet cuts untouched
    - \c partial , \c cell : Combine intra-cell faces
    - \c full , \c diagcell : Perform \c partial and remove face-diagonal
      points
    - \c clean : Perform \c full and eliminate open edges as well.
      (<b>May cause excessive erosion!</b>)
    .

SourceFiles
    isoSurfaceParams.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_isoSurfaceParams_H
#define Foam_isoSurfaceParams_H

#include "boundBox.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class isoSurfaceSelector Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceParams
{
public:

    // Data Types

        //- The algorithm types
        enum algorithmType : uint8_t
        {
            ALGO_DEFAULT = 0,   //!< Use current 'standard' algorithm
            ALGO_TOPO,
            ALGO_CELL,
            ALGO_POINT
        };

        //- The filtering (regularization) to apply
        enum class filterType : uint8_t
        {
            NONE = 0,   //!< No filtering
            CELL,       //!< Remove pyramid edge points
            DIAGCELL,   //!< Remove pyramid edge points, face-diagonals
            NONMANIFOLD, //!< Remove pyramid edge points, face-diagonals
                         //!< and non-manifold faces

            PARTIAL = CELL,     //!< Same as CELL
            FULL = DIAGCELL,    //!< Same as DIAGCELL
            CLEAN = NONMANIFOLD //!< Same as NONMANIFOLD
        };


private:

    // Private Data

        //- Algorithm type
        algorithmType algo_;

        //- Filtering for iso-surface faces/points
        filterType filter_;

        //- Point snapping enabled
        bool snap_;

        //- Merge tolerance for cell/point (default: 1e-6)
        scalar mergeTol_;

        //- Optional bounding box for clipping (default: inverted)
        boundBox clipBounds_;


public:

    // Public Data

        //- Names for the iso-surface algorithms
        static const Enum<algorithmType> algorithmNames;

        //- Names for the filtering types
        static const Enum<filterType> filterNames;


    // Static Member Functions

        //- Get 'isoMethod' or 'isoAlgorithm' as enumeration
        static algorithmType getAlgorithmType
        (
            const dictionary& dict,
            const algorithmType deflt
        );

        //- Get 'regularise' as bool or enumeration
        static filterType getFilterType
        (
            const dictionary& dict,
            const filterType deflt
        );


    // Constructors

        //- Default construct, or with specified algorithm
        explicit isoSurfaceParams
        (
            const algorithmType algo = algorithmType::ALGO_DEFAULT,
            const filterType filter = filterType::DIAGCELL
        ) noexcept;

        //- Default construct, setting parameters from dictionary
        explicit isoSurfaceParams
        (
            const dictionary& dict,
            const isoSurfaceParams& params = isoSurfaceParams()
        );

        //- Default construct, setting parameters from dictionary
        explicit isoSurfaceParams
        (
            const dictionary& dict,
            const algorithmType algo,
            const filterType filter = filterType::DIAGCELL
        );


    // Member Functions

        //- Get current algorithm
        algorithmType algorithm() const noexcept
        {
            return algo_;
        }

        //- Set algorithm
        void algorithm(algorithmType algo) noexcept
        {
            algo_ = algo;
        }

        //- Get current filter type
        filterType filter() const noexcept
        {
            return filter_;
        }

        //- Set filter type
        void filter(filterType fltr) noexcept
        {
            filter_ = fltr;
        }

        //- Get point snapping flag
        bool snap() const noexcept
        {
            return snap_;
        }

        //- Set point snapping flag
        void snap(bool on) noexcept
        {
            snap_ = on;
        }

        //- Get current merge tolerance
        scalar mergeTol() const noexcept
        {
            return mergeTol_;
        }

        //- Set merge tolerance (cell/point algo)
        void mergeTol(const scalar relTol) noexcept
        {
            mergeTol_ = relTol;
        }

        //- Get optional clipping bounding box
        const boundBox& getClipBounds() const noexcept
        {
            return clipBounds_;
        }

        //- Access optional clipping bounding box
        boundBox& getClipBounds() noexcept
        {
            return clipBounds_;
        }

        //- Set optional clipping bounding box
        void setClipBounds(const boundBox& bb);


    // Information

        //- Print information about the settings
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
