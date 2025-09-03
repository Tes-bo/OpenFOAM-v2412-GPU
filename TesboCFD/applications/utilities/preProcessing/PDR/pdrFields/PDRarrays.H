/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::PDRarrays

Description
    Work array definitions for PDR fields

SourceFiles
    PDRarrays.C
    PDRarraysCalc.C

\*---------------------------------------------------------------------------*/

#ifndef PDRarrays_H
#define PDRarrays_H

#include "symmTensor.H"
#include "symmTensor2D.H"
#include "SquareMatrix.H"
#include "IjkField.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class PDRblock;
class PDRmeshArrays;
class PDRobstacle;
class PDRpatchDef;


/*---------------------------------------------------------------------------*\
                          Class PDRarrays Declaration
\*---------------------------------------------------------------------------*/

class PDRarrays
{
    //- Reference to PDRblock
    std::reference_wrapper<const PDRblock> pdrBlock_;

public:

    // Data Members
    // Entries used for analysis and when writing fields

    //- Volume blockage
    IjkField<scalar> v_block;

    //- Surface area in cell
    IjkField<scalar> surf;

    //- Obstacle size in cell
    IjkField<scalar> obs_size;

    //- Summed area blockage (directional) from sharp obstacles
    IjkField<vector> area_block_s;

    //- Summed area blockage (directional) from round obstacles
    IjkField<vector> area_block_r;

    //- A total directional blockage in the cell
    IjkField<Vector<bool>> dirn_block;

    //- Face area blockage for face,
    //- summed from cell centre-plane to cell centre-plane
    IjkField<vector> face_block;

    //- Longitudinal area blockage from obstacles that extend all the way
    //- through the cell in a given direction.
    IjkField<vector> along_block;

    IjkField<vector> betai_inv1;

    //- Number of obstacles in cell.
    //  Can be non-integer if an obstacle does not pass all way through cell
    IjkField<scalar> obs_count;

    //- Number of obstacles parallel to specified direction
    IjkField<vector> sub_count;

    //- Addition to count to account for grating comprises many bars
    //- (to get Lobs right)
    IjkField<vector> grating_count;

    //- Tensorial drag from sharp obstacles
    IjkField<symmTensor> drag_s;

    //- Directional drag from round obstacles
    IjkField<vector> drag_r;


    // Next arrays are for 2D calculations of intersection

    // One-dimensional scratch areas for cell overlaps
    Vector<List<scalar>> overlap_1d;

    // In two dimensions, area of cell covered by circle
    SquareMatrix<scalar> aboverlap;

    // In two dimensions, length of perimeter of circle witthin cell
    SquareMatrix<scalar> abperim;

    // For offset cells, i.e. face blockage
    SquareMatrix<scalar> a_lblock, b_lblock;

    // For centred cells
    SquareMatrix<scalar> ac_lblock, bc_lblock;

    // The count in the cells
    SquareMatrix<scalar> c_count;

    //- Cell-centred drag
    SquareMatrix<symmTensor2D> c_drag;

    //- Face field for (directional) for patch Id
    IjkField<labelVector> face_patch;

    //- Face field for (directional) hole in face
    IjkField<Vector<bool>> hole_in_face;


    // Constructors

        //- Construct null
        PDRarrays();

        //- Construct and reset
        explicit PDRarrays(const PDRblock& pdrBlock);


    //- Destructor
    ~PDRarrays() = default;


    // Member Functions

        //- Reset PDRblock reference, resize and zero arrays
        void reset(const PDRblock& pdrBlock);

        //- Reference to PDRblock
        const PDRblock& block() const
        {
            return pdrBlock_.get();
        }

        //- Summary of the blockages
        //  For diagnostics and general overview
        void blockageSummary() const;

        //- Add cylinder blockage
        void addCylinder(const PDRobstacle& obs);

        //- Add general (non-cylinder) blockage
        void addBlockage
        (
            const PDRobstacle& obs,
            DynamicList<PDRpatchDef>& patches,
            const int volumeSign
        );


        static void calculateAndWrite
        (
            PDRarrays& arr,
            const PDRmeshArrays& meshIndexing,
            const fileName& casepath,
            const UList<PDRpatchDef>& patches
        );

        void calculateAndWrite
        (
            const fileName& casepath,
            const PDRmeshArrays& meshIndexing,
            const UList<PDRpatchDef>& patches
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
