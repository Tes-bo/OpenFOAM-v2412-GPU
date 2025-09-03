/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::PDRmeshArrays

Description
    OpenFOAM/PDRblock addressing information

    Provides mapping for a rectilinear OpenFOAM mesh in terms of
    i-j-k indices for faces and cells.

    The mesh points are first binned according to their i-j-k locations.
    Next the faces are classified according to their lowest x/y/z
    coordinates and the face orientation as x/y/z.
    Orientation in the sense +x or -x is not noted.
    The cell faces are then examined to determine the appropriate i-j-k
    location.

SourceFiles
    PDRmeshmeshArraysIO.C

\*---------------------------------------------------------------------------*/

#ifndef PDRmeshArrays_H
#define PDRmeshArrays_H

#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class PDRblock;
class polyMesh;
class Time;

/*---------------------------------------------------------------------------*\
                          Class PDRmeshArrays Declaration
\*---------------------------------------------------------------------------*/

class PDRmeshArrays
{
public:

    //- Relative tolerance when matching grid points. Default = 0.02
    static scalar gridPointRelTol;

    //- The cell i-j-k addressing range
    labelVector cellDims;

    //- The face i-j-k addressing range
    labelVector faceDims;

    //- For each cell, the corresponding i-j-k address.
    List<labelVector> cellIndex;

    //- For each face, the corresponding i-j-k address.
    List<labelVector> faceIndex;

    //- For each face, the x/y/z orientation
    List<direction> faceOrient;


    // Constructors

        //- Default construct
        PDRmeshArrays() = default;


    //- Destructor
    ~PDRmeshArrays() = default;


    // Member Functions

        //- The number of cells
        label nCells() const
        {
            return cellIndex.size();
        }

        //- The number of faces
        label nFaces() const
        {
            return faceIndex.size();
        }


        //- Determine i-j-k indices for faces/cells
        void classify(const polyMesh& mesh, const PDRblock& pdrBlock);

        //- Read OpenFOAM mesh and determine i-j-k indices for faces/cells
        void read(const Time& runTime, const PDRblock& pdrBlock);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
