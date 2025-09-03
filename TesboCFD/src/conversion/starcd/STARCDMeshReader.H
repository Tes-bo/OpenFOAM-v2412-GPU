/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileFormats::STARCDMeshReader

Description
    Read PROSTAR vrt/cel/bnd files.
    The protected data in meshReader are filled.

    Starting with PROSTAR version 4, the files have become easier to read.
    - vertices are space-delimited.
    - the cell format is logical.
    - trimmed and degenerate cells are saved as polyhedral.
    - the boundaries corresponds to cells and their faces.

Notes
    - Shell entries are read, yet ignored.

SourceFiles
    STARCDMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshReader_H
#define STARCDMeshReader_H

#include "meshReader.H"
#include "STARCDCore.H"
#include "boundaryRegion.H"
#include "cellShape.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::STARCDMeshReader Declaration
\*---------------------------------------------------------------------------*/

class STARCDMeshReader
:
    public meshReader,
    protected STARCDCore
{
    // Private data

        //- Retain solid cell types
        bool keepSolids_;


protected:

    // Protected Data

        //- Cell shapes
        cellShapeList cellShapes_;

        //- Point labels (imported Point numbering not necessarily contiguous)
        labelList mapToFoamPointId_;

        //- Cell labels (imported Cell numbering not necessarily contiguous)
        labelList mapToFoamCellId_;

        //- Boundary region data
        boundaryRegion boundaryRegion_;


    // Protected Member Functions

        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);

        //- Read points from file, return the max prostar id used.
        label readPoints(const fileName&, const scalar scaleFactor);

        //- Read cell connectivities from file
        virtual void readCells(const fileName&);

        //- Remove unused points
        void cullPoints();

        //- Read boundary (cell/face) definitions
        void readBoundary(const fileName&);

        //- Read auxiliary data from constant/{boundaryRegion,cellTable}
        void readAux(const objectRegistry&);


public:

    // Constructors

        //- Construct from case name
        STARCDMeshReader
        (
            const fileName& prefix,
            const objectRegistry& registry,
            const scalar scaleFactor = 1.0,
            const bool keepSolids = false
        );


    //- Destructor
    virtual ~STARCDMeshReader() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
