/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::thresholdCellFaces

Description
    Selects the mesh cell faces specified by a threshold value.
    Non-triangulated by default.

SourceFiles
    thresholdCellFaces.C

\*---------------------------------------------------------------------------*/

#ifndef thresholdCellFaces_H
#define thresholdCellFaces_H

#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class thresholdCellFaces Declaration
\*---------------------------------------------------------------------------*/

class thresholdCellFaces
:
    public meshedSurface
{
    // Private typedefs for convenience
    typedef meshedSurface Mesh;


    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- For each face, the original cell in mesh
        labelList meshCells_;


    // Private Member Functions

        void calculate
        (
            const scalarField& field,
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate
        );

public:

    //- Runtime type information
    TypeName("thresholdCellFaces");


    // Constructors

        //- Construct from mesh, field and threshold values
        thresholdCellFaces
        (
            const polyMesh& mesh,
            const scalarField& field,
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate = false
        );


    // Member Functions

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            return meshCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
