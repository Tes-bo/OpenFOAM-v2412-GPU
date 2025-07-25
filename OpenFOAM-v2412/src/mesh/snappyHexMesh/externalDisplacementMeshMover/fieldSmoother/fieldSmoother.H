/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::fieldSmoother

Description
    Utility functions for mesh motion solvers

SourceFiles
    fieldSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef fieldSmoother_H
#define fieldSmoother_H

#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldSmoother Declaration
\*---------------------------------------------------------------------------*/

class fieldSmoother
{
    // Private data

        //- Reference to the poly mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        fieldSmoother(const fieldSmoother&) = delete;

        //- No copy assignment
        void operator=(const fieldSmoother&) = delete;


public:

    // Run-time type information
    TypeName("fieldSmoother");


    // Constructors

        //- Construct from a polyMesh
        fieldSmoother(const polyMesh&);


    //- Destructor
    virtual ~fieldSmoother();


    // Member Functions

        //- Smooth interior normals
        void smoothNormals
        (
            const label nIter,
            const bitSet& isMeshMasterPoint,
            const bitSet& isMeshMasterEdge,
            const labelList& fixedPoints,
            pointVectorField& normals
        ) const;

        //- Smooth patch normals
        void smoothPatchNormals
        (
            const label nIter,
            const bitSet& isPatchMasterPoint,
            const bitSet& isPatchMasterEdge,
            const indirectPrimitivePatch& adaptPatch,
            pointField& normals
        ) const;

        //- Smooth a scalar field towards, but not beyond, a minimum value
        template <class Type>
        void minSmoothField
        (
            const label nIter,
            const bitSet& isPatchMasterPoint,
            const bitSet& isPatchMasterEdge,
            const indirectPrimitivePatch& adaptPatch,
            const scalarField& fieldMin,
            Field<Type>& field
        ) const;

        //- Smooth and then un-smooth a displacement
        void smoothLambdaMuDisplacement
        (
            const label nIter,
            const bitSet& isMeshMasterPoint,
            const bitSet& isMeshMasterEdge,
            const bitSet& isSmoothable,
            vectorField& displacement
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
