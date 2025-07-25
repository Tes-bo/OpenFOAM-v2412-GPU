/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::laplacianConstraintPointSmoother

Description
    Laplacian point smoothing. Points are moved towards the average of the
    surrounding points. In case of constraints (boundary points) only
    use average of equal or higher constrained points.

SourceFiles
    laplacianConstraintPointSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef laplacianConstraintPointSmoother_H
#define laplacianConstraintPointSmoother_H

#include "pointSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointSmoothers
{

/*---------------------------------------------------------------------------*\
              Class laplacianConstraintPointSmoother Declaration
\*---------------------------------------------------------------------------*/

class laplacianConstraintPointSmoother
:
    public pointSmoother
{
public:

    //- Runtime type information
    TypeName("laplacianConstraint");


    // Constructors

        //- Construct from a dictionary and a polyMesh
        laplacianConstraintPointSmoother
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~laplacianConstraintPointSmoother() = default;


    // Member Functions

        //- Calculate the point displacement
        virtual void calculate
        (
            const labelList& facesToMove,
            const pointField& oldPoints,
            const pointField& currentPoints,
            const pointField& faceCentres,
            const vectorField& faceAreas,
            const pointField& cellCentres,
            const scalarField& cellVolumes,
            vectorField& pointDisplacement
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointSmoothers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
