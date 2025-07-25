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
    Foam::equipotentialPointSmoother

Description
    Equipotential point smoothing. Points are moved towards the centroid of the
    surrounding cell volumes. This method tends to equilise cell volumes, and
    can generate distorted cells around refinement patterns.

SourceFiles
    equipotentialPointSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef equipotentialPointSmoother_H
#define equipotentialPointSmoother_H

#include "pointSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointSmoothers
{

/*---------------------------------------------------------------------------*\
                   Class equipotentialPointSmoother Declaration
\*---------------------------------------------------------------------------*/

class equipotentialPointSmoother
:
    public pointSmoother
{
public:

    //- Runtime type information
    TypeName("equipotential");


    // Constructors

        //- Construct from a dictionary and a polyMesh
        equipotentialPointSmoother
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~equipotentialPointSmoother() = default;


    // Member Functions

        //- Calculate the point displacements
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
