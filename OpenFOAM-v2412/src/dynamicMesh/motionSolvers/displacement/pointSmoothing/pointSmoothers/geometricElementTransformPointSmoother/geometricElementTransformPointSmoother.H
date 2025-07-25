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
    Foam::pointSmoothers::geometricElementTransformPointSmoother

Description
    Geometric Element Transformation Method (GETMe) point smoother. Points are
    moved in a direction normal to the face of the dual element. This tends to
    orthogonalise elements. This method can "push" as well as "pull".

    The original references for this method only formulate it for a limited
    variety of specific polyhedra. This version is written for arbitrary
    polyhedra, and most closely follows the following reference:

    \verbatim
        "Fast smoothing of mixed volume meshes based on the effective geometric
        element transformation method"
        Dimitris Vartziotis, Joachim Wipper,
        Comput. Methods Appl. Mech. Engrg. 201-204 (2012) 65-81
    \endverbatim

    This implementation does not include the various specific measures employed
    in the referred article to improve the quality of types of irregular
    polyhedra. It also does not use the same quality criteria as these are only
    suitable for cell vertices with exactly three connecting edges.

SourceFiles
    geometricElementTransformPointSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef geometricElementTransformPointSmoother_H
#define geometricElementTransformPointSmoother_H

#include "pointSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellPointConnectivity;

namespace pointSmoothers
{

/*---------------------------------------------------------------------------*\
             Class geometricElementTransformPointSmoother Declaration
\*---------------------------------------------------------------------------*/

class geometricElementTransformPointSmoother
:
    public pointSmoother
{
private:

    // Private Data

        //- Transformation parameter
        const scalar transformationParameter_;


    // Private Member Functions

        static scalar cellQuality
        (
            const polyMesh& mesh,
            const pointField& currentPoints,
            const cellPointConnectivity& connectivity,
            const label celli
        );


public:

    //- Runtime type information
    TypeName("geometricElementTransform");


    // Constructors

        //- Construct from a dictionary and a mesh
        geometricElementTransformPointSmoother
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~geometricElementTransformPointSmoother() = default;


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

        static tmp<scalarField> cellQuality
        (
            const polyMesh& mesh,
            const pointField& currentPoints
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointSmoothers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
