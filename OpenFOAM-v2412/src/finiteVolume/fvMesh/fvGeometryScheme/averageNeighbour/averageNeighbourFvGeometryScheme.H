/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::averageNeighbourFvGeometryScheme

Description
    Geometry calculation scheme to minimise non-orthogonality/

SourceFiles
    averageNeighbourFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef averageNeighbourFvGeometryScheme_H
#define averageNeighbourFvGeometryScheme_H

#include "highAspectRatioFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;

/*---------------------------------------------------------------------------*\
               Class averageNeighbourFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class averageNeighbourFvGeometryScheme
:
    public highAspectRatioFvGeometryScheme
{
private:

        //- No copy construct
        averageNeighbourFvGeometryScheme
        (
            const averageNeighbourFvGeometryScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const averageNeighbourFvGeometryScheme&) = delete;


protected:

        //- Number of averaging iterations
        const label nIters_;

        //- Blending between old-iteration cell centres and current average
        const scalar relax_;

        //- Clipping for pyramid heights - allowable shrinkage as fraction
        //  of original
        const scalar minRatio_;


        //- Clip face-centre correction vector if new triangle area
        //- would get below min. Return number of clipped faces.
        label clipFaceTet
        (
            const scalar minRatio,
            const vectorField& faceCentres,
            const vectorField& faceNormals,
            vectorField& faceCorrection
        ) const;

        //- Calculate pyramid heights
        void makePyrHeights
        (
            const pointField& cellCentres,
            const vectorField& faceCentres,
            const vectorField& faceNormals,

            scalarField& ownHeight,
            scalarField& neiHeight
        ) const;

        //- Clip correction vector if new pyramid height would get below min.
        //- Return number of clipped cells.
        label clipPyramids
        (
            const pointField& cellCentres,
            const vectorField& faceCentres,
            const vectorField& faceNormals,

            const scalarField& minOwnHeight,
            const scalarField& minNeiHeight,

            vectorField& correction
        ) const;

        //- Average neighbouring cell centres to minimise non-ortho
        tmp<pointField> averageNeighbourCentres
        (
            const pointField& cellCentres,
            const vectorField& faceNormals,
            const scalarField& faceWeights
        ) const;

        tmp<pointField> averageCentres
        (
            const pointField& cellCentres,
            const pointField& faceCentres,
            const vectorField& faceNormals
        ) const;

        //- Make weights based on non-orthogonality
        void makeNonOrthoWeights
        (
            const pointField& cellCentres,
            const vectorField& faceNormals,

            scalarField& cosAngles,
            scalarField& faceWeights
        ) const;


public:

    //- Runtime type information
    TypeName("averageNeighbour");


    // Constructors

        //- Construct from mesh
        averageNeighbourFvGeometryScheme
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~averageNeighbourFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Calculate geometry quantities using mesh topology and provided
        //- points. If oldPoints provided only does local update. Returns
        //- true if anything changed, false otherwise
        virtual bool updateGeom
        (
            const pointField& points,
            const refPtr<pointField>& oldPoints,    // optional old points
            pointField& faceCentres,
            vectorField& faceAreas,
            pointField& cellCentres,
            scalarField& cellVolumes
        ) const
        {
            NotImplemented;
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
