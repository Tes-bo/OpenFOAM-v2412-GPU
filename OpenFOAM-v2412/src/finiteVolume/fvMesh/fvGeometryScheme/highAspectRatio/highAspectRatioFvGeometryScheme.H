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
    Foam::highAspectRatioFvGeometryScheme

Description
    Geometry calculation scheme with automatic stabilisation for high-aspect
    ratio cells.

SourceFiles
    highAspectRatioFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef highAspectRatioFvGeometryScheme_H
#define highAspectRatioFvGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
               Class highAspectRatioFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class highAspectRatioFvGeometryScheme
:
    public basicFvGeometryScheme
{
protected:

        const scalar minAspect_;

        const scalar maxAspect_;


    // Protected Member Functions

        //- Calculate cell and face weight. Is 0 for cell < minAspect, 1 for
        //  cell > maxAspect
        void calcAspectRatioWeights
        (
            scalarField& cellWeight,
            scalarField& faceWeight
        ) const;

        //- Helper : calculate (weighted) average face and cell centres
        static void makeAverageCentres
        (
            const polyMesh& mesh,
            const pointField& points,
            const pointField& faceAreas,
            const scalarField& magFaceAreas,
            pointField& faceCentres,
            pointField& cellCentres
        );


private:

    // Private Member Functions

        //- No copy construct
        highAspectRatioFvGeometryScheme
        (
            const highAspectRatioFvGeometryScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const highAspectRatioFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("highAspectRatio");


    // Constructors

        //- Construct from mesh
        highAspectRatioFvGeometryScheme
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~highAspectRatioFvGeometryScheme() = default;


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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
