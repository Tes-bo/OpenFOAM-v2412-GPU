/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceInterpolation_H
#define Foam_surfaceInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class fvGeometryScheme;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
{
    // Private Data

        // Reference to fvMesh
        const fvMesh& mesh_;


        // Demand-driven data

            //- Geometry calculation
            mutable tmp<fvGeometryScheme> geometryPtr_;

            //- Linear difference weighting factors
            mutable autoPtr<surfaceScalarField> weights_;

            //- Cell-centre difference coefficients
            mutable autoPtr<surfaceScalarField> deltaCoeffs_;

            //- Non-orthogonal cell-centre difference coefficients
            mutable autoPtr<surfaceScalarField> nonOrthDeltaCoeffs_;

            //- Non-orthogonality correction vectors
            mutable autoPtr<surfaceVectorField> nonOrthCorrectionVectors_;


protected:

    // Protected Member Functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and its debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    //- Destructor
    virtual ~surfaceInterpolation();


    // Member Functions

        //- Return reference to geometry calculation scheme
        virtual const fvGeometryScheme& geometry() const;

        //- Set geometry calculation scheme
        void geometry(tmp<fvGeometryScheme>&);

        //- Return reference to linear difference weighting factors
        virtual const surfaceScalarField& weights() const;

        //- Return reference to cell-centre difference coefficients
        virtual const surfaceScalarField& deltaCoeffs() const;

        //- Return reference to non-orthogonal cell-centre difference
        //  coefficients
        virtual const surfaceScalarField& nonOrthDeltaCoeffs() const;

        //- Return reference to non-orthogonality correction vectors
        virtual const surfaceVectorField& nonOrthCorrectionVectors() const;

        //- Do what is necessary if the mesh has moved
        virtual bool movePoints();

        //- Update all geometric data
        virtual void updateGeom();


    // Storage Management

        //- Has weights
        bool hasWeights() const noexcept { return bool(weights_); }

        //- Update mesh for topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
