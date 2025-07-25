/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::edgeInterpolation

Description
    Face to edge interpolation scheme. Included in faMesh.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    edgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeInterpolation_H
#define Foam_edgeInterpolation_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class edgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

class edgeInterpolation
{
    // Private Data

        // Reference to faMesh
        const faMesh& faMesh_;


        // Demand-driven data

            //- Geodesic distance between centroids of neighbour finite areas
            mutable std::unique_ptr<edgeScalarField> lPNptr_;

            //- Central-differencing weighting factors
            mutable std::unique_ptr<edgeScalarField> weightingFactorsPtr_;

            //- Face-gradient difference factors
            mutable std::unique_ptr<edgeScalarField> differenceFactorsPtr_;

            //- Non-orthogonality correction vectors
            mutable std::unique_ptr<edgeVectorField> correctionVectorsPtr_;

            //- Skew correction vectors
            mutable std::unique_ptr<edgeVectorField> skewCorrectionVectorsPtr_;

            //- Is mesh orthogonal
            mutable bool orthogonal_;

            //- Is mesh skew
            mutable bool skew_;


    // Private Member Functions

        //- Return skewness correction per edge
        const vector& skewCorr(const label edgeI) const;

        //- Construct geodesic distance between P and N
        void makeLPN() const;

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeCorrectionVectors() const;

        //- Construct skewness correction vectors
        void makeSkewCorrectionVectors() const;


public:

    // Declare name of the class and it's debug switch
    ClassName("edgeInterpolation");


    // Generated Methods

        //- No copy construct
        edgeInterpolation(const edgeInterpolation&) = delete;

        //- No copy assignment
        void operator=(const edgeInterpolation&) = delete;


    // Constructors

        //- Construct given an faMesh
        explicit edgeInterpolation(const faMesh&);


    //- Destructor
    ~edgeInterpolation() = default;


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const noexcept { return faMesh_; }

        //- Return reference to PN geodesic distance
        const edgeScalarField& lPN() const;

        //- Return reference to weighting factors array
        const edgeScalarField& weights() const;

        //- Return reference to difference factors array
        const edgeScalarField& deltaCoeffs() const;

        //- Return reference to non-orthogonality correction vectors array
        const edgeVectorField& correctionVectors() const;

        //- Return reference to skew vectors array
        const edgeVectorField& skewCorrectionVectors() const;

        //- Return whether mesh is orthogonal or not
        bool orthogonal() const;

        //- Return whether mesh is skew or not
        bool skew() const;


    // Mesh Motion

        //- Do what is necessary if the mesh has moved
        bool movePoints() const;


    // Storage Management

        //- True if weights exist
        bool hasWeights() const noexcept { return bool(weightingFactorsPtr_); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
