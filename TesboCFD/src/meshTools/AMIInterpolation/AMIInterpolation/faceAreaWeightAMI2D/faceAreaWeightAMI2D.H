/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::faceAreaWeightAMI2D

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method that performs
    the intersection of src and tgt face area in 2D.

SourceFiles
    faceAreaWeightAMI2D.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMI2D_H
#define faceAreaWeightAMI2D_H

#include "advancingFrontAMI.H"
#include "AABBTree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faceAreaWeightAMI2D Declaration
\*---------------------------------------------------------------------------*/

class faceAreaWeightAMI2D
:
    public advancingFrontAMI
{
protected:

    // Protected Data

        //- Face bounding box factor
        scalar Cbb_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const faceAreaWeightAMI2D&) = delete;

        //- Helper function to write non-matched source faces to the set
        //- of candidate faces
        void writeNoMatch
        (
            const label srcFacei,
            const labelList& tgtFaceCandidates,
            const boundBox& srcFaceBb
        ) const;


        // Evaluation

            //- Calculate and store the area of intersection between source and
            //- target faces
            void storeInterArea
            (
                const label srcFacei,
                const label tgtFacei,
                DynamicList<label>& srcAddr,
                DynamicList<scalar>& srcWght,
                DynamicList<vector>& srcCtr,
                DynamicList<label>& tgtAddr,
                DynamicList<scalar>& tgtWght
            ) const;


            //- Return the set of tgt face IDs that overlap the src face bb
            labelList overlappingTgtFaces
            (
                const AABBTree<face>& tree,
                const List<boundBox>& tgtFaceBbs,
                const boundBox& srcFaceBb
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMI2D");


    // Constructors

        //- Construct from dictionary
        faceAreaWeightAMI2D
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        faceAreaWeightAMI2D
        (
            const bool requireMatch,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh,
            const bool restartUncoveredSourceFace = true
        );

        //- Construct as copy
        faceAreaWeightAMI2D(const faceAreaWeightAMI2D& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>(new faceAreaWeightAMI2D(*this));
        }


    //- Destructor
    virtual ~faceAreaWeightAMI2D() = default;


    // Member Functions

        //- Update addressing, weights and (optional) centroids
        virtual bool calculate
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        );

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
