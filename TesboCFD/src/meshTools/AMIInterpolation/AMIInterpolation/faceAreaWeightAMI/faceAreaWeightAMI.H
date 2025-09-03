/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::faceAreaWeightAMI

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method

    Searching is performed using an advancing front.

SourceFiles
    faceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMI_H
#define faceAreaWeightAMI_H

#include "advancingFrontAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

class faceAreaWeightAMI
:
    public advancingFrontAMI
{
private:

    // Private data

        //- Flag to restart uncovered source faces
        const bool restartUncoveredSourceFace_;


protected:

    // Protected Member Functions

        //- No copy assignment
        void operator=(const faceAreaWeightAMI&) = delete;


        // Marching front

            //- Calculate addressing, weights and centroids using temporary
            //- storage
            virtual void calcAddressing
            (
                List<DynamicList<label>>& srcAddress,
                List<DynamicList<scalar>>& srcWeights,
                List<DynamicList<point>>& srcCentroids,
                List<DynamicList<label>>& tgtAddress,
                List<DynamicList<scalar>>& tgtWeights,
                label srcFacei,
                label tgtFacei
            );

            //- Determine overlap contributions for source face srcFacei
            virtual bool processSourceFace
            (
                const label srcFacei,
                const label tgtStartFacei,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

            //- Attempt to re-evaluate source faces that have not been included
            virtual void restartUncoveredSourceFace
            (
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

            //- Set the source and target seed faces
            virtual bool setNextFaces
            (
                label& startSeedi,
                label& srcFacei,
                label& tgtFacei,
                const bitSet& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                const bool errorOnNotFound = true
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            virtual void calcInterArea
            (
                const label srcFacei,
                const label tgtFacei,
                scalar& area,
                vector& centroid
            ) const;

            //- Return true if faces overlap
            virtual bool overlaps
            (
                const label srcFacei,
                const label tgtFacei,
                const scalar threshold
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMI");


    // Constructors

        //- Construct from dictionary
        faceAreaWeightAMI
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        faceAreaWeightAMI
        (
            const bool requireMatch,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh,
            const bool restartUncoveredSourceFace = true
        );

        //- Construct as copy
        faceAreaWeightAMI(const faceAreaWeightAMI& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>(new faceAreaWeightAMI(*this));
        }


    //- Destructor
    virtual ~faceAreaWeightAMI() = default;


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
