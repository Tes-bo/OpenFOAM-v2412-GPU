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
    Foam::FaceInteraction

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Face zone-based particle interactions.

    \verbatim
    faceInteraction1
    {
        type            faceInteraction;

        // List of (faceZone interactionType)
        faceZones
        (
            (faceZone1      stick)
            (faceZone2      escape)
            (faceZone3      rebound)
        );

        // Optional limiting to diameter range
        dMin            0;
        dMax            1;

        writeToFile     yes; // default = yes
    }
    \endverbatim

    File written per faceZone as:

      postProcessing/lagrangian/cloudName/modelName/time/modelName_faceZone.dat

SourceFiles
    FaceInteraction.C

\*---------------------------------------------------------------------------*/

#ifndef FaceInteraction_H
#define FaceInteraction_H

#include "CloudFunctionObject.H"
#include "Enum.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FaceInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FaceInteraction
:
    public CloudFunctionObject<CloudType>,
    public functionObjects::writeFile
{
public:

    // Public Data Types

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Enumeration defining the interaction types
        enum class interactionType
        {
            ESCAPE,  //!< particles escape/are removed
            STICK,   //!< particles stick to the faceZone faces
            REBOUND  //!< particles rebound from the faceZone faces
        };

        //- Names for the interaction types
        static const Enum<interactionType> interactionTypeNames_;


private:

    // Private Data

        //- Face zone IDs
        labelList faceZoneIDs_;

        //- Face zone bounding boxes
        List<boundBox> faceZoneBBs_;

        //- Face zone particle interaction type
        List<interactionType> faceZoneInteraction_;

        //- File per zone
        PtrList<OFstream> filePtrs_;

        //- Number of escape particles per zone
        labelList nEscapeParticles_;

        //- Number of stick particles per zone
        labelList nStickParticles_;

        //- Number of rebound particles per zone
        labelList nReboundParticles_;

        //- Minimum diameter threshold
        scalar dMin_;

        //- Maximum diameter threshold
        scalar dMax_;


protected:

    // Protected Member Functions

        //- Return true if this particle will be assessed
        bool processParticle(const parcelType& p, const label localZonei);

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("faceInteraction");


    // Constructors

        //- Construct from dictionary
        FaceInteraction
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        FaceInteraction(const FaceInteraction<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new FaceInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FaceInteraction() = default;


    // Member Functions

        //- Post-face hook
        virtual bool postFace
        (
            const parcelType& p,
            const typename parcelType::trackingData& td
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FaceInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
