/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::FacePostProcessing

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Records particle face quantities on used-specified face zone

    Currently supports:
        accumulated mass
        average mass flux

SourceFiles
    FacePostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef FacePostProcessing_H
#define FacePostProcessing_H

#include "CloudFunctionObject.H"
#include "faceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FacePostProcessing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FacePostProcessing
:
    public CloudFunctionObject<CloudType>
{
        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


    // Private Data

        //- Face zone IDs
        labelList faceZoneIDs_;

        //- Surface output format
        const word surfaceFormat_;

        //- Should be reset/cleared on writing?
        bool resetOnWrite_;

        //- Flag to indicate whether data should be written to file
        bool logToFile_;

        //- Total time
        scalar totalTime_;

        //- Mass storage
        List<scalarField> mass_;

        //- Mass total storage
        List<scalarField> massTotal_;

        //- Mass flow rate storage
        List<scalarField> massFlowRate_;

        //- Output file pointer per zone
        PtrList<OFstream> outputFilePtr_;

        //- Last calculation time
        scalar timeOld_;


    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile
        (
            const word& zoneName,
            const label zoneI,
            const label nFaces,
            const scalar totArea
        );


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("facePostProcessing");


    // Constructors

        //- Construct from dictionary
        FacePostProcessing
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        FacePostProcessing(const FacePostProcessing<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new FacePostProcessing<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FacePostProcessing() = default;


    // Member Functions

        //- Should data be reset on write?
        inline bool resetOnWrite() const;

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

#include "FacePostProcessingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FacePostProcessing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
