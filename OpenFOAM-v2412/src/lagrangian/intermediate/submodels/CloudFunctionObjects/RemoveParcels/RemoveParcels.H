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
    Foam::RemoveParcels

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Removes parcels that hit user-specified face zone faces

    Example usage:
    \verbatim
    removeParcels1
    {
        type            removeParcels;
        log             yes;
        resetOnWrite    no;
        resetOnStart    no;
        faceZones       (cycLeft cycRight);
    }
    \endverbatim

SourceFiles
    RemoveParcels.C

\*---------------------------------------------------------------------------*/

#ifndef RemoveParcels_H
#define RemoveParcels_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                     Class RemoveParcels Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RemoveParcels
:
    public CloudFunctionObject<CloudType>
{
        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


    // Private Data

        //- Face zone IDs
        labelList faceZoneIDs_;

        //- Number of parcels removed per face zone
        List<label> nParcels_;

        //- Mass of parcels removed per face zone
        List<scalar> mass_;

        //- Type ID of parcels to remove; all selected if not set (<0)
        label typeId_;

        //- Flag to indicate whether data should be written to file
        bool logToFile_;

        //- Flag to reset counters on write
        bool resetOnWrite_;

        //- Flag to reset on start
        bool resetOnStart_;

        //- Output file pointer per zone
        PtrList<OFstream> outputFilePtr_;


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
    TypeName("removeParcels");


    // Constructors

        //- Construct from dictionary
        RemoveParcels
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        RemoveParcels(const RemoveParcels<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new RemoveParcels<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RemoveParcels() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve
        (
            const typename parcelType::trackingData& td
        );

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
    #include "RemoveParcels.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
