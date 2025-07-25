/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::ParticlePostProcessing

Description
    Writes out various standard Lagrangian data elements of
    particles hitting on a given list of patches or face zones.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | \<case\>/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

    The output file contains columns depending on the \c fields entry.

Usage
    Minimal example by using
    \c constant/reactingCloud1Properties.cloudFunctions:
    \verbatim
    ParticlePostProcessing1
    {
        // Mandatory entries
        type                 particlePostProcessing;
        maxStoredParcels     <scalar>;

        // Optional entries
        fields               (<wordRes>);

        // Conditional entries

            // Option-1
            patches          (<wordRes>);

            // Option-2
            faceZones        (<wordRes>);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type   | Reqd  | Deflt
      type         | Type name: particlePostProcessing  | word   |  yes  | -
      maxStoredParcels | Maximum number of parcels to process | label | yes | -
      fields       | Names of standard Lagrangian fields | wordRes | no  | all
      patches      | Names of operand patches          | wordRes |  choice  | -
      faceZones    | Names of operand face zones       | wordRes |  choice  | -
    \endtable

    The inherited entries are elaborated in:
      - \link CloudFunctionObject.H \endlink
      - \link writeFile.H \endlink

Note
  - The underlying type of \c maxStoredParcels is set as a scalar for I/O.
  - Empty files indicate that no particles hit selected surfaces.

SourceFiles
    ParticlePostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ParticlePostProcessing_H
#define Foam_ParticlePostProcessing_H

#include "CloudFunctionObject.H"
#include "writeFile.H"
#include "cloudFunctionObjectTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class ParticlePostProcessing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticlePostProcessing
:
    public CloudFunctionObject<CloudType>,
    public functionObjects::writeFile
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::particleType parcelType;

        //- Collector surfaces
        cloudFunctionObjectTools::collector collector_;

        //- Maximum number of parcels to store - set as a scalar for I/O
        scalar maxStoredParcels_;

        //- Field header
        string header_;

        //- Field name filters
        wordRes fields_;

        //- List of times for each data record
        List<DynamicList<scalar>> times_;

        //- List of output data per surface
        List<DynamicList<string>> data_;


    // Private Member Functions

        //- Write output file header
        void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("particlePostProcessing");


    // Constructors

        //- Construct from dictionary
        ParticlePostProcessing
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        ParticlePostProcessing(const ParticlePostProcessing<CloudType>& ppp);

        //- No copy assignment
        void operator=(const ParticlePostProcessing<CloudType>&) = delete;

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticlePostProcessing<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticlePostProcessing() = default;


    // Member Functions

    // Access

        //- Return maximum number of parcels to store per surface
        label maxStoredParcels() const noexcept { return maxStoredParcels_; }


    // Evaluation

        //- Post-patch hook
        virtual bool postPatch
        (
            const parcelType& p,
            const polyPatch& pp,
            const typename parcelType::trackingData& td
        );

        //- Post-face hook
        virtual bool postFace
        (
            const parcelType& p,
             const typename parcelType::trackingData& td
        );


    // I-O

        //- Write post-processing info
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticlePostProcessing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
