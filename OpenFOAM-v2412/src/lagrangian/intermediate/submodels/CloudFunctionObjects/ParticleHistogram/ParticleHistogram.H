/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::ParticleHistogram

Description
    Computes a histogram for the distribution of particle diameters
    and corresponding number of particles hitting on a given list of
    patches or face zones.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | \<case\>/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

    The output file contains four columns, from the first to the last:
    \vartable
        dEdge1      | Left bin-edge of particle diameter
        dEdge2      | Right bin-edge of particle diameter
        nParticles  | Number of particles falling into the corresponding bin
        nParticlesCumulative | Cumulative number of binned particles
    \endvartable

Usage
    Minimal example by using
    \c constant/reactingCloud1Properties.cloudFunctions:
    \verbatim
    particleHistogram1
    {
        // Mandatory entries
        type                 particleHistogram;
        nBins                <label>;
        min                  <scalar>;
        max                  <scalar>;
        maxStoredParcels     <scalar>;

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
      type         | Type name: particleHistogram       | word   |  yes  | -
      nBins        | Number of histogram bins           | label  |  yes  | -
      max          | Maximum value of histogram data    | scalar |  yes  | -
      min          | Minimum value of histogram data    | scalar |  yes  | -
      maxStoredParcels | Maximum number of parcels to process | label | yes | -
      patches      | Names of operand patches          | wordRes |  choice  | -
      faceZones    | Names of operand face zones       | wordRes |  choice  | -
    \endtable

    The inherited entries are elaborated in:
      - \link CloudFunctionObject.H \endlink
      - \link writeFile.H \endlink

Note
  - The underlying type of \c maxStoredParcels is set as a scalar for I/O.

SourceFiles
    ParticleHistogram.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ParticleHistogram_H
#define Foam_ParticleHistogram_H

#include "CloudFunctionObject.H"
#include "writeFile.H"
#include "MinMax.H"
#include "cloudFunctionObjectTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ParticleHistogram Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleHistogram
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

        //- Number of data bins
        const label nBins_;

        //- Maximum number of parcels to store - set as a scalar for I/O
        const scalar maxStoredParcels_;

        //- Min/Max value of histogram data
        const scalarMinMax range_;

        //- Bin edges of histogram
        scalarField binEdges_;

        //- Accumulated number of particles per surface
        //- binned according to the histogram settings
        List<scalarList> nParticlesCumulative_;

        // List of surface-hit particle diameters
        List<DynamicList<scalar>> dParticles_;

        // List of number of surface-hit particles
        List<DynamicList<scalar>> nParticles_;


    // Private Member Functions

        //- Write output file header
        void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("particleHistogram");


    // Constructors

        //- Construct from dictionary
        ParticleHistogram
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        ParticleHistogram(const ParticleHistogram<CloudType>& ph);

        //- No copy assignment
        void operator=(const ParticleHistogram<CloudType>&) = delete;

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleHistogram<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleHistogram() = default;


    // Member Functions

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
    #include "ParticleHistogram.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
