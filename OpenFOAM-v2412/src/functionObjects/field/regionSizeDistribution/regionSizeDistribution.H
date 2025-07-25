/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::functionObjects::regionSizeDistribution

Group
    grpFieldFunctionObjects

Description
    Creates a droplet size distribution via interrogating a continuous phase
    fraction field.

    Looks up a phase-fraction (alpha) field and splits the mesh into regions
    based on where the field is below the threshold value.  These
    regions ("droplets") can now be analysed.

    Regions:
    - print the regions connected to a user-defined set of patches.
      (in spray calculation these form the liquid core)
    - print the regions with too large volume.  These are the 'background'
      regions.
    - (debug) write regions as a volScalarField
    - (debug) print for all regions the sum of volume and alpha*volume

    Output (volume scalar) fields include:
    - alpha_liquidCore : alpha with outside liquid core set to 0
    - alpha_background : alpha with outside background set to 0.

    Histogram:
    - determine histogram of diameter (given minDiameter, maxDiameter, nBins)
    - write graph of number of droplets per bin
    - write graph of sum, average and deviation of droplet volume per bin
    - write graph of sum, average and deviation of user-defined fields.  For
      volVectorFields these are those of the 3 components and the magnitude.
    - (optional) write graph of histogram of centroids on iso planes
      downstream of the injector determined by origin, direction and maxDiameter
      up to maxDownstream

    Operands:
    \table
      Operand           | Type              | Location
      input             | -                 | -
      output file       | -                 <!--
                    --> | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<files\>
      output field      | -                 | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    regionSizeDistribution1
    {
        // Mandatory entries (unmodifiable)
        type            regionSizeDistribution;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        field           alpha;
        patches         (inlet);
        fields          (p U);
        threshold       0.4;
        maxDiameter     5e-5;
        nBins           100;
        setFormat       gnuplot;

        // Optional entries (runtime modifiable)
        minDiameter     0.0;
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type    axes;
                e3      (0 0 1);
                e1      (1 0 0);
            }
        }

        // Optional downstream iso-plane bins
        isoPlanes       true;

            // Mandatory entries if isoPlanes=true   (runtime modifiable)
            // Plane normal and point definition
            origin          (1e-4 0 5e-4);
            direction       (1 0 1);

            // Maximum diameter of the cylinder formed by the origin point
            // and direction
            maxD            3e-4;

            // Maximum downstream distance
            maxDownstream   6e-4;

            // Number of iso-plane bins
            nDownstreamBins 20;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: regionSizeDistribution  | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Phase field to interrogate         | word |  yes  | -
      patches      | Patches wherefrom the liquid core is identified <!--
                                                    --> | wordList | yes | -
      fields     | Fields to sample                     | wordList | yes | -
      threshold  | Phase fraction applied to delimit regions | scalar | yes | -
      maxDiameter  | Maximum droplet diameter           | scalar | yes | -
      minDiameter  | Minimum droplet diameter           | scalar | no  | 0.0
      nBins        | Number of bins for histogram       | label  | yes | -
      setFormat    | Output format                      | word   | yes | -
      isoPlanes    | Flag for isoPlanes                 | bool | no    | false
      origin       | Origin of the plane when isoPlanes is used <!--
                                                    --> | vector | yes | -
      direction   <!--
      --> | Direction of the plane when isoPlanes is used | vector | yes | -
      maxD        | Maximum diameter of the sampling <!--
                --> cylinder when isoPlanes is used     | vector | yes | -
      nDownstreamBins  <!--
      --> | Number of bins when isoPlanes is used       | label | yes  | -
      maxDownstream  <!--
      --> | Maximum distance from origin when isoPlanes is used <!--
                                                    --> | scalar | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::regionSizeDistribution

SourceFiles
    regionSizeDistribution.C
    regionSizeDistributionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_regionSizeDistribution_H
#define Foam_functionObjects_regionSizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordSetWriter.H"
#include "Map.H"
#include "volFieldsFwd.H"
#include "wordRes.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class regionSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class regionSizeDistribution
:
    public functionObjects::fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Number of bins
        label nBins_;

        //- Name of field
        word alphaName_;

        //- Clip value
        scalar threshold_;

        //- Maximum droplet diameter
        scalar maxDiam_;

        //- Minimum droplet diameter
        scalar minDiam_;

        //- Patches to walk from
        wordRes patchNames_;

        //- Names of fields to sample on regions
        wordRes fields_;

        //- Output formatter to write
        mutable autoPtr<coordSetWriter> formatterPtr_;

        //- Optional coordinate system
        autoPtr<coordinateSystem> csysPtr_;

        // Optional extra definition of bins on planes downstream to the origin
        // point and maximum diameter

            //- Optional origin point
            vector origin_;

            //- Optional plane normal direction
            vector direction_;

            //- Optional maximum diameter on plane
            scalar maxDiameter_;

            //- Optional number of bins for
            scalar nDownstreamBins_;

            //- Optional maximum downstream coordinate from origin
            scalar maxDownstream_;

            //- Switch to enable iso-planes sampling
            bool isoPlanes_;


    // Private Member Functions

        //- Write volfields with the parts of alpha which are not
        //- droplets (liquidCore, backGround)
        void writeAlphaFields
        (
            const regionSplit& regions,
            const labelHashSet& keepRegions,
            const Map<scalar>& regionVolume,
            const volScalarField& alpha
        ) const;

        //- Mark all regions starting at patches
        labelHashSet findPatchRegions(const regionSplit&) const;

        //- Helper: divide if denom != 0
        static tmp<scalarField> divide(const scalarField&, const scalarField&);

        //- Given per-region data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const scalarField& sortedField,     // per region field data

            const labelList& indices,           // index of bin for each region
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;

        //- Given per-cell data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const scalarField& cellField,       // per cell field data
            const regionSplit& regions,         // per cell the region(=droplet)
            const labelList& sortedRegions,     // valid regions in sorted order
            const scalarField& sortedNormalisation,

            const labelList& indices,           // index of bin for each region
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;


public:

    //- Runtime type information
    TypeName("regionSizeDistribution");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        regionSizeDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        regionSizeDistribution(const regionSizeDistribution&) = delete;

        //- No copy assignment
        void operator=(const regionSizeDistribution&) = delete;


    // Destructor
    virtual ~regionSizeDistribution() = default;


    // Member Functions

        //- Read the regionSizeDistribution data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate the regionSizeDistribution and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
