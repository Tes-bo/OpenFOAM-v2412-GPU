/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::heatTransferCoeffModels::faceZoneReferenceTemperature

Description
    Heat transfer coefficient calculation that employs the area-average
    temperature of a specified face zone as the reference temperature.

    The heat transfer coefficient is calculated by:

    \f[
        h = \frac{q}{T_{ref} - T_p}
    \f]

    where
    \vartable
      h       | Heat transfer coefficient                         [W/m^2/K]
      q       | Heat flux                                         [W/m^2]
      T_{ref} | Area average of reference face zone temperature   [K]
      T_p     | Temperature field of current patch                [K]
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Inherited entries
        ...

        // Mandatory entries
        htcModel            faceZoneReferenceTemperature;
        referenceFaceZone   <word>;

        // Optional entries
        referenceRegion     <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      htcModel  | Model name: faceZoneReferenceTemperature  | word | yes  | -
      referenceFaceZone | Name of reference face zone    | word | yes  | -
      referenceRegion | Name of region that reference face zone resides  <!--
                -->                                      | word | no   | region0
    \endtable

SourceFiles
    faceZoneReferenceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_faceZoneReferenceTemperature_H
#define heatTransferCoeffModels_faceZoneReferenceTemperature_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                Class faceZoneReferenceTemperature Declaration
\*---------------------------------------------------------------------------*/

class faceZoneReferenceTemperature
:
    public heatTransferCoeffModel
{
    // Private Data

        //- Patch index of reference face zone
        label faceZonei_;

        //- Name of region that the reference face zone resides
        word refRegionName_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;


    // Private Member Functions

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces(const dictionary& dict);

        //- Return area average of face-zone temperature field
        scalar faceZoneAverageTemperature();

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );


public:

    //- Runtime type information
    TypeName("faceZoneReferenceTemperature");


    // Constructors

        //- Construct from components
        faceZoneReferenceTemperature
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        faceZoneReferenceTemperature(const faceZoneReferenceTemperature&) = delete;

        //- No copy assignment
        void operator=(const faceZoneReferenceTemperature&) = delete;


    //- Destructor
    virtual ~faceZoneReferenceTemperature() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
