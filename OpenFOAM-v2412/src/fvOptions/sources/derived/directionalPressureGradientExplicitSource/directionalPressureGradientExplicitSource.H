/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::fv::directionalPressureGradientExplicitSource

Group
    grpFvOptionsSources

Description
    Applies an explicit pressure gradient source in such a way to deflect the
    flow towards an specific direction (\c flowDir). Alternatively add an extra
    pressure drop in the \c flowDir direction using a model.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    directionalPressureGradientExplicitSource1
    {
        // Mandatory entries (unmodifiable)
        type            directionalPressureGradientExplicitSource;

        // Mandatory entries (unmodifiable)
        model       <modelName>;
        fields      (<fieldName>);

        // Mandatory entries (runtime modifiable)
        flowDir     (1 1 0);
        faceZone    <faceZoneName>;

        // Conditional mandatory entries (unmodifiable)

            // when <timePath>/uniform/<name>Properties file exists
            gradient    <vectorField>; // reading from the aforementioned file

            // when model=DarcyForchheimer
            // deltaP = (D*mu + 0.5*rho*magUn)*magUn*length
            D           5e7;
            I           0;
            length      1e-3;

            // when model=constant
            pressureDrop    40;

            // when model=volumetricFlowRateTable
            outOfBounds     clamp;
            fileName        "volFlowRateTable";

        // Optional entries (runtime modifiable)
        relaxationFactor    0.3;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: directionalPressureGradientExplicitSource <!--
                -->                                 | word   | yes   | -
      model     | Pressure drop model [Pa]          | word   | yes   | -
      fields    | Name of operand field             | word   | yes   | -
      gradient  | Initial pressure gradient field   | vectorField | cndtnl | -
      flowDir   | Deflection flow direction         | vector | yes   | -
      faceZone  | Name of upstream faceZone         | word   | yes   | -
      relaxationFactor | Relaxation factor for flow deflection | scalar | no | 0.3
      D         | Darcy pressure loss coefficient   | scalar | cndtnl | -
      I         | Inertia pressure lost coefficient | scalar | cndtnl | -
      length    | Porous media length               | scalar | cndtnl | -
      presssureDrop | Constant pressure drop        | scalar | cndtnl | -
      fileName  | Interpolation table for volumetric flow rate <!--
                -->                     | interpolationTable | cndtnl | -
    \endtable

    Options for the \c model entry:
    \verbatim
      volumetricFlowRateTable | Pressure-gradient file
      constant                | Constant pressure drop
      DarcyForchheimer        | Darcy-Forchheimer model
    \endverbatim

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note
  - In order to obtain the upwind velocities this function loops over
    the slaves cells of the faceZone specified in the dictionary, on the other
    hand, the cellZone to which this source term is applied should be composed
    of the master cells and they should be 'downwind' the faceZone.

SourceFiles
    directionalPressureGradientExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef directionalPressureGradientExplicitSource_H
#define directionalPressureGradientExplicitSource_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvOption.H"
#include "cellSetOption.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
          Class directionalPressureGradientExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class directionalPressureGradientExplicitSource
:
    public fv::cellSetOption
{
public:

    //  Public Enumeration

        //- Modes of pressure drop
        enum pressureDropModel
        {
            pVolumetricFlowRateTable,
            pConstant,
            pDarcyForchheimer
        };


private:

    // Private Data

        static const Enum<pressureDropModel> pressureDropModelNames_;

        //- Pressure drop model
        pressureDropModel model_;

        //- Pressure gradient before correction
        vectorField gradP0_;

        //- Change in pressure gradient
        vectorField dGradP_;

        //- Pressure drop due to porous media
        vectorField gradPporous_;

        //- Flow direction
        vector flowDir_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> invAPtr_;

        //- Darcy pressure loss coefficient
        scalar D_;

        //- Inertia pressure lost coefficient
        scalar I_;

        //- Porous media length
        scalar length_;

        //- Constant pressure drop
        scalar pressureDrop_;

        //- Volumetric flow rate vs pressure drop table
        interpolationTable<scalar> flowRate_;

        //- Name of the faceZone at the heat exchange inlet
        word faceZoneName_;

        //- Id for the face zone
        label zoneID_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Cells faces mapping
        labelList cellFaceMap_;


    // Private Member Functions

        //- Init
        void initialise();

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const vectorField& gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("directionalPressureGradientExplicitSource");


    // Constructors

        //- Construct from explicit source name and mesh
        directionalPressureGradientExplicitSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        directionalPressureGradientExplicitSource
        (
            const directionalPressureGradientExplicitSource&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const directionalPressureGradientExplicitSource&
        ) = delete;


    //- Destructor
    virtual ~directionalPressureGradientExplicitSource() = default;


    // Member Functions

        // Evaluation

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream& os) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
