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
    Foam::fv::tabulatedNTUHeatTransfer

Group
    grpFvOptionsSources

Description
    Applies a tabulated heat transfer model for inter-region heat exchanges.
    The heat flux is calculated based on the Number of Thermal Units (NTU).

    A 2-D table of NTU as functions of the primary
    and secondary mass flow rates is required.

    The exchanger geometry can be specified using either:
    - \c calculated:
        - inlet area of each region and core volume determined by interrogating
          mesh patches, and mesh-to-mesh interpolation volume
        - calculated core volume can be partially blocked by specifying a
          \c coreBlockageRatio [0-1] entry
    - \c user:
        - inlet area of each region provided by the user
        - core volume automatically calculated by the mesh-to-mesh
          interpolation volume if not provided by user

    Heat transfer coefficient calculated by:

    \f[
        htc = C_{min} \frac{NTU}{V_{core}}
    \f]

    Where \f$ C_{min} \f$ is given by:

    \f[
        C_{min} = min \left(Cp_1 \dot{m}_1, Cp_2 \dot{m}_2 \right)
    \f]

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    tabulatedNTUHeatTransfer1
    {
        // Mandatory entries (unmodifiable)
        type            tabulatedNTUHeatTransfer;

        // Mandatory entries (runtime modifiable)
        geometryMode    <geometryMode>;
        outOfBounds     clamp;
        file            "ntuTable";

        // Optional entries (runtime modifiable)
        U               <Uname>;
        Unbr            <UnbrName>;
        rho             <rhoName>;
        rhoNbr          <rhoNbrName>;

        // Conditional mandatory entries (runtime modifiable)

            // when geometryMode=user
            Ain             0.01728;
            AinNbr          0.3456;

            // when geometryMode=calculated
            inletPatch      inlet_HWK;
            inletPatchNbr   inlet_air;
            inletBlockageRatio 0.10;
            inletBlockageRatioNbr 0.04;
            coreBlockageRatio 0;

        // Conditional optional entries (runtime modifiable)

            // when geometryMode=user
            Vcore           0.01244;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type  | Reqd | Dflt
      type      | Type name: tabulatedNTUHeatTransfer | word  | yes  | -
      geometryMode | Geometry mode type               | word  | yes  | -
      file  | Heat transfer coefficient table  | interpolate2DTable | yes | -
      U         | Name of operand velocity field   | word  | no   | U
      Unbr      | Name of operand neighbour velocity field | word | no | U
      Ain       | Inlet area [m2]                     | scalar | cndtnl | -
      AinNbr    | Neighbour region inlet area [m2]    | scalar | cndtnl | -
      Vcore     | Heat exchanger core volume          | scalar | cndtnl | -
      inletPatch | Name of inlet patch                | word   | cndtnl | -
      inletNbrPatch | Name of inlet patch neighbour   | word   | cndtnl | -
      inletBlockageRatio | Inlet patch blockage ratio [0, 1]   <!--
                    -->                               | scalar | cndtnl | -
      inletBlockageRatioNbr | Inlet neighbour patch blockage ratio [0, 1] <!--
                    -->                               | scalar | cndtnl | -
      coreBlockageRatio  | Core volume blockage ratio [0, 1] <!--
                    -->                               | scalar | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
     - \link interRegionHeatTransfer.H \endlink
     - \link interpolate2DTable.H \endlink

    Options for the \c geometryMode entry:
    \verbatim
      calculated    | Use settings computed internally
      user          | Use settings provided by the user
    \endverbatim

    Example usage:
    \verbatim
    coolerToAir
    {
        type            tabulatedNTUHeatTransfer;
        active          yes;

        tabulatedNTUHeatTransferCoeffs
        {
            interpolationMethod cellVolumeWeight;
            nbrRegion       air;
            master          true;

            fields          (h);
            outOfBounds     clamp;
            file            "ntuTable";
            nbrModel        airToCooler;
            semiImplicit    no;


            geometryMode    user;
            Ain             0.01728;
            AinNbr          0.3456;
            Vcore           0.01244;  // Optional

            // geometryMode    calculated;
            // inletPatch      inlet_HWK;
            // inletPatchNbr   inlet_air;
            // inletBlockageRatio 0.10;
            // inletBlockageRatioNbr 0.04;
            // coreBlockageRatio 0;
        }
    }
    \endverbatim

See also
  - Foam::fv::interRegionHeatTransferModel
  - Foam::fv::constantHeatTransfer
  - Foam::fv::tabulatedHeatTransfer
  - Foam::fv::variableHeatTransfer

SourceFiles
    tabulatedNTUHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedNTUHeatTransfer_H
#define tabulatedNTUHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"
#include "Enum.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedNTUHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedNTUHeatTransfer
:
    public interRegionHeatTransferModel
{
public:

    // Public Enumerations

        //- Options for the geometry mode type
        enum geometryModeType
        {
            gmCalculated,
            gmUser
        };

        //- Names for geometryModeType
        static const Enum<geometryModeType> geometryModelNames_;


private:

    // Private Data

        //- Name of operand velocity field
        word UName_;

        //- Name of operand neighbour velocity field
        word UNbrName_;

        //- Name of operand density field
        word rhoName_;

        //- Name of operand neighbour density field
        word rhoNbrName_;

        //- Pointer to 2-D look-up table of NTU f(mDot1, mDot2)
        autoPtr<interpolation2DTable<scalar>> ntuTable_;

        //- Geometry input mode
        geometryModeType geometryMode_;

        //- Inlet area [m2]
        scalar Ain_;

        //- Neighbour region inlet area [m2]
        scalar AinNbr_;

        //- Heat exchanger core volume
        scalar Vcore_;


    // Private Member Functions

        //- NTU table helper
        const interpolation2DTable<Foam::scalar>& ntuTable();

        //- Thermophysical properties helper
        const basicThermo& thermo(const fvMesh& mesh) const;

        //- Initialise geometry
        void initialiseGeometry();


public:

    //- Runtime type information
    TypeName("tabulatedNTUHeatTransfer");


    // Constructors

        //- Construct from components
        tabulatedNTUHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        tabulatedNTUHeatTransfer(const tabulatedNTUHeatTransfer&) = delete;

        //- No copy assignment
        void operator=(const tabulatedNTUHeatTransfer&) = delete;


    //- Destructor
    virtual ~tabulatedNTUHeatTransfer() = default;


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc();

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
