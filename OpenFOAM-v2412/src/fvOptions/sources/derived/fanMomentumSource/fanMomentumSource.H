/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Louis Vittoz, SimScale GmbH
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::fv::fanMomentumSource

Group
    grpFvOptionsSources

Description
    This source term models the action of a fan on the flow. It calculates
    flow rate through a set of upstream faces of encompassing the cell zone. The
    set of upstream faces is automatically calculated based on the flow
    direction and the surrounding face zone. Based on the flow rate, the
    pressure gradient is calculated based on the fan pressure curve and the
    thickness of the fan (measurement section). The pressure gradient is then
    added to the momentum equation.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <name>
    {
        // Mandatory entries
        type                fanMomentumSource;
        fanCurve            <Function1<scalar>>;
        flowDir             <vector>;
        thickness           <scalar>;
        cellZone            <word>;
        faceZone            <word>;

        // Optional entries
        gradient            <scalar>;
        rho                 <scalar>;
        U                   <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Deflt
      type      | Type name: fanMomentumSource      | word  | yes  | -
      fanCurve  | Pressure drop vs flow-rate  | Function1\<scalar\>| yes | -
      flowDir   | Direction of the flow through the fan | vector | yes | -
      cellZone  | Cell zone representing the fan    | word   | yes | -
      faceZone  | Face zone around the cell zone    | word   | yes | -
      thickness | Thickness of the fan [m]          | scalar | yes | -
      gradient  | Initial pressure gradient         | scalar | no  | -
      rho       | Reference density for incompressible flow | scalar | no  | -
      U         | Name of velocity field            | word   | no  | U
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link Function1.H \endlink

Note
  - The fan curve needs to be provded for a pressure drop in [Pa] and
    is specified as a function of the volumetric flow rate in [m^3/s].

SourceFiles
    fanMomentumSource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_fanMomentumSource_H
#define Foam_fv_fanMomentumSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class fanMomentumSource Declaration
\*---------------------------------------------------------------------------*/

class fanMomentumSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Local list of pairs of upstream patch IDs and upstream face IDs
        //  For internal faces, the upstream patch ID is -1
        List<labelPair> upstreamPatchFaceInfo_;

        //- Cells in zone
        labelHashSet cellsInZones_;

        //- Volumetric flow rate [m3] vs. pressure drop [Pa] table
        autoPtr<Function1<scalar>> fanCurvePtr_;

        //- Name of the velocity field this function object operates on
        word UName_;

        //- Direction of flow through the fan
        vector flowDir_;

        //- Thickness of the fan, in [m]
        scalar thickness_;

        //- Pressure drop
        scalar gradPFan_;

        //- Id for the surrounding face zone
        label surroundingFaceZoneID_;

        //- Reference density for incompressible cases
        scalar rho_;

        //- Whether to use reference density
        bool useRefRho_;


    // Private Member Functions

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP, const scalar flowRate) const;

        //- Collect all faces upstream of
        //- the centre of gravity of the cell zone
        void initUpstreamFaces();

        //- Return the volumetric flow rate given flux and face-interpolated
        //  density. For incompressible cases, use surfaceScalarField::null()
        //  for rhof
        scalar calcFlowRate
        (
            const surfaceScalarField& phi,
            const surfaceScalarField& rhof
        ) const;


public:

    //- Runtime type information
    TypeName("fanMomentumSource");


    // Constructors

        //- Construct from explicit source name and mesh
        fanMomentumSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        fanMomentumSource(const fanMomentumSource&) = delete;

        //- No copy assignment
        void operator=(const fanMomentumSource&) = delete;


    //- Destructor
    virtual ~fanMomentumSource() = default;


    // Member Functions

    // Evaluation

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


    // IO

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
