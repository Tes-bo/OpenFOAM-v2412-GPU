/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd
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
    Foam::compressible::thermalBaffleFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a coupled temperature condition between
    multiple mesh regions.

    The regions are generally referred to as the:
      - primary region,
      - and baffle region.

    The primary region creates the baffle region and evolves its energy
    equation either:
      - 1-D, normal to each patch face
      - 3-D, normal and tangential components

    The thermodynamic properties of the baffle material are specified via
    dictionary entries on the master patch.

    The underlying polyPatch on the primary region needs to be defined
    as a complete mappedWall. i.e:

    \verbatim
        type            mappedWall;
        sampleMode      nearestPatchFace;
        samplePatch     <patchName>;
        sampleRegion    <regionName>;

        coupleGroup     <coupleGroupName>;

    \endverbatim

    The /keyword{internal} is true when the thermal baffle is used on internal
    faces and the baffle is coupled on both sides (bottom and top).

    If this baffle is used on a patch set internal = false. By default is
    true.

    \table
   Property         | Description             | Required    | Default value
    Tnbr            | name of the field       | no | T
    qrNbr           | name of the radiative flux in the nbr region | no | none
    qr              | name of the radiative flux in this region | no | none
    thicknessLayers | inherited from temperatureCoupledBase | inherited
    kappaLayers     | inherited from temperatureCoupledBase | inherited
    kappaMethod     | inherited from temperatureCoupledBase | inherited
    kappa           | inherited from temperatureCoupledBase | inherited
    thermalInertia  | Add thermal inertia to wall node      | no | false
    region          | Name of the solid region              | yes
    active          | Active                                | yes
    internal        | Is the baffle an internal regiion     | no  | true
    thermoType      | Solid thermo type                     | yes
    mixture         | Mixture solid                         | yes
    radiation       | Radiation model in the solid          | yes
    extrudeModel    | Extrude model                         | yes
    nLayers         | Number of cell across the ragion      | yes
    expansionRatio  | Expansion ratio for the solid region  | yes
    columnCells     | true for 1D or false for 3D region    | yes

    \endtable

Usage
    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        type                compressible::thermalBaffle;

        // Underlying coupled boundary condition
        Tnbr               T;
        kappaMethod        fluidThermo; // or solidThermo
        kappa              none;
        qrNbr              qr;//or none.Name of qr field on neighbourregion
        qr                 none;// or none.Name of qr field on localregion
        value              uniform 300;

        // Baffle region name
        region              baffleRegion;
        active              yes;
        internal            yes;

        // Solid thermo in solid region
        thermoType
        {
            type            heSolidThermo;
            mixture         pureMixture;
            transport       constIso;
            thermo          hConst;
            equationOfState rhoConst;
            specie          specie;
            energy          sensibleEnthalpy;
        }

        mixture
        {
            specie
            {
                molWeight       20;
            }
            transport
            {
                kappa           0.01;
            }
            thermodynamics
            {
                Hf              0;
                Cp              15;
            }
            density
            {
                rho             80;
            }
        }

        radiation
        {
            radiationModel  opaqueSolid;
            absorptionEmissionModel none;
            scatterModel    none;
        }

        // Extrude model for new region

        extrudeModel        linearNormal;
        nLayers             50;
        expansionRatio      1;
        columnCells         false; //3D or 1D
        linearNormalCoeffs
        {
            thickness           0.02;
        }

    }
    \endverbatim

    If the thermal baffle is used on internal faces of the primary region,
    then the slave in this region is defined as :

    Slave patch on primary region:
    \verbatim
    <slavePatchName>
    {
        type                compressible::thermalBaffle;
        kappaMethod         fluidThermo;
        kappa               none;
        value               uniform 300;
    \endverbatim

    Patches on baffle region:
    \verbatim
    bottom
    {
        type                compressible::thermalBaffle;
        kappaMethod         solidThermo;
        kappa               none;
        value               uniform 300;
    }

    top
    {
        type                compressible::thermalBaffle;
        kappaMethod         solidThermo;
        kappa               none;
        value               uniform 300;
    }
    \endverbatim

See also
    Foam::compressible::turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    Foam::regionModels::thermalBaffleModels::thermalBaffleModel

SourceFiles
    thermalBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffleFvPatchScalarField_H
#define thermalBaffleFvPatchScalarField_H


#include "autoPtr.H"
#include "regionModel.H"
#include "thermalBaffleModel.H"
#include "extrudePatchMesh.H"
#include "turbulentTemperatureRadCoupledMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class thermalBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffleFvPatchScalarField
:
    public turbulentTemperatureRadCoupledMixedFvPatchScalarField
{
    // Typedefs

        //- The baffle region model
        typedef
            regionModels::thermalBaffleModels::thermalBaffleModel
            baffleType;


    // Private Data

        //- Enumeration of patch IDs
        enum patchID
        {
            bottomPatchID,
            topPatchID,
            sidePatchID
        };

        //- Is the baffle owner
        bool owner_;

        //- Is the baffle internal
        bool internal_;

        //- Thermal baffle
        autoPtr<baffleType> baffle_;

        //- Dictionary
        dictionary dict_;

        //- Auto pointer to extrapolated mesh from patch
        autoPtr<extrudePatchMesh> extrudeMeshPtr_;


    // Private Member Functions

        //- Extrude mesh
        void createPatchMesh();

public:

    //- Runtime type information
    TypeName("compressible::thermalBaffle");


    // Constructors

        //- Construct from patch and internal field
        thermalBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalBaffleFvPatchScalarField onto a new patch
        thermalBaffleFvPatchScalarField
        (
            const thermalBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        thermalBaffleFvPatchScalarField
        (
            const thermalBaffleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
