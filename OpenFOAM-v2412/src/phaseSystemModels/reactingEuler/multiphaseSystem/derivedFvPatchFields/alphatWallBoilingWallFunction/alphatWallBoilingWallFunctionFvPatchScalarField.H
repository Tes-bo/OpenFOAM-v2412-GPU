/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd
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
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField

Description
    A thermal wall function for simulation of boiling wall.

    This alpha wall function can handle the following regimes:
      - single phase
      - subcooled nucleate wall boiling
      - transitional boiling
      - film boiling

    The wall function uses a partition method to transfer heat either
    to the liquid or vapor phase. At the moment, this function works
    in a wall temperature fixed mode, i.e. there is no consideration
    for the sudden change of heat transfer coefficient (htc) after
    reaching TDBN (deviation from nucleate boiling temperature).

    References:
    \verbatim
        Srinivasan, V., Moon, K. M., Greif, D.,
        Wang, D. M., & Kim, M. H. (2010).
        Numerical simulation of immersion quenching
        process of an engine cylinder head.
        Applied Mathematical Modelling, 34(8), 2111-2128.
        DOI:10.1016/j.apm.2009.10.023
    \endverbatim


    For the single phase non-boiling regime the standard
    \c JayatillekeWallFunction is used.

    For the sub-cool nucleate boiling regime the following runtime
    selectable submodels are used:
      - nucleation site density
      - bubble departure frequency
      - bubble departure diameter

    Implements a version of the well-known RPI wall boiling model
    (Kurul & Podowski, 1991). The model implementation is similar to the model
    described by Peltola & Pättikangas (2012) but has been extended with the
    wall heat flux partitioning models.

    References:
    \verbatim
        Kurul, N., & Podowski, M. Z. (1991).
        On the modeling of multidimensional effects in boiling channels.
        Proceedings of the 27th National Heat Transfer Conference.
        Minneapolis, Minn, USA, July 28-31, 1991.
        ISBN: 0-89448-162-1, pp. 30-40

        Peltola, J., & Pättikangas, T. (2012).
        Development and validation of a boiling model
        for OpenFOAM multiphase solver.
        Proceedings of the CFD4NRS-4. p. 59.
        Daejeon, Democratic People's Republic of Korea, September 10-12, 2012.
    \endverbatim


    Alternatively a correlation can be used instead of the RPI wall boiling model.
    If the keyword nucleatingModel a model is provided the BC uses it
    instead of the RPI model.

    The transition boiling regime flux (TBF) is modelled following
    a temperature based linear interpolation between the critical heat flux
    (CHF) and the minimum heat flux (MHF) in such a way that when the wall
    temperature is between the range of TDBN and the Leidenfrost temperature
    (TLeiden) a linear interpolation is used between CHF and MHF.

    Thus, the following models are required:
      - LeidenfrostModel
      - CHFModel
      - CHFSubCoolModel
      - MHFModel
      - TDNBModel
      - filmBoilingModel

    The linear interpolation is as follows:

    \f[
        TBF = CHF*\phi + (1 - \phi)*MHF
    \f]

    with
    \f[
        \phi = w_p*(T_w - T_{DNB})/(T_{Leiden} - T_{DNB})
    \f]

    where:
    \vartable
        w_p | Model constant
        T_w | Wall temperature [K]
    \endvartable


    The film boiling regime is applied when \f$T_w\f$ is larger than
    \f$T_{Leiden}\f$. In this regime the correlation from the
    \c filmBoilingModel is used for calculating the cht from the wall.

    The \c filmBoilingModel is needed in the vapor field in order to calculate
    the heat transfer to the vapor phase in film boiling regime.


Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            compressible::alphatWallBoilingWallFunction;
        phaseType       <word>;
        otherPhase      <word>;
        relax           <Function1<scalar>>;

        partitioningModel
        {
            type        Lavieville;
            alphaCrit   0.2;
        }

        // Conditional entries

            // Option-1: phaseType=vapor

                // Optional entries
                LeidenfrostModel
                {
                    type        Spiegler;
                    Tcrit       647;
                }

                filmBoilingModel
                {
                    type        Bromley;
                }


            // Option-2: phaseType=liquid
            nucleationSiteModel
            {
                type        LemmertChawla;
            }

            departureDiamModel
            {
                type        TolubinskiKostanchuk;
            }

            departureFreqModel
            {
                type        Cole;
            }

                // Optional entries
                LeidenfrostModel
                {
                    type        Spiegler;
                    Tcrit       647;
                }

                CHFModel
                {
                    type        Zuber;
                }

                CHFSubCoolModel
                {
                    type        HuaXu;
                    Kburn       0.5;
                }

                MHFModel
                {
                    type        Jeschar;
                    Kmhf        1;
                }

                TDNBModel
                {
                    type        Schroeder;
                }

                filmBoilingModel
                {
                    type        Bromley;
                }

                dDep    <scalarField>;
                K       <scalar>;
                wp      <scalar>;
                qQuenching <scalarField>;


        // Optional entries
        alphatConv      <scalarField>;

        //Inherited entries
        ...
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type      | compressible::alphatWallBoilingWallFunction | word | yes | -
      phaseType | Name of phase type                     | word | yes  | -
      otherPhase | Name of other phase                   | word | yes  | -
      relax     | Relaxation factor for dmdt | Function1\<scalar\> <!--
                -->                                              | yes | -
      alphatConv | Convective turbulent thermal diffusivity    <!--
                -->                                 | scalarField | no | 0
      partitioningModel | Run-time selected heat flux partitioning model <!--
                -->                                      | dict | yes  | -
    \endtable

    Options for the \c phaseType and \c otherPhase entries:
    \verbatim
      vapor       | Vapor phase
      liquid      | Liquid phase
    \endverbatim

    when \c phaseType=liquid:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      nucleationSiteModel | Nucleation site density model | dict | yes  | -
      departureDiamModel  | Bubble departure diameter model     <!--
                        -->                               | dict | yes  | -
      departureFreqModel | Bubble departure frequency model | dict | yes | -
      LeidenfrostModel   | Leidenfrost temperature model  | dict | no   | -
      CHFModel | Critical heat flux model                 | dict | no   | -
      CHFSubCoolModel | CHF sub-cool model                | dict | no   | -
      MHFModel        | Minium heat flux model            | dict | no   | -
      TDNBModel       | Departure from nulceate boiling model | dict | no | -
      filmBoilingModel | Film boiling model               | dict | no   | -
      K | Model constant for area of bubbles            | scalar | no   | 4.0
      wp | Wetting parameter for transient boiling      | scalar | no   | 1.0
    \endtable

    The inherited entries are elaborated in:
     -\link alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H\endlink

Notes
  - Runtime selectabale submodels may require model specific entries
  - \c phaseType and \c otherPhase entries should be the opposite of each other.

See also
    Foam::alphatPhaseChangeJayatillekeWallFunctionFvPatchField

SourceFiles
    alphatWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H
#define compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H

#include "Function1.H"

#include "alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"
#include "nucleateFluxModel.H"

#include "LeidenfrostModel.H"
#include "filmBoilingModel.H"
#include "CHFModel.H"
#include "CHFSubCoolModel.H"
#include "MHFModel.H"
#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class alphatWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
{
public:

    // Public Enumerations

        //- Enumeration listing the possible operational modes
        enum phaseType
        {
            vaporPhase,
            liquidPhase
        };


private:

    // Private Data

        //- Enumeration of regimes per face
        enum regimeType
        {
            subcool,
            transient,
            film,
            nonBoiling
        };

        //- Name of the other phase (vapor/liquid phase)
        word otherPhaseName_;

        //- Names of heat source types
        static const Enum<phaseType> phaseTypeNames_;

        //- Heat source type
        phaseType phaseType_;

        //- Relaxation factor for dmdt
        autoPtr<Function1<scalar>> relax_;

        //- Patch face area by cell volume
        scalarField AbyV_;

        // Sub-cooling nucleating boiling

            //- Convective turbulent thermal diffusivity
            scalarField alphatConv_;

            //- Departure diameter field
            scalarField dDep_;

            //- Quenching surface heat flux
            scalarField qq_;

            //- Model constant for area of bubbles
            scalar K_;

            //- Run-time selected heat flux partitioning model
            autoPtr<wallBoilingModels::partitioningModel>
                partitioningModel_;

            //- Run-time selected nucleation site density model
            autoPtr<wallBoilingModels::nucleationSiteModel>
                nucleationSiteModel_;

            //- Run-time selected bubble departure diameter model
            autoPtr<wallBoilingModels::departureDiameterModel>
                departureDiamModel_;

            //- Run-time selected bubble departure frequency model
            autoPtr<wallBoilingModels::departureFrequencyModel>
                departureFreqModel_;

            //- Run-time sub-cooling heat flux correlatiom
            autoPtr<wallBoilingModels::nucleateFluxModel>
                nucleatingModel_;


        // Film boiling model

            //- Run-time selected for filmBoiling model
            autoPtr<wallBoilingModels::filmBoilingModel>
                filmBoilingModel_;


        // Transition boiling model

            //- Run-time selected for Leidenfrost temperature
            autoPtr<wallBoilingModels::LeidenfrostModel>
                LeidenfrostModel_;

            //- Run-time selected for CHF
            autoPtr<wallBoilingModels::CHFModel> CHFModel_;

            //- Run-time selected for CHF sub-cool
            autoPtr<wallBoilingModels::CHFSubCoolModel> CHFSoobModel_;

            //- Run-time selected for MHF
            autoPtr<wallBoilingModels::MHFModel> MHFModel_;

             //- Run-time selected for MHF
            autoPtr<wallBoilingModels::TDNBModel> TDNBModel_;

            //- Wetting parameter for transient boiling
            scalar wp_;

            //- Use Liquid temperature at y+=250
            bool liquidTatYplus_;

            //- Face regime
            labelField regimeTypes_;


public:

    //- Runtime type information
    TypeName("compressible::alphatWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- alphatWallBoilingWallFunctionFvPatchScalarField
        //- onto a new patch
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
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

        using alphatPhaseChangeWallFunctionFvPatchScalarField::dmdt;

        //- Is there phase change mass transfer for this phasePair
        virtual bool activePhasePair(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& dmdt(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& mDotL(const phasePairKey&) const;

        //- Return the departure diameter field
        const scalarField& dDeparture() const
        {
            return dDep_;
        }

        //- Return the quenching surface heat flux [W/m2]
        const scalarField& qq() const
        {
            return qq_;
        }

        //- Return the evaporation surface heat flux [W/m2]
        tmp<scalarField> qe() const
        {
            return mDotL_/AbyV_;
        }

        //- Return const reference to the face regime
        const labelField& regimeTypes() const noexcept
        {
            return regimeTypes_;
        }

        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
