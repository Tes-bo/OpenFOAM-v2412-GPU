/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::HeatTransferModel

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Templated class to calculate the fluid-particle heat transfer
    coefficients based on a specified Nusselt-number model.

    \f[
        h = \frac{\mathrm{Nu} \, \kappa}{d_p}
    \f]
    where

    \vartable
        h            | Convective heat transfer coefficient of the flow
        \mathrm{Nu}  | Nusselt number
        \kappa       | Thermal conductivity of carrier in the film
        d_p          | Particle diameter
    \endvartable

    Optionally, Bird-Stewart-Lightfoot correction can be applied
    to correct the heat transfer coefficient for evaporation:

    \f[
        h_{corr} = h \, \frac{\beta}{ \exp(\beta) + 1 }
    \f]
    with

    \f[
        \beta = \frac{N \, C_p \, W}{h}
    \f]
    where

    \vartable
        \beta   | Correction factor
        N       | Molar flux
        C_p     | Specific heat capacity
        W       | Molecular weight
    \endvartable

    Reference:
    \verbatim
        Bird, R. B., Stewart, W. E., & Lightfoot, E. N. (1960).
        Transport phenomena.
        John Wiley & Sons., New York.
        DOI:10.1002/aic.690070245
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        heatTransferModel    <model>;

        <model>Coeffs
        {
            BirdCorrection   true;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description                | Type | Reqd   | Deflt
      heatTransferModel  | Type name: \<model\>       | word | yes    | -
      \<model\>Coeffs    | Model properties           | dict | cndtnl | -
      BirdCorrection     | Flag to apply Bird-Stewart-Lightfoot's <!--
      --> correction to the heat transfer coefficient | bool | cndtnl | -
    \endtable

    Options for the \c \<model\> entry:
    \verbatim
      RanzMarshall   | Ranz-Marshall correlation for Nusselt number
      none           | No active model
    \endverbatim

SourceFiles
    HeatTransferModel.C
    HeatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferModel_H
#define HeatTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class HeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HeatTransferModel
:
    public CloudSubModelBase<CloudType>
{
    // Private Data

        //- Flag to apply Bird-Stewart-Lightfoot's correction to the htc
        const Switch BirdCorrection_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HeatTransferModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Generated Methods

        //- No copy assignment
        void operator=(const HeatTransferModel<CloudType>&) = delete;


    // Constructors

        //- Construct null from owner
        HeatTransferModel(CloudType& owner);

        //- Construct from dictionary
        HeatTransferModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Copy construct
        HeatTransferModel(const HeatTransferModel<CloudType>& htm);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~HeatTransferModel() = default;


    //- Selector
    static autoPtr<HeatTransferModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- The Bird HTC correction flag
        bool BirdCorrection() const noexcept
        {
            return BirdCorrection_;
        }


        // Evaluation

            //- Return Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const = 0;

            //- Return heat transfer coefficient
            virtual scalar htc
            (
                const scalar dp,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHeatTransferModel(CloudType)                                       \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::HeatTransferModel<thermoCloudType>,                              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            HeatTransferModel<thermoCloudType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeHeatTransferModelType(SS, CloudType)                               \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);         \
                                                                               \
    Foam::HeatTransferModel<thermoCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType>>             \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeatTransferModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
