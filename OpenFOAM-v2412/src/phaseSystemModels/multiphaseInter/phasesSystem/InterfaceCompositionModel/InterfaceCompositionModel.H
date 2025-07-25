/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::InterfaceCompositionModel

Description
    Base class for interface composition models, templated on the two
    thermodynamic models either side of the interface.

SourceFiles
    InterfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef InterfaceCompositionModel_H
#define InterfaceCompositionModel_H

#include "interfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

using namespace multiphaseInter;

// Forward Declarations
class phaseModel;
class phasePair;
template <class ThermoType> class pureMixture;
template <class ThermoType> class multiComponentMixture;

/*---------------------------------------------------------------------------*\
                  Class InterfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class InterfaceCompositionModel
:
    public interfaceCompositionModel
{
protected:

    // Private data

        //- Thermo (from)
        const Thermo& fromThermo_;

        //- Other Thermo (to)
        const OtherThermo& toThermo_;

        //- Lewis number
        const dimensionedScalar Le_;


    // Private member functions

        //- Get a reference to the local thermo for a pure mixture
        template<class ThermoType>
        const typename pureMixture<ThermoType>::thermoType&
        getLocalThermo
        (
            const word& speciesName,
            const pureMixture<ThermoType>& globalThermo
        ) const;

        //- Get a reference to the local thermo for a multi component mixture
        template<class ThermoType>
        const typename multiComponentMixture<ThermoType>::thermoType&
        getLocalThermo
        (
            const word& speciesName,
            const multiComponentMixture<ThermoType>& globalThermo
        ) const;

        //- Return mass fraction for a pureMixture equal to one
        template<class ThermoType>
        tmp<volScalarField> getSpecieMassFraction
        (
            const word& speciesName,
            const pureMixture<ThermoType>& thermo
        ) const;

        //- Return mass fraction for speciesName
        template<class ThermoType>
        tmp<volScalarField> getSpecieMassFraction
        (
            const word& speciesName,
            const multiComponentMixture<ThermoType>& thermo
        ) const;

        //- Return moleculas weight of the mixture for pureMixture [Kg/mol]
        template<class ThermoType>
        tmp<volScalarField> MwMixture
        (
            const pureMixture<ThermoType>& thermo
        ) const;

        //- Return moleculas weight of the mixture for multiComponentMixture
        // [Kg/mol]
        template<class ThermoType>
        tmp<volScalarField> MwMixture
        (
            const multiComponentMixture<ThermoType>&
        ) const;


public:

    //- Construct from components
    InterfaceCompositionModel(const dictionary& dict, const phasePair& pair);

    //- Destructor
    virtual ~InterfaceCompositionModel() = default;


    // Member Functions

        //- Mass fraction difference between the interface and the field
        virtual tmp<volScalarField> dY
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- Reference mass fraction for species based models
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- Specie mass diffusivity for pure mixture
        virtual tmp<volScalarField> Dfrom
        (
            const word& speciesName
        ) const;

        //- Specie mass diffusivity for specie in a multicomponent
        virtual tmp<volScalarField> Dto
        (
             const word& speciesName
        ) const;

        //- Latent heat (to - from)(thermo - otherThermo)
        virtual tmp<volScalarField> L
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiation for multi-component (from) to single-component (to)
#define makeInterfaceDispSpecieMixtureType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, SpecieMixture<Mix<Phys>>>                             \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, OtherMix<OtherPhys>>                        \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addInterfaceCompositionToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )


// Instantiation for single-component (from) to multi-component (to)
#define makeInterfaceContSpecieMixtureType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, Mix<Phys>>                                            \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, SpecieMixture<OtherMix<OtherPhys>>>         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addInterfaceCompositionToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )


// Instantiation for single-component-single-component pairs
#define makeInterfacePureType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, Mix<Phys>>                                            \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, OtherMix<OtherPhys>>                        \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addInterfaceCompositionToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )


// Instantiation for multi-component-multi-component pairs
#define makeSpecieInterfaceSpecieMixtures(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, SpecieMixture<Mix<Phys>>>                             \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, SpecieMixture<OtherMix<OtherPhys>>>         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addInterfaceCompositionToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )


// Addition to the run-time selection table
#define addInterfaceCompositionToRunTimeSelectionTable(Type, Thermo, OtherThermo)\
                                                                               \
    typedef Type<Thermo, OtherThermo>                                          \
        Type##Thermo##OtherThermo;                                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Type##Thermo##OtherThermo,                                             \
        (                                                                      \
            word(Type##Thermo##OtherThermo::typeName_()) + "<"                 \
          + word(Thermo::typeName) + ","                                       \
          + word(OtherThermo::typeName) + ">"                                  \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        interfaceCompositionModel,                                             \
        Type##Thermo##OtherThermo,                                             \
        dictionary                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InterfaceCompositionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
