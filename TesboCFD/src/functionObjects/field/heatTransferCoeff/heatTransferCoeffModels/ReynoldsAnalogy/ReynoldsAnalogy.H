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
    Foam::heatTransferCoeffModels::ReynoldsAnalogy

Description
    Heat transfer coefficient calculation based on Reynolds Analogy,
    which is used to relate turbulent momentum and heat transfer.

    The heat transfer coefficient is derived
    from the skin friction coefficient:

    \f[
        C_f = \frac{\tau_w}{0.5 \rho_{ref} |U|^2}
    \f]

    as:

    \f[
        h = 0.5 \rho_{ref} c_{p,ref} |U_{ref}| C_f
    \f]

    where
    \vartable
      h       | Heat transfer coefficient                          [W/m^2/K]
      \rho_{ref}  | Reference fluid density                        [kg/m^3]
      c_{p,ref} | Reference specific heat capacity at constant pressure [J/kg/K]
      U_{ref}     | Reference velocity                             [m/s]
      C_f         | Skin friction coefficient                      [-]
      \tau_w      | Wall shear stress                              [m^2/s^2]
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Inherited entries
        ...

        // Mandatory entries
        htcModel        ReynoldsAnalogy;
        UInf            <vector>;

        // Optional entries
        U               <word>;
        Cp              <word>;
        rho             <word>;

        // Conditional mandatory entries

            // when Cp == CpInf
            CpInf       <scalar>;

            // when rho == rhoInf
            rhoInf      <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type    | Model name: ReynoldsAnalogy              | word | yes  | -
      UInf    | Reference velocity                     | vector | yes  | -
      U       | Name of velocity field                   | word | no   | U
      Cp      | Name of reference specific heat capacity | word | no   | Cp
      CpInf   | Reference specific heat capacity value | scalar | choice | -
      rho     | Name of fluid density field              | word | no   | rho
      rhoInf  | Reference fluid density value          | scalar | choice | -
    \endtable

Note
  - In order to use a reference \c Cp, set \c Cp to \c CpInf.
  - In order to use a reference \c rho, set \c rho to \c rhoInf.

SourceFiles
    ReynoldsAnalogy.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_ReynoldsAnalogy_H
#define heatTransferCoeffModels_ReynoldsAnalogy_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                        Class ReynoldsAnalogy Declaration
\*---------------------------------------------------------------------------*/

class ReynoldsAnalogy
:
    public heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Name of velocity field
        word UName_;

        //- Reference velocity
        vector URef_;

        //- Name of fluid density field
        word rhoName_;

        //- Reference fluid density
        scalar rhoRef_;

        //- Name of specific heat capacity field
        word CpName_;

        //- Reference specific heat capacity
        scalar CpRef_;


    // Protected Member Functions

        //- Return fluid density field [kg/m^3]
        virtual tmp<scalarField> rho(const label patchi) const;

        //- Return heat capacity at constant pressure [J/kg/K]
        virtual tmp<scalarField> Cp(const label patchi) const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return skin friction coefficient field [-]
        tmp<FieldField<Field, scalar>> Cf() const;

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );


public:

    //- Runtime type information
    TypeName("ReynoldsAnalogy");


    // Constructors

        //- Construct from components
        ReynoldsAnalogy
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        ReynoldsAnalogy(const ReynoldsAnalogy&) = delete;

        //- No copy assignment
        void operator=(const ReynoldsAnalogy&) = delete;


    //- Destructor
    virtual ~ReynoldsAnalogy() = default;


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
