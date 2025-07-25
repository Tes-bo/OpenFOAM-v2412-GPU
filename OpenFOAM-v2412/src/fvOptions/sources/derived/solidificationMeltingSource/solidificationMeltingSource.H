/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::fv::solidificationMeltingSource

Group
    grpFvOptionsSources

Description
    This source is designed to model the effect of solidification and
    melting processes, e.g. windhield defrosting, within a specified region.
    The phase change occurs at the melting temperature, \c Tmelt.

    The presence of the solid phase in the flow field is incorporated into the
    model as a momentum porosity contribution; the energy associated with the
    phase change is added as an enthalpy contribution.

    References:
    \verbatim
        Voller, V. R., & Prakash, C. (1987).
        A fixed grid numerical modelling methodology for
        convection-diffusion mushy region phase-change problems.
        International Journal of Heat and Mass Transfer, 30(8), 1709-1719.
        DOI:10.1016/0017-9310(87)90317-6

        Swaminathan, C. R., & Voller, V. R. (1992).
        A general enthalpy method for modeling solidification processes.
        Metallurgical transactions B, 23(5), 651-664.
        DOI:10.1007/BF02649725
    \endverbatim

    The model generates a field \c \<name\>:alpha1 which can be visualised to
    to show the melt distribution as a fraction [0-1].

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    solidificationMeltingSource1
    {
        // Mandatory entries (unmodifiable)
        type           solidificationMeltingSource;

        // Mandatory entries (runtime modifiable)
        Tmelt          273;
        L              334000;
        thermoMode     <thermoModeName>;
        rhoRef         800;
        beta           5e-6;

        // Optional entries (runtime modifiable)
        relax          0.9;
        T              <Tname>;
        rho            <rhoName>;
        U              <Uname>;
        phi            <phiName>;
        Cu             1e5;
        q              1e-2;

        // Conditional optional entries (runtime modifiable)

            // when thermoMode=lookup
            Cp         Cp;

        // Conditional mandatory entries (runtime modifiable)

            // when Cp=CpRef
            CpRef  1000;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd | Dflt
      type      | Type name: solidificationMeltingSource  | word   | yes   | -
      Tmelt     | Melting temperature [K]           | scalar | yes  | -
      L         | Latent heat of fusion [J/kg]      | scalar | yes  | -
      thermoMode  | Thermo mode                     | word   | yes  | -
      rhoRef    | Reference (solid) density         | scalar | yes  | -
      beta      | Thermal expansion coefficient [1/K] | scalar | yes | -
      relax     | Relaxation factor [0-1]           | scakar | no   | 0.9
      T         | Name of operand temperature field | word   | no   | T
      rho       | Name of operand density field     | word   | no   | rho
      U         | Name of operand velocity field    | word   | no   | U
      phi       | Name of operand flux field        | word   | no   | phi
      Cu        | Mushy region momentum sink coefficient [1/s] <!--
                -->                                 | scalar | no   | 1e5
      q         | Coefficient used in porosity calc | scalar | no   | 1e-2
      Cp        | Name of specific heat capacity field | word | cndtnl | Cp
      CpRef     | Specific heat capacity value      | scalar | cndtnl  | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

    Options for the \c thermoMode entry:
    \verbatim
      thermo    | Access Cp information from database
      lookup    | Access Cp information by looking up from dictionary
    \endverbatim

SourceFiles
    solidificationMeltingSource.C
    solidificationMeltingSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationMeltingSource_H
#define solidificationMeltingSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class solidificationMeltingSource Declaration
\*---------------------------------------------------------------------------*/

class solidificationMeltingSource
:
    public fv::cellSetOption
{
public:

    // Public Enumerations

        //- Options for the thermo mode specification
        enum thermoMode
        {
            mdThermo,
            mdLookup
        };

        //- Names for thermoMode
        static const Enum<thermoMode> thermoModeTypeNames_;


private:

    // Private Data

        //- Temperature at which melting occurs [K]
        scalar Tmelt_;

        //- Latent heat of fusion [J/kg]
        scalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Thermodynamics mode
        thermoMode mode_;

        //- Reference density - typically the solid density
        scalar rhoRef_;

        //- Name of operand temperature field
        word TName_;

        //- Name of specific heat capacity field
        word CpName_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand flux field
        word phiName_;

        //- Mushy region momentum sink coefficient [1/s]
        scalar Cu_;

        //- Coefficient used in porosity calculation
        scalar q_;

        //- Thermal expansion coefficient [1/K]
        scalar beta_;

        //- Phase fraction indicator field
        volScalarField alpha1_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Temperature change cached for source calculation when alpha1 updated
        scalarField deltaT_;


    // Private Member Functions

        //- Return the specific heat capacity field
        tmp<volScalarField> Cp() const;

        //- Update the model
        void update(const volScalarField& Cp);

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn);


public:

    //- Runtime type information
    TypeName("solidificationMeltingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        solidificationMeltingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        solidificationMeltingSource
        (
            const solidificationMeltingSource&
        ) = delete;

        //- No copy assignment
        void operator=(const solidificationMeltingSource&) = delete;


    //- Destructor
    ~solidificationMeltingSource() = default;


    // Member Functions

        //- Add explicit contribution to enthalpy equation
        virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi);

        //- Add implicit contribution to momentum equation
        virtual void addSup(fvMatrix<vector>& eqn, const label fieldi);

        //- Add explicit contribution to compressible enthalpy equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidificationMeltingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
