/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::functionObjects::proudmanAcousticPower

Group
    grpFieldFunctionObjects

Description
    Computes the acoustic power due to the volume of isotropic turbulence
    using Proudman's formula.

    The acoustic power, i.e. \f$ P_A \f$ [\f$W/m^3\f$], in terms of turbulent
    kinetic energy, i.e. \f$ k \f$, and turbulent kinetic energy dissipation
    rate, i.e. \f$ \epsilon \f$, is given as:

        \f[
            P_A = \alpha_\epsilon \rho \epsilon M_t^5
        \f]

    where \f$ \alpha_\epsilon = 0.1 \f$ is a constant and

        \f[
            M_t = \frac{\sqrt{2 k}}{a_0}
        \f]

    with \f$ a_0 \f$ the speed of sound.  The acoustic power is also output in
    dB using:

        \f[
            L_P = 10 \log \frac{P_A}{P_{ref}}
        \f]

    where \f$ P_{ref} = 1e^{-12} \f$ [\f$W/m^3\f$] is a constant.

    Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    proudmanAcousticPower1
    {
        // Mandatory entries (unmodifiable)
        type        proudmanAcousticPower;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        alphaEps    0.1;
        // For incompressible flow simulations
        rhoInf      1.225;
        aRef        340;
        // Turbulence field names (if not retrieved from the turb model)
        k           kMean;
        epsilon     epsilonMean;
        omega       none; // omegaMean

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                         | Type | Req'd | Dflt
      type         | Type name: proudmanAcousticPower    | word |  yes  | -
      libs         | Library name: fieldFunctionObjects  | word |  yes  | -
      rhoInf       | Freestream density (for incompressible) | scalar <!--
               --> | conditional | -
      aRef         | Speed of sound (incompressible)     | scalar <!--
               --> | conditional | -
      alphaEps     | Empirical model coefficient         | scalar | no  | 0.1
      k            | Turbulence k field name             | word | no    | none
      epsilon      | Turbulence epsilon field name       | word | no    | none
      omega        | Turbulence omega field name         | word | no    | none
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    The freestream density and reference speed of sound are only necessary
    when a thermodynamics package is unavailable, typically for incompressible
    cases.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::proudmanAcousticPower

SourceFiles
    proudmanAcousticPower.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_proudmanAcousticPower_H
#define functionObjects_proudmanAcousticPower_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class proudmanAcousticPower Declaration
\*---------------------------------------------------------------------------*/

class proudmanAcousticPower
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Empirical model coefficient
        scalar alphaEps_;

        //- Freestream density (incompressible calcs only)
        dimensionedScalar rhoInf_;

        //- Reference speed of sound (incompressible calcs only)
        dimensionedScalar aRef_;

        //- Name of turbulence k field; default = none
        word kName_;

        //- Name of turbulence epsilon field; default = none
        word epsilonName_;

        //- Name of turbulence omega field; default = none
        word omegaName_;


    // Private Member Functions

        //- Multiply the field by density and return
        tmp<volScalarField> rhoScale(const tmp<volScalarField>& fld) const;

        //- Speed of sound
        tmp<volScalarField> a() const;

        //- Turbulence kinetic energy dissipation rate
        tmp<volScalarField> k() const;

        //- Turbulence dissipation
        tmp<volScalarField> epsilon() const;


public:

    //- Runtime type information
    TypeName("proudmanAcousticPower");


    // Constructors

        //- Construct from Time and dictionary
        proudmanAcousticPower
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        proudmanAcousticPower(const proudmanAcousticPower&) = delete;

        //- No copy assignment
        void operator=(const proudmanAcousticPower&) = delete;


    //- Destructor
    virtual ~proudmanAcousticPower() = default;


    // Member Functions

        //- Read the Proudman acoustic power data
        virtual bool read(const dictionary&);

        //- Calculate the Proudman acoustic power
        virtual bool execute();

        //- Write the Proudman acoustic power
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
