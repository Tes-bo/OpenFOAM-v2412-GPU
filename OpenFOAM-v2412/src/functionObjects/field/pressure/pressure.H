/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::pressure

Group
    grpFieldFunctionObjects

Description
    Provides several methods to convert an input pressure
    field into derived forms, including:

    - static pressure
        \f[
            p_s = p_{ref} + \rho p_k
        \f]
    - total pressure
        \f[
            p_0 = p_{ref} + p + 0.5 \rho |\vec U|^2
        \f]
    - isentropic pressure
        \f[
            p_i = p*(1 + ((\gamma-1)*M^2)/2)^{(\gamma/(\gamma - 1))}
        \f]
    - static pressure coefficient
        \f[
            Cp = \frac{p_s - p_{\inf}}{0.5 \rho_{\inf} |\vec U_{\inf}|^2}
        \f]
    - total pressure coefficient
        \f[
            Cp_0 = \frac{p_0 - p_{\inf}}{0.5 \rho_{\inf} |\vec U_{\inf}|^2}
        \f]

    where
    \vartable
      \rho        | Density [kg/m3]
      \vec U      | Velocity [m/s]
      \rho_{\inf} | Freestream density [kg/m3]
      p_{\inf}    | Freestream pressure [Pa]
      U_{\inf}    | Freestream velocity [m/s]
      p_k         | Kinematic pressure (p/rho)[m2/s2]
      p_s         | Static pressure [Pa]
      p_0         | Total pressure [Pa]
      p_{ref}     | Reference pressure level [Pa]
      p_i         | Total isentropic pressure
      Cp          | Pressure coefficient
      Cp_0        | Total pressure coefficient
      \gamma      | Specific heat ratio
    \endvartable

    The function object will operate on both kinematic (\f$ p_k \f$) and static
    pressure (\f$ p \f$) fields.

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
    pressure1
    {
        // Mandatory entries (unmodifiable)
        type            pressure;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        mode            \<option\>;

        // Optional entries (runtime modifiable)
        field           \<pName\>;
        U               \<UName\>;
        rho             <rhoName>;
        rhoInf          1.0; // enabled if rho=rhoInf
        pRef            0.0;
        hydroStaticMode none;
        g               (0 -9.81 0); // enabled if hydroStaticMode != none
        hRef            0.0;         // enabled if hydroStaticMode != none
        pInf            0.0;
        UInf            (1 0 0);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: pressure                | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      mode         | Calculation mode (see below)       | word |  yes  | -
      field        | Name of the pressure field         | word |  no   | p
      U            | Name of the velocity field         | word |  no   | U
      rho          | Name of the density field          | word |  no   | rho
      rhoInf | Freestream density for coefficient calculation | scalar <!--
         --> | conditional| -
      pRef         | Reference pressure for total pressure    | scalar | no | 0
      hydrostaticMode | Hydrostatic contributions (see below) | word | no | none
      g            | Gravity vector (see below)               | vector | no | -
      hRef         | Reference height (see below)             | scalar | no | -
      pInf | Freestream pressure for coefficient calculation  | scalar | no | -
      UInf | Freestream velocity for coefficient calculation  | vector | no | -
    \endtable

    Options for the \c mode entry:
    \verbatim
        static      | static pressure
        total       | total pressure
        isentropic  | isentropic pressure
        staticCoeff | static pressure coefficient
        totalCoeff  | total pressure coefficient
    \endverbatim

    The optional \c hydrostaticMode entry provides handling for the term
    \f$ \rho (\vec{g} \dot \vec{h})\f$ where options include
    \verbatim
        none     | not included
        add      | add the term, e.g. to convert from p_rgh to p
        subtract | subtract the term, e.g. to convert from p to p_rgh
    \endverbatim

    If the \c hydrostaticMode is active, values are also required for
    gravity, \c g, and reference height, \c hRef.  By default these will be
    retrieved from the database.  When these values are not available
    the user must provide them, e.g.
    \verbatim
        g           (0 -9.81 0);
        hRef        0;
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::pressure

SourceFiles
    pressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressure_H
#define functionObjects_pressure_H

#include "fieldExpression.H"
#include "volFieldsFwd.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class pressure Declaration
\*---------------------------------------------------------------------------*/

class pressure
:
    public fieldExpression
{
public:

    // Public Data Types

        //- Enumeration for pressure calculation mode
        enum mode : unsigned
        {
            STATIC       = (1 << 0), //!< Static pressure
            TOTAL        = (1 << 1), //!< Total pressure
            ISENTROPIC   = (1 << 2), //!< Isentropic pressure
            COEFF        = (1 << 3), //!< Coefficient manipulator
            STATIC_COEFF = (STATIC | COEFF),
            TOTAL_COEFF  = (TOTAL | COEFF),
        };

        static const Enum<mode> modeNames;

        //- Enumeration for hydrostatic contributions
        enum hydrostaticMode : unsigned
        {
            NONE = 0,
            ADD,
            SUBTRACT
        };

        static const Enum<hydrostaticMode> hydrostaticModeNames;


private:

    // Private Data

        //- Calculation mode
        mode mode_;

        //- Hydrostatic constribution mode
        hydrostaticMode hydrostaticMode_;

        //- Name of velocity field
        word UName_;

        //- Name of density field
        word rhoName_;


        // Total pressure calculation

            //- Reference pressure level
            scalar pRef_;


        // Pressure coefficient calculation

            //- Freestream pressure
            scalar pInf_;

            //- Freestream velocity
            vector UInf_;

            //- Freestream density
            scalar rhoInf_;

            //- Flag to show whether rhoInf has been initialised
            bool rhoInfInitialised_;


        //- p +/- rgh calculation

            //- Gravity vector
            mutable dimensionedVector g_;

            //- Flag to show whether g has been initialised
            bool gInitialised_;

            //- Reference height
            mutable dimensionedScalar hRef_;

            //- Flag to show whether hRef has been initialised
            bool hRefInitialised_;


    // Private Member Functions

        //- Return the name of the derived pressure field
        word resultName() const;

        //- Multiply the static pressure p by rhoInf if necessary and return
        tmp<volScalarField> rhoScale(const volScalarField& p) const;

        //- Multiply the given field by rho or rhoInf as appropriate and return
        tmp<volScalarField> rhoScale
        (
            const volScalarField& p,
            const tmp<volScalarField>& tsf
        ) const;

        //- Add the hydrostatic contribution
        void addHydrostaticContribution
        (
            const volScalarField& p,
            volScalarField& prgh
        ) const;

        //- Calculate and return the pressure
        tmp<volScalarField> calcPressure
        (
            const volScalarField& p,
            const tmp<volScalarField>& tp
        ) const;

        //- Convert to coeff by applying the freestream dynamic pressure scaling
        tmp<volScalarField> coeff(const tmp<volScalarField>& tp) const;

        //- Calculate the derived pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pressure");


    // Constructors

        //- Construct from Time and dictionary
        pressure
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        pressure(const pressure&) = delete;

        //- No copy assignment
        void operator=(const pressure&) = delete;


    //- Destructor
    virtual ~pressure() = default;


    // Member Functions

        //- Read the pressure data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
