/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Thorsten Zirwes
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::BilgerMixtureFraction

Group
    grpThermophysicalFunctionObjects

Description
    Calculates the Bilger mixture-fraction field (i.e.
    \f$f_{\mathrm{Bilger}}\f$) based on the elemental composition
    of the mixture. Elements \c C, \c H, \c S and \c O are considered.

    \f$f_{\mathrm{Bilger}}\f$ is the mass mixing ratio of fuel and oxidiser
    in [kg fuel / kg mixture], and is invariant to the reaction progress of
    the mixture (e.g. the same for unburnt and burnt mixtures).

    \f$f_{\mathrm{Bilger}}\f$ equals to the unity
    for pure fuel and to zero for pure oxidiser.

    The Bilger mixture-fraction field is computed based on the following:

    \f[
        f_{\mathrm{Bilger}} =
            \frac{\beta - \beta_{\mathrm{ox}}}{\beta_{\mathrm{fuel}}
          - \beta_{\mathrm{ox}}}
    \f]

    with

    \f[
        \beta =
            2\frac{Y_C}{W_C}
          + 2\frac{Y_S}{W_S}
          + \frac{1}{2}\frac{Y_H}{W_H}
          - \frac{Y_O}{W_O}
    \f]

    where
    \vartable
      \beta               | Coupling function                        [kmol/kg]
      Y_e                 | Elemental mass fraction of element, e
      W_e                 | Atomic weight of element, e
      {.}_{\mathrm{ox}}   | Subscript to denote oxidiser composition
      {.}_{\mathrm{fuel}} | Subscript to denote fuel composition
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    References:
    \verbatim
        Original method:
            Bilger, R. W. (1979).
            Turbulent jet diffusion flames.
            Energy and Combustion Science, p. 109-131.
            DOI:10.1016/B978-0-08-024780-9.50011-3

        Implementation:
            Zirwes, T., Zhang, F., Habisreuther, P., Hansinger, M.,
            Bockhorn, H., Pfitzner, M., & Trimis, D. (2019).
            Quasi-DNS dataset of a piloted flame
            with inhomogeneous inlet conditions.
            Flow, Turbulence and Combustion, vol 104, p. 997-1027.
            DOI:10.1007/s10494-019-00081-5
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    BilgerMixtureFraction1
    {
        // Mandatory entries (unmodifiable)
        type                 BilgerMixtureFraction;
        libs                 (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fuel
        {
            // Optional entries (runtime modifiable)
            fractionBasis    mass;

            // Conditional mandatory entries (runtime modifiable)
            CH4              1;
        }

        oxidiser
        {
            // Optional entries (runtime modifiable)
            fractionBasis    mole;

            // Conditional mandatory entries (runtime modifiable)
            O2               0.23;
            N2               0.77;
        }

        // Optional entries (runtime modifiable)
        phase                <phaseName>;
        result               <resultName>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Dflt
      type     | Type name: BilgerMixtureFraction    | word | yes  | -
      libs     | Library name: fieldFunctionObjects  | word | yes  | -
      fuel     | Dictionary for fuel composition     | dict | yes  | -
      oxidiser | Dictionary for oxidiser composition | dict | yes  | -
      phase    | Name of phase (e.g. "gas")          | word | no   | ""
      result   | Name of resulting field             | word | no   | f_Bilger
      fractionBasis | Species-fraction interpretation method | word | no | mass
    \endtable

    Options for the \c fractionBasis entry:
    \verbatim
      mass   | Interpret species fractions as mass fractions
      mole   | Interpret species fractions as mole fractions
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

Note
  - The mole or mass fractions are automatically normalized to the unity.

See also
  - Foam::functionObject
  - Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    BilgerMixtureFraction.C

\*---------------------------------------------------------------------------*/

#ifndef BilgerMixtureFraction_H
#define BilgerMixtureFraction_H

#include "fvMeshFunctionObject.H"
#include "specieElement.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
/*---------------------------------------------------------------------------*\
                    Class BilgerMixtureFraction Declaration
\*---------------------------------------------------------------------------*/

class BilgerMixtureFraction
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the phase (e.g. "gas" for multiphase applications)
        word phaseName_;

        //- Name of the resulting mixture-fraction field
        word resultName_;

        //- Reference to thermo object
        const basicSpecieMixture& thermo_;

        //- Number of species
        label nSpecies_;

        // Amount of oxygen required to fully oxidise the oxidiser
        scalar o2RequiredOx_;

        // Amount of oxygen required to oxidise the fuel minus the oxidiser
        scalar o2RequiredFuelOx_;

        //- Number of carbon atoms for each species
        labelField nAtomsC_;

        //- Number of sulphur atoms for each species
        labelField nAtomsS_;

        //- Number of hydrogen atoms for each species
        labelField nAtomsH_;

        //- Number of oxygen atoms for each species
        labelField nAtomsO_;

        //- Mass fractions of species in the oxidiser
        scalarField Yoxidiser_;

        //- Mass fractions of species in the fuel
        scalarField Yfuel_;


    // Private Member Functions

        //- Calculate the Bilger mixture-fraction
        void calcBilgerMixtureFraction();

        //- Read composition of fuel and oxidiser from subdictionary
        bool readComposition
        (
            const dictionary& subDict,
            scalarField& comp
        );

        //- Compute amount of oxygen required to oxidise a mixture
        scalar o2Present(const scalarField&) const;

        //- Compute amount of oxygen present in a mixture
        scalar o2Required(const scalarField&) const;


public:

    //- Runtime type information
    TypeName("BilgerMixtureFraction");


    // Constructors

        //- Construct from Time and dictionary
        BilgerMixtureFraction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        BilgerMixtureFraction(const BilgerMixtureFraction&) = delete;

        //- No copy assignment
        void operator=(const BilgerMixtureFraction&) = delete;


    //- Destructor
    virtual ~BilgerMixtureFraction() = default;


    // Member Functions

        //- Read the BilgerMixtureFraction data
        virtual bool read(const dictionary&);

        //- Calculate the Bilger mixture-fraction field
        virtual bool execute();

        //- Clear the Bilger mixture-fraction field from registry
        virtual bool clear();

        //- Write Bilger mixture-fraction field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
