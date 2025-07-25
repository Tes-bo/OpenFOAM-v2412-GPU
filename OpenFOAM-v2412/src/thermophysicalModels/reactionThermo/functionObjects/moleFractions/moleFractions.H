/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::moleFractions

Group
    grpThermophysicalFunctionObjects

Description
    Calculates mole-fraction fields from the mass-fraction
    fields of the psi/rhoReactionThermo and caches them
    for output and further post-processing.

    The names of the mole-fraction fields are obtained from
    the corresponding mass-fraction fields prepended by "X_".

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    moleFractions1
    {
        // Conditional mandatory entries (unmodifiable)
        // Either of the below depending on
        // the thermodynamics package used in the solver.

            // Option-1
            type                 psiReactionThermoMoleFractions;

            // Option-2
            type                 rhoReactionThermoMoleFractions;

        // Mandatory entries (unmodifiable)
        libs                 (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        phase                <phaseName>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Dflt
      type     | Type name: psiReactionThermoMoleFractions or <!--
                --> rhoReactionThermoMoleFractions   | word | yes  | -
      libs     | Library name: fieldFunctionObjects  | word | yes  | -
      phase    | Name of phase (e.g. "gas")          | word | no   | ""
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    moleFractions.C
    moleFractionsFunctionObjects.C

\*---------------------------------------------------------------------------*/

#ifndef moleFractions_H
#define moleFractions_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class moleFractions Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class moleFractions
:
    public functionObjects::fvMeshFunctionObject
{
    // Private Data

        //- Species mole fractions
        PtrList<volScalarField> X_;

        //- Name of phase
        word phaseName_;


    // Private Member Functions

        //- Calculate the mole fraction fields
        virtual void calcMoleFractions();


public:

    //- Runtime type information
    TypeName("moleFractions");


    // Constructors

        //- Construct from Time and dictionary
        moleFractions
        (
            const word& name,
            const Time& t,
            const dictionary& dict
        );

        //- No copy construct
        moleFractions(const moleFractions&) = delete;

        //- No copy assignment
        void operator=(const moleFractions&) = delete;


    //- Destructor
    virtual ~moleFractions() = default;


    // Member Functions

        //- Read the moleFractions data
        virtual bool read(const dictionary& dict);

        //- Calculate the mole-fraction fields
        virtual bool execute();

        //- The mole-fraction fields auto-write - no-op
        virtual bool write()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "moleFractions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
