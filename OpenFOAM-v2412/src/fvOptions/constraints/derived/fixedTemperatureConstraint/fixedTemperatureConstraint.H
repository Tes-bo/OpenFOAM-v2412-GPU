/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::fv::fixedTemperatureConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain temperature equation (i.e. \c T) with
    a given set of fixed values within a specified region.

    Constraints applied to:
    \verbatim
      T    | Temperature                   [K]
    \endverbatim

    Required fields:
    \verbatim
      T    | Temperature                   [K]
      e/h  | Internal energy/Enthalphy     [m2/s2]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    fixedTemperatureConstraint1
    {
        // Mandatory entries (unmodifiable)
        type            fixedTemperatureConstraint;
        mode            <mode>;

        // Conditional mandatory/optional entries (runtime modifiable)

            // when mode=uniform (mandatory)
            temperature     constant 500; // fixed temperature with time [K]

            // when mode=lookup (optional)
            T            <Tname>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd  | Dflt
      type      | Type name: fixedTemperatureConstraint | word   | yes   | -
      mode      | Temperature mode - see below          | word   | yes   | -
      temperature | Function1 type temperature field    | Function1 | cndtnl | -
      T         | Name of operand temperature field     | word   | cndtnl | T
    \endtable

    Options for the \c mode entry:
    \verbatim
      uniform | Use Function1 type temperature field
      lookup  | Use temperature field given field name
    \endverbatim

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note:
  - The \c uniform option allows the use of a time-varying
    uniform temperature by means of the \c Function1 type.

See also
  - Foam::fv::limitTemperature
  - Foam::Function1

SourceFiles
    fixedTemperatureConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureConstraint_H
#define fixedTemperatureConstraint_H

#include "cellSetOption.H"
#include "Enum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class fixedTemperatureConstraint Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureConstraint
:
    public fv::cellSetOption
{
public:

    //- Temperature mode
    enum temperatureMode
    {
        tmUniform,
        tmLookup
    };

    //- String representation of temperatureMode enums
    static const Enum<temperatureMode> temperatureModeNames_;


protected:

    // Protected Data

        //- Operation mode
        temperatureMode mode_;

        //- Uniform temperature [K]
        autoPtr<Function1<scalar>> Tuniform_;

        //- Temperature field name
        word TName_;


public:

    //- Runtime type information
    TypeName("fixedTemperatureConstraint");


    // Constructors

        //- Construct from components
        fixedTemperatureConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        fixedTemperatureConstraint(const fixedTemperatureConstraint&) = delete;

        //- No copy assignment
        void operator=(const fixedTemperatureConstraint&) = delete;


    //- Destructor
    virtual ~fixedTemperatureConstraint() = default;


    // Member Functions

        //- Constrain energy equation to fix the temperature
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
