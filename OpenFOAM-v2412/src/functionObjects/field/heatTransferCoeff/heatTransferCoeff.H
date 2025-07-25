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
    Foam::functionObjects::heatTransferCoeff

Group
    grpFieldFunctionObjects

Description
    Computes the heat transfer coefficient [W/(m^2 K)]
    as a \c volScalarField for a given set of patches.

    The heat transfer coefficient is definition dependent;
    therefore, different coefficient values could be
    obtained based on different definitions.

    Optionally, the Nusselt number (i.e. the ratio of convective
    to conductive heat transfer at a boundary in a fluid) can be output:

    \f[
        Nu = \frac{h L}{\kappa}
    \f]

    where
    \vartable
      Nu | Nusselt number
      h  | Heat transfer coefficient
      L  | Characteristic length that defines the scale of the physical system
      \kappa | Thermal conductivity of the fluid
    \endvartable

    Operands:
    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | -              | -
      output field      | volScalarField | \<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory entries
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       <word>;
        patches     (<wordRes>);
        htcModel    <word>;

        // Optional entries
        qr          <word>;
        L           <scalar>;
        kappa       <scalar>;

        // Conditional entries based on selected <htcModel>
        ...

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                        | Type | Reqd | Deflt
      type        | Type name: heatTransferCoeff       | word | yes  | -
      libs        | Library name: fieldFunctionObjects | word | yes  | -
      field       | Name of operand field              | word | yes  | -
      patches     | Names of operand patches        | wordRes | yes  | -
      htcModel    | Heat transfer coefficient model <!--
                -->                        - see below | word | yes  | -
      qr          | Name of radiative heat flux field  | word | no   | qr
      L           | Characteristic length that defines <!--
                --> the scale of the physical system | scalar | no   | 1
      kappa       | Thermal conductivity of fluid    | scalar | no   | 1
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link fieldExpression.H \endlink
      - \link heatTransferCoeffModel.H \endlink

    Options for the \c htcModel entry:
    \verbatim
      ReynoldsAnalogy           | Reynold's analogy
      localReferenceTemperature | Local reference temperature
      fixedReferenceTemperature | Fixed reference temperature
      faceZoneReferenceTemperature | Face-zone reference temperature
    \endverbatim

SourceFiles
    heatTransferCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_heatTransferCoeff_H
#define functionObjects_heatTransferCoeff_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class heatTransferCoeffModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class heatTransferCoeff Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoeff
:
    public fieldExpression
{
    // Private Data

        //- Characteristic length that defines the scale of the physical system
        scalar L_;

        //- Thermal conductivity of the fluid
        scalar kappa_;

        //- Heat transfer coefficient model
        autoPtr<heatTransferCoeffModel> htcModelPtr_;


protected:

        //- Calculate the heat transfer coefficient field and return true
        //- if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("heatTransferCoeff");


    // Constructors

        //- Construct from Time and dictionary
        heatTransferCoeff
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        heatTransferCoeff(const heatTransferCoeff&) = delete;

        //- No copy assignment
        void operator=(const heatTransferCoeff&) = delete;


    //- Destructor
    virtual ~heatTransferCoeff() = default;


    // Member Functions

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
