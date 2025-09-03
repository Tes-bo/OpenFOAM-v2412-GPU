/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::functionObjects::reactingEulerHtcModel

Description
    A heat transfer coefficient for reactingEuler solvers.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    reactingEulerHtcModel1
    {
        // Mandatory entries
        type        reactingEulerHtcModel;
        libs        (fieldFunctionObjects);

        field       <word>;
        patches     (<wordRes>);
        htcModel    <word>;

        // Optional entries
        qr          <word>;

        // Conditional entries based on selected <htcModel>
        ...

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                        | Type | Reqd | Deflt
      type        | Type name: reactingEulerHtcModel   | word | yes  | -
      libs        | Library name: fieldFunctionObjects | word | yes  | -
      field       | Name of operand field              | word | yes  | -
      patches     | Names of operand patches        | wordRes | yes  | -
      htcModel    | Heat transfer coefficient model <!--
                -->                        - see below | word | yes  | -
      qr          | Name of radiative heat flux field  | word | no   | qr
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
    reactingEulerHtcModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_reactingEulerHtcModel_H
#define Foam_functionObjects_reactingEulerHtcModel_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class heatTransferCoeffModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class reactingEulerHtcModel Declaration
\*---------------------------------------------------------------------------*/

class reactingEulerHtcModel
:
    public fieldExpression
{
    // Private Data

        //- Heat transfer coefficient model
        autoPtr<heatTransferCoeffModel> htcModelPtr_;


protected:

    // Protected Member Functions

        //- Calculate the heat transfer coefficient field
        //  \return true on success
        virtual bool calc();

        //- Calculate heat flux
        tmp<FieldField<Field, scalar>> q() const;

        //- No copy construct
        reactingEulerHtcModel(const reactingEulerHtcModel&) = delete;

        //- No copy assignment
        void operator=(const reactingEulerHtcModel&) = delete;


public:

    //- Runtime type information
    TypeName("reactingEulerHtcModel");


    // Constructors

        //- Construct from components
        reactingEulerHtcModel
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reactingEulerHtcModel() = default;


    // Member Functions

        //- Read the heatTransferCoeff data
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
