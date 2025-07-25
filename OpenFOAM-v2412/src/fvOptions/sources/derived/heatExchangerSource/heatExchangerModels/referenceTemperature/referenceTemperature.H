/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::heatExchangerModels::referenceTemperature

Description
    A heat exchanger model where heat exchange
    is calculated via a temperature table.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <name>
    {
        // Inherited entries
        ...

        // Mandatory entries
        model                    referenceTemperature;
        targetQdot               <Function1<scalar>>;

        // Conditional mandatory entries

            // Option-1
            Tref                     <scalar>;

            // Option-2
            file                     "<word>";
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Deflt
      model     | Model name:effectivenessTable     | word  | yes  | -
      targetQdot | Target heat rejection [J/s]                    <!--
                -->                    | Function1\<scalar\> | yes | -
      Tref      | Reference temperature [K]         | scalar | choice | -
      file      | 2D reference temperature table = function of primary    <!--
                --> mass flow rate [kg/s] and primary flow temperature    <!--
                -->                                 | file   | choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link heatExchangerSource.H \endlink
      - \link interpolation2DTable.H \endlink

SourceFiles
    referenceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_heatExchangerModels_referenceTemperature_H
#define Foam_heatExchangerModels_referenceTemperature_H

#include "heatExchangerModel.H"
#include "interpolation2DTable.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatExchangerModels
{

/*---------------------------------------------------------------------------*\
                    Class referenceTemperature Declaration
\*---------------------------------------------------------------------------*/

class referenceTemperature
:
    public heatExchangerModel
{
    // Private Data

        //- Target heat rejection [J/s]
        autoPtr<Function1<scalar>> targetQdotPtr_;

        //- 2D reference temperature table
        //- Function of primary mass flow rate [kg/s]
        //- and primary flow temperature [K]
        autoPtr<interpolation2DTable<scalar>> TrefTablePtr_;

        //- Net mass flux [kg/s]
        scalar sumPhi_;

        //- Total heat exchange [J/s]
        scalar Qt_;

        //- Reference temperature [K]
        scalar Tref_;


    // Private Member Functions

    // Evaluation

        //- Return primary-flow net mass flux
        scalar primaryNetMassFlux() const;

        //- Return primary-flow inlet temperature
        scalar primaryInletTemperature() const;

        //- Return temperature differences
        scalarField temperatureDifferences(const labelList& cells) const;

        //- Return interim weights
        scalar weight
        (
            const labelList& cells,
            const scalarField& deltaTCells
        ) const;


    // I-O

        //- Write file header information
        void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("referenceTemperature");


    // Constructors

        //- Construct from components
        referenceTemperature
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& coeffs
        );

        //- No copy construct
        referenceTemperature(const referenceTemperature&) = delete;

        //- No copy assignment
        void operator=(const referenceTemperature&) = delete;


    //- Destructor
    virtual ~referenceTemperature() = default;


    // Member Functions

    // Evaluation

        //- Initialise data members of the model
        virtual void initialise();

        //- Return energy density per unit length [J/m3/m]
        virtual tmp<scalarField> energyDensity(const labelList& cells);


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Write data to stream and files
        virtual void write(const bool log);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatExchangerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
