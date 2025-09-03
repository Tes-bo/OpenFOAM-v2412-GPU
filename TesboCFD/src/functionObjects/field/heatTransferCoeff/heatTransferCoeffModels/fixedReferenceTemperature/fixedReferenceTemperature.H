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
    Foam::heatTransferCoeffModels::fixedReferenceTemperature

Description
    Heat transfer coefficient calculation
    that employs a fixed reference temperature.

    The heat transfer coefficient is calculated by:

    \f[
        h = \frac{q}{T_{ref} - T_p}
    \f]

    where
    \vartable
      h       | Heat transfer coefficient                         [W/m^2/K]
      q       | Heat flux                                         [W/m^2]
      T_{ref} | User-specified reference temperature              [K]
      T_p     | Patch temperature field                           [K]
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Inherited entries
        ...

        // Mandatory entries
        htcModel        fixedReferenceTemperature;
        Tref            <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type      | Model name: fixedReferenceTemperature  | word | yes  | -
      Tref      | User-specied reference temperature   | scalar | yes  | -
    \endtable

SourceFiles
    fixedReferenceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_fixedReferenceTemperature_H
#define heatTransferCoeffModels_fixedReferenceTemperature_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                   Class fixedReferenceTemperature Declaration
\*---------------------------------------------------------------------------*/

class fixedReferenceTemperature
:
    public heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Reference temperature [K]
        scalar TRef_;


    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );


public:

    //- Runtime type information
    TypeName("fixedReferenceTemperature");


    // Constructors

        //- Construct from components
        fixedReferenceTemperature
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        fixedReferenceTemperature(const fixedReferenceTemperature&) = delete;

        //- No copy assignment
        void operator=(const fixedReferenceTemperature&) = delete;


    //- Destructor
    virtual ~fixedReferenceTemperature() = default;


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
