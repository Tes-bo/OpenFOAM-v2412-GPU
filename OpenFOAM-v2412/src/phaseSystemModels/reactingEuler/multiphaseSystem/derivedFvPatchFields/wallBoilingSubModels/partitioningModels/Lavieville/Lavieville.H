/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::wallBoilingModels::partitioningModels::Lavieville

Description
    A model for wall heat flux partitioning based on
    Lavieville et al. (2006) for boiling flows.

    Reference:
    \verbatim
        Lavieville, J., Quemerais, E., Mimouni, S.,
        Boucker, M., & Mechitoua, N. (2006).
        NEPTUNE CFD V1.0 theory manual.
        NEPTUNE report Nept_2004_L1, 2(3).
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    partitioningModel
    {
        // Mandatory entries
        type            Lavieville;
        alphaCrit       <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                       | Type | Reqd | Deflt
      type          | Type name: Lavieville             | word | yes  | -
      alphaCrit     | Critical liquid fraction        | scalar | yes  | -
    \endtable

SourceFiles
    Lavieville.C

\*---------------------------------------------------------------------------*/

#ifndef Lavieville_H
#define Lavieville_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class Lavieville Declaration
\*---------------------------------------------------------------------------*/

class Lavieville
:
    public partitioningModel
{
    // Private data

        //- Critical liquid fraction
        scalar alphaCrit_;


public:

    //- Runtime type information
    TypeName("Lavieville");


    // Constructors

        //- Construct from a dictionary
        Lavieville(const dictionary& dict);


    //- Destructor
    virtual ~Lavieville();


    // Member Functions

        //- Liquid blending function
        virtual tmp<scalarField> fLiquid(const scalarField& alphaLiquid) const;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
