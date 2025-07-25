/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::LeidenfrostModel

Description
    Base class for Leidenfrost-effect models.

SourceFiles
    LeidenfrostModel.C

\*---------------------------------------------------------------------------*/

#ifndef LeidenfrostModel_H
#define LeidenfrostModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class LeidenfrostModel Declaration
\*---------------------------------------------------------------------------*/

class LeidenfrostModel
{
public:

    //- Runtime type information
    TypeName("LeidenfrostModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        LeidenfrostModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Generated Methods

        //- Default construct
        LeidenfrostModel() = default;

        //- Destructor
        virtual ~LeidenfrostModel() = default;


    // Selectors

        //- Select default constructed
        static autoPtr<LeidenfrostModel> New(const dictionary& dict);


    // Member Functions

        //- Calculate and return the Leidenfrost temperature
        virtual tmp<scalarField> TLeid
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const = 0;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
