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
    Foam::CHFModels::CHFSubCoolModel

Description
    Base class for critical heat flux (CHF)
    sub-cooling correlation models for boiling flows.

SourceFiles
    CHFSubCoolModel.C

\*---------------------------------------------------------------------------*/

#ifndef CHFSubCoolModel_H
#define CHFSubCoolModel_H

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
                         Class CHFSubCoolModel Declaration
\*---------------------------------------------------------------------------*/

class CHFSubCoolModel
{
public:

    //- Runtime type information
    TypeName("CHFSubCoolModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        CHFSubCoolModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Generated Methods

        //- Default construct
        CHFSubCoolModel() = default;

        //- Destructor
        virtual ~CHFSubCoolModel() = default;


    // Selectors

        //- Select default constructed
        static autoPtr<CHFSubCoolModel> New(const dictionary& dict);


    // Member Functions

        //- Calculate temperature
        virtual tmp<scalarField> CHFSubCool
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
