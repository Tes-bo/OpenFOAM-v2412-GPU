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
    Foam::wallBoilingModels::nucleationSiteModel

Description
    Base class for nucleation site density models.

SourceFiles
    nucleationSiteModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationSiteModel_H
#define nucleationSiteModel_H

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
                         Class nucleationSiteModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationSiteModel
{
public:

    //- Runtime type information
    TypeName("nucleationSiteModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationSiteModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Default construct
        nucleationSiteModel() = default;


    // Selectors

        //- Select default constructed
        static autoPtr<nucleationSiteModel> New(const dictionary& dict);


    //- Destructor
    virtual ~nucleationSiteModel() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> N
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
