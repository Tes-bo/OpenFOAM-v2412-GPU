/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::standardPhaseChange

Description
    Standard phase change model with modification for boiling

SourceFiles
    standardPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_standardPhaseChange_H
#define Foam_standardPhaseChange_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class standardPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class standardPhaseChange
:
    public phaseChangeModel
{
protected:

    // Protected Data

        //- Minimum film height for model to be active
        const scalar deltaMin_;

        //- Length scale / [m]
        const scalar L_;

        //- Boiling temperature factor / []
        //  Used to set max limit on temperature to Tb*TbFactor
        const scalar TbFactor_;

        //- Switch to treat YInf as zero
        bool YInfZero_;


    // Protected member functions

        //- Return Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;

        template<class YInfType>
        void correctModel
        (
            const scalar dt,
            scalarField& availableMass,
            scalarField& dMass,
            scalarField& dEnergy,
            const YInfType& YInf
        );

        //- No copy construct
        standardPhaseChange(const standardPhaseChange&) = delete;

        //- No copy assignment
        void operator=(const standardPhaseChange&) = delete;


public:

    //- Runtime type information
    TypeName("standardPhaseChange");


    // Constructors

        //- Construct from surface film model
        standardPhaseChange
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~standardPhaseChange() = default;


    // Member Functions

        // Evolution

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
