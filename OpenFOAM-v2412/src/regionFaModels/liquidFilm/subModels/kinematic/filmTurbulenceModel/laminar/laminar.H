/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::laminar

Description
    Film laminar turbulence model.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_surfaceFilmModels_laminar_H
#define regionFaModels_surfaceFilmModels_laminar_H

#include "filmTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                           Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public filmTurbulenceModel
{
    // Private Member Functions

        //- No copy construct
        laminar(const laminar&) = delete;

        //- No copy assignment
        void operator=(const laminar&) = delete;


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from surface film model
        laminar(liquidFilmBase& film, const dictionary& dict);


    //- Destructor
    virtual ~laminar() = default;


    // Member Functions

        // Turbulence

            //- Wall friction
            tmp<faVectorMatrix> wallFriction(areaVectorField& U) const;

            //- Return the film turbulence viscosity
            virtual tmp<areaScalarField> mut() const;

            //- Return the source for the film momentum equation
            virtual tmp<faVectorMatrix> Su(areaVectorField& U) const;


        // Evaluation

            //- Correct/update the model
            virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
