/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::regionFaModels::kinematicThinFilm

Description
    Thin film model.

Usage
    Example of the boundary condition specification:
    \verbatim
    {
        // Mandatory entries
        liquidFilmModel         kinematicThinFilm;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type  | Reqd | Deflt
      liquidFilmModel  | Type name: kinematicThinFilm | dict  | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link liquidFilmModel.H \endlink

SourceFiles
    kinematicThinFilm.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicThinFilm_H
#define kinematicThinFilm_H

#include "volFieldsFwd.H"
#include "liquidFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class kinematicThinFilm Declaration
\*---------------------------------------------------------------------------*/

class kinematicThinFilm
:
    public liquidFilmModel
{
public:

    //- Runtime type information
    TypeName("kinematicThinFilm");


    // Constructors

        //- Construct from components and dict
        kinematicThinFilm
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        kinematicThinFilm(const kinematicThinFilm&) = delete;

        //- No copy assignment
        void operator=(const kinematicThinFilm&) = delete;


    //- Destructor
    virtual ~kinematicThinFilm() = default;


    // Member Functions

        // Evolution

            //- Pre-evolve film
            virtual void preEvolveRegion();

            //- Evolve the film
            virtual void evolveRegion();

            //- Post-evolve film
            virtual void postEvolveRegion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
