/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

Namespace
    Foam::filmSeparationModels

Description
    A namespace for various \c filmSeparation model implementations.

Class
    Foam::filmSeparationModel

Description
    A base class for \c filmSeparation models.

SourceFiles
    filmSeparationModel.C
    filmSeparationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparationModel_H
#define Foam_filmSeparationModel_H

#include "liquidFilmBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class filmSeparationModel Declaration
\*---------------------------------------------------------------------------*/

class filmSeparationModel
{
    // Private Data

        //- Const reference to the film properties
        const regionModels::areaSurfaceFilmModels::liquidFilmBase& film_;


public:

    //- Runtime type information
    TypeName("filmSeparationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            filmSeparationModel,
            dictionary,
            (
                const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Selectors

        //- Return a reference to the selected filmSeparation model
        static autoPtr<filmSeparationModel> New
        (
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        filmSeparationModel(const filmSeparationModel&) = delete;

        //- No copy assignment
        void operator=(const filmSeparationModel&) = delete;


    // Constructors

        //- Construct from the base film model and dictionary
        filmSeparationModel
        (
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmSeparationModel() = default;


    // Member Functions

    // Access

        //- Return const access to the film properties
        const regionModels::areaSurfaceFilmModels::liquidFilmBase& film() const
        {
            return film_;
        }

        //- Return const access to the finite-area mesh
        const faMesh& mesh() const noexcept { return film_.regionMesh(); }


    // Evaluation

        //- Calculate the mass ratio of film separation
        virtual tmp<scalarField> separatedMassRatio() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
