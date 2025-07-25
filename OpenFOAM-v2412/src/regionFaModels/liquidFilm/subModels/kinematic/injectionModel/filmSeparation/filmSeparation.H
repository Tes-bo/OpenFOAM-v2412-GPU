/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::filmSeparation

Description
    The \c filmSeparation is a collection of curvature thin-film separation
    models designed to predict the onset of film separation and mass separation
    in geometries featuring sharp and/or rounded corners.

Usage
    Minimal example by using boundary-condition files:
    \verbatim
    injectionModels
    {
        // Mandatory entries
        filmSeparation
    }

    filmSeparationCoeffs
    {
        model           <word>;

        // Conditional entries

            // Option-1: when model == OwenRyley

            // Option-2: when model == Friedrich

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      model        | Name of the filmSeparation model   | word | yes  | -
    \endtable

    Options for the \c model entry:
    \verbatim
      OwenRyley    | Model proposed by Owen-Ryley (1985)
      Friedrich    | Model proposed by Friedrich et al. (2008)
    \endverbatim

    The inherited entries are elaborated in:
    - \link injectionModel.H \endlink

SourceFiles
    filmSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparation_H
#define Foam_filmSeparation_H

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class filmSeparationModel;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class filmSeparation Declaration
\*---------------------------------------------------------------------------*/

class filmSeparation
:
    public injectionModel
{
    // Private Data

        //- Film-separation model
        autoPtr<filmSeparationModel> filmSeparationModelPtr_;


public:

    //- Runtime type information
    TypeName("filmSeparation");


    // Constructors

        //- Construct from base film model and dictionary
        filmSeparation
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmSeparation();


    // Member Functions

        //- Correct film properties due to the film separation
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
