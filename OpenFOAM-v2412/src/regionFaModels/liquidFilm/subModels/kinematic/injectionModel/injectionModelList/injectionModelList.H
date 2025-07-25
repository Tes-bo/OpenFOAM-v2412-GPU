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
    Foam::regionModels::surfaceFilmModels::injectionModelList

Description
    List container for film injection models

Usage
    Example of the model specification:
    \verbatim
    {
        // Mandatory entries
        injectionModels
        (
            <injectionModel1>
            <injectionModel2>
            ...
        );

        <injectionModel1>Coeffs
        {
            ...
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property        | Description                   | Type  | Reqd | Deflt
      injectionModels | List of injection models      | list  | yes  | -
    \endtable

SourceFiles
    injectionModelList.C

\*---------------------------------------------------------------------------*/

#ifndef injectionModelList_H
#define injectionModelList_H

#include "injectionModel.H"
#include "filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class injectionModelList Declaration
\*---------------------------------------------------------------------------*/

class injectionModelList
:
    public PtrList<injectionModel>,
    public filmSubModelBase
{
    // Private Data

        //- List of mass injected per patch
        scalar massInjected_;


    // Private Member Functions

        //- No copy construct
        injectionModelList(const injectionModelList&) = delete;

        //- No copy assignment
        void operator=(const injectionModelList&) = delete;


public:

    // Constructors

        //- Construct null
        injectionModelList(liquidFilmBase& film);

        //- Construct from type name, dictionary and surface film model
        injectionModelList
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~injectionModelList();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                volScalarField& massToInject,
                volScalarField& diameterToInject
            );


        // I-O

            //- Provide some info
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
