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
    Foam::regionModels::surfaceFilmModels::forceList

Description
    List container for film sources

Usage
    Example of the model specification:
    \verbatim
    {
        // Mandatory entries
        forces
        (
            <forceModel1>
            <forceModel2>
            ...
        );

        <forceModel1>Coeffs
        {
            ...
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                   | Type  | Reqd | Deflt
      forces       | List of force models          | list  | yes  | -
    \endtable

SourceFiles
    forceList.C

\*---------------------------------------------------------------------------*/

#ifndef forceList_H
#define forceList_H

#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class forceList Declaration
\*---------------------------------------------------------------------------*/

class forceList
:
    public PtrList<force>
{
public:

    // Constructors

        //- Construct null
        forceList(liquidFilmBase& film);

        //- Construct from type name, dictionary and surface film model
        forceList
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~forceList();


    // Member Functions

        //- Return (net) force system
        tmp<faVectorMatrix> correct(areaVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
