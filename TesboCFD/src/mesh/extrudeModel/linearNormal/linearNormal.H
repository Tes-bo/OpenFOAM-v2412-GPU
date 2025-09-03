/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::extrudeModels::linearNormal

Description
    Extrudes by transforming points normal to the surface by a given distance.

    The linearNormalCoeffs dictionary entries
    \table
        Property  | Description                             | Required | Default
        thickness | Extrusion thickness                     | yes |
        firstCellThickness | Thickness of the first layer   | no  | 0
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef linearNormal_H
#define linearNormal_H

#include "point.H"
#include "extrudeModel.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                 Class extrudeModels::linearNormal Declaration
\*---------------------------------------------------------------------------*/

class linearNormal
:
    public extrudeModel
{
    // Private data

        //- Layer thickness
        scalar thickness_;

        //- First cell thickness
        scalar firstCellThickness_;

        //- Layer cell distribution
        scalarList layerPoints_;


public:

    //- Runtime type information
    TypeName("linearNormal");

    // Constructors

        //- Construct from dictionary
        explicit linearNormal(const dictionary& dict);


    //- Destructor
    virtual ~linearNormal() = default;


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
