/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::extrudeModels::linearRadial

Description
    Extrudes by transforming radially from the surface.

    The linearRadialCoeffs dictionary entries
    \table
        Property  | Description                             | Required | Default
        R         | The radial thickness                    | yes |
        Rsurface  | Surface inner radius                    | no  | -1
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef linearRadial_H
#define linearRadial_H

#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                 Class extrudeModels::linearRadial Declaration
\*---------------------------------------------------------------------------*/

class linearRadial
:
    public extrudeModel
{
    // Private data

        scalar R_;
        scalar Rsurface_;


public:

    //- Runtime type information
    TypeName("linearRadial");

    // Constructors

        //- Construct from dictionary
        explicit linearRadial(const dictionary& dict);


    //- Destructor
    virtual ~linearRadial() = default;


    // Member Operators

        //- Return point
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
