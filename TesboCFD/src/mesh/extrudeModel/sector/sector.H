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
    Foam::extrudeModels::sector

Description
    Extrudes by rotating a surface around an axis
    - extrusion is opposite the surface/patch normal so inwards the source
      mesh
    - axis direction has to be consistent with this.
    - use -mergeFaces option if doing full 360 and want to merge front and back
    - note direction of axis. This should be consistent with rotating against
      the patch normal direction. If you get it wrong you'll see all cells
      with extreme aspect ratio and internal faces wrong way around in
      checkMesh

    The sectorCoeffs dictionary entries
    \table
        Property  | Description                             | Required | Default
        point     | A point on the axis                     | yes |
        axis      | Axis direction, normalized on input     | yes |
        angle     | Sector angle (degrees)                  | yes |
    \endtable

Note
    For compatibility, accepts the entry \c axisPt (1812 and earlier)
    as equivalent to \c point.

\*---------------------------------------------------------------------------*/

#ifndef sector_H
#define sector_H

#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                    Class extrudeModels::sector Declaration
\*---------------------------------------------------------------------------*/

class sector
:
    public extrudeModel
{
    // Private data

        //- Point on axis
        const point refPoint_;

        //- Normalized axis direction
        const vector axis_;

        //- Overall angle (radians), converted from degrees on input
        const scalar angle_;


public:

    //- Runtime type information
    TypeName("sector");

    // Constructors

        //- Construct from dictionary
        explicit sector(const dictionary& dict);


    //- Destructor
    virtual ~sector() = default;


    // Member Operators

        //- Return corresponding point
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
