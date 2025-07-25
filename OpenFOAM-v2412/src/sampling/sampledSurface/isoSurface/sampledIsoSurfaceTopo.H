/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::sampledIsoSurfaceTopo

Description
    A sampledSurface defined by a surface of iso value using a
    \em topo algorithm.
    It only recalculates the iso-surface if time changes.
    To be used in sampleSurfaces / functionObjects.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type     isoSurfaceTopo;
            isoField T;
            isoValue 373;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | isoSurfaceTopo                           | yes |
        isoField | field name for obtaining iso-surface     | yes |
        isoValue | value of iso-surface                     | yes |
        average  | cell values from averaged point values   | no  | false
        bounds   | limit with bounding box                  | no  |
        zone     | limit to cell zone (name or regex)       | no  |
        zones    | limit to cell zones (names, regexs)      | no  |
        regularise | filter faces (bool or enum)            | no  | true
        triangulate | triangulate faces (requires regularise) | no  | false
    \endtable

SourceFiles
    sampledIsoSurfaceTopo.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurfaceTopo_H
#define sampledIsoSurfaceTopo_H

#include "sampledIsoSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledIsoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class sampledIsoSurfaceTopo
:
    public sampledIsoSurface
{
public:

    //- Runtime type information
    TypeNameNoDebug("sampledIsoSurfaceTopo");


    // Constructors

        //- Construct from dictionary
        sampledIsoSurfaceTopo
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledIsoSurfaceTopo() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
