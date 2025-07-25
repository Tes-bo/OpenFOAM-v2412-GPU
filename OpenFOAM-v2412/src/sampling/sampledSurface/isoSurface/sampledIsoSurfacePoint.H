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
    Foam::sampledIsoSurfacePoint

Description
    A sampledSurface defined by a surface of iso value using a
    \em point algorithm (always triangulated!).
    It only recalculates the iso-surface if time changes.
    To be used in sampleSurfaces / functionObjects.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type     isoSurfacePoint;
            isoField T;
            isoValue 373;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | isoSurfacePoint                          | yes |
        isoField | field name for obtaining iso-surface     | yes |
        isoValue | value of iso-surface                     | yes |
        average  | cell values from averaged point values   | no  | false
        bounds   | limit with bounding box                  | no  |
        zone     | limit to cell zone (name or regex)       | no  |
        zones    | limit to cell zones (names, regexs)      | no  |
        regularise | point snapping                         | yes |
        mergeTol | tolerance for merging points             | no  | 1e-6
    \endtable

Note
    For the isoMethod \b point should use a "cellPoint" sampleScheme
    since the triangles generated with this algorithm span across cells.

SourceFiles
    sampledIsoSurfacePoint.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurfacePoint_H
#define sampledIsoSurfacePoint_H

#include "sampledIsoSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledIsoSurfacePoint Declaration
\*---------------------------------------------------------------------------*/

class sampledIsoSurfacePoint
:
    public sampledIsoSurface
{
public:

    //- Runtime type information
    TypeNameNoDebug("sampledIsoSurfacePoint");


    // Constructors

        //- Construct from dictionary
        sampledIsoSurfacePoint
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledIsoSurfacePoint() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
