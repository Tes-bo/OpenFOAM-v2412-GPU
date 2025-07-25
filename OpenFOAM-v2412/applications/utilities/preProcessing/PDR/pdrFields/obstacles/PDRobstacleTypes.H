/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Namespace
    Foam::PDRobstacles

Description
    Reader classes for concrete PDRsetFields obstacle types.

SourceFiles
    PDRobstacleTypes.C

\*---------------------------------------------------------------------------*/

#ifndef PDRobstacleTypes_H
#define PDRobstacleTypes_H

#include "PDRobstacle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRobstacles
{

/*---------------------------------------------------------------------------*\
                          Class cylinder Declaration
\*---------------------------------------------------------------------------*/

/*! \brief A cylinder, selectable as \c cyl or \c cylinder

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The end centre-point                      | yes |
    length      | The cylinder length                       | yes |
    diameter    | The cylinder diameter                     | yes |
    direction   | The x/y/z direction                       | yes |
\endtable

Example,
\verbatim
cyl
{
     point       (0 0 0);
     length      0.95;
     diameter    0.05;
     direction   x;
}
\endverbatim

\note Does not have porosity.
*/
struct cylinder : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::CYLINDER;
    static void read(PDRobstacle& obs, const dictionary& dict);
};



/*---------------------------------------------------------------------------*\
                          Class diagbeam Declaration
\*---------------------------------------------------------------------------*/

/*! \brief A diagonal beam, which is cylinder-like,
 selectable as \c diag or \c diagbeam

 If the beam angle is close to zero, it will be changed to a box
 (PDRobstacles::cuboid)

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The end centre-point                      | yes |
    length      | The beam length                           | yes |
    width       | The beam cross-dimensions                 | yes |
    angle       | The beam angle (degrees)                  | yes |
    direction   | The x/y/z direction                       | yes |
\endtable

Example,
\verbatim
diag
{
     point       (0 0 0);
     length      0.95;
     width       (0.05 0.01);
     angle       25;
     direction   x;
}
\endverbatim

\note Does not have porosity.
*/
struct diagbeam : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::DIAG_BEAM;
    static void read(PDRobstacle& obs, const dictionary& dict);
};


/*---------------------------------------------------------------------------*\
                           Class cuboid Declaration
\*---------------------------------------------------------------------------*/

/*! \brief A cuboid, selectable as \c box

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The lower left corner                     | yes |
    size        | The (x y z) dimensions                    | yes |
    porosity    | The volumetric porosity                   | no  | 0
    porosities  | The directional porosities                | no  | (0 0 0)
\endtable

Example,
\verbatim
box
{
     point       (0 0 0);
     size        (0.05 0.05 2);
}
\endverbatim
*/
struct cuboid : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::CUBOID;
    static void read(PDRobstacle& obs, const dictionary& dict);
};



/*---------------------------------------------------------------------------*\
                           Class wallbeam Declaration
\*---------------------------------------------------------------------------*/

/*! \brief A wallbeam, selectable as \c wallbeam which is currently identical
 to a box (PDRobstacles::cuboid)

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The lower left corner                     | yes |
    size        | The (x y z) dimensions                    | yes |
\endtable

Example,
\verbatim
wallbeam
{
     point       (0 0 0);
     size        (0.05 0.05 2);
}
\endverbatim

\note Does not have porosity.
*/
struct wallbeam : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::WALL_BEAM;
    static void read(PDRobstacle& obs, const dictionary& dict);
};


/*---------------------------------------------------------------------------*\
                           Class grating Declaration
\*---------------------------------------------------------------------------*/

/*! \brief A grating, selectable as \c grate or \c grating

The dimensions must include one component that is zero,
which establishes the orientation.

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The lower left corner                     | yes |
    size        | The (x y z) dimensions                    | yes |
    slats       | The slat width                            | no  | 0
\endtable

Example,
\verbatim
grating
{
     point       (0 0 0);
     size        (0.1 0.1 0);
     slats       0.005;
}
\endverbatim
*/
struct grating : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::GRATING;
    static void read(PDRobstacle& obs, const dictionary& dict);
};


/*---------------------------------------------------------------------------*\
                           Class louver Declaration
\*---------------------------------------------------------------------------*/

/*! \brief Louver blowoff, selectable as \c louver or \c louvre

Dictionary controls
\table
    Property    | Description                               | Required | Default
    point       | The lower left corner                     | yes |
    size        | The (x y z) dimensions                    | yes |
    pressure    | The blowoff pressure (bar)                | yes |
    time        | The blowoff time                          | no  | 0
    type        | The blowoff type (1, 2)                   | no  | 1
\endtable

Example,
\verbatim
louver
{
     point       (0 0 0);
     size        (0.1 0.1);
     pressure    3;
     type        1;
}
\endverbatim
*/
struct louver : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::LOUVER_BLOWOFF;
    static void read(PDRobstacle& obs, const dictionary& dict);
};


/*---------------------------------------------------------------------------*\
                           Class patch Declaration
\*---------------------------------------------------------------------------*/

/*! \brief Rectangular patch, selectable as \c patch

Dictionary controls
\table
    Property    | Description                               | Required | Default
    name        | The patch name corner                     | yes |
    point       | The lower left corner                     | yes |
    size        | The (x y z) dimensions                    | yes |
    direction   | The direction (x/y/z, _x/_y/_z or "-x"/"-y"/"-z" | yes |
\endtable

Example,
\verbatim
patch
{
     name        inlet;
     direction   _x;
     point       (0 0 0);
     size        (0 0.05 1.0);
}
\endverbatim
*/
struct patch : public PDRobstacle
{
    static constexpr int enumTypeId = legacyTypes::RECT_PATCH;
    static void read(PDRobstacle& obs, const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDRobstacles
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
