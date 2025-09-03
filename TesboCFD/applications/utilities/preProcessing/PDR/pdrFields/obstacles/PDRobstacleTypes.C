/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "PDRobstacleTypes.H"
#include "PDRobstacleTypes.H"
#include "Enum.H"
#include "unitConversion.H"
#include "addToMemberFunctionSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addObstacleReader(obsType, obsName)                                   \
    namespace Foam                                                            \
    {                                                                         \
    namespace PDRobstacles                                                    \
    {                                                                         \
        addNamedToMemberFunctionSelectionTable                                \
        (                                                                     \
            PDRobstacle,                                                      \
            obsType,                                                          \
            read,                                                             \
            dictionary,                                                       \
            obsName                                                           \
        );                                                                    \
    }                                                                         \
    }


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read porosity, change to blockage. Clamp values [0-1] silently

// Volume porosity -> blockage
inline scalar getPorosity(const dictionary& dict)
{
    return 1 - clamp(dict.getOrDefault<scalar>("porosity", 0), zero_one{});
}

// Direction porosities -> blockage
inline vector getPorosities(const dictionary& dict)
{
    vector blockage(vector::one);

    if (dict.readIfPresent("porosities", blockage))
    {
        for (scalar& val : blockage)
        {
            val = 1 - clamp(val, zero_one{});
        }
    }

    return blockage;
}


// Check for "porosity", or "porosities"
// inline static bool hasPorosity(const dictionary& dict)
// {
//     return dict.found("porosity") || dict.found("porosities");
// }


static const Foam::Enum<Foam::vector::components>
vectorComponentsNames
({
    { vector::components::X, "x" },
    { vector::components::Y, "y" },
    { vector::components::Z, "z" },
});


enum inletDirnType
{
    _X = -1,    // -ve x
    _Y = -2,    // -ve y
    _Z = -3,    // -ve z
    X = 1,      // +ve x
    Y = 2,      // +ve y
    Z = 3,      // +ve z
};

static const Foam::Enum<inletDirnType>
inletDirnNames
({
    { inletDirnType::_X, "-x" },
    { inletDirnType::_Y, "-y" },
    { inletDirnType::_Z, "-z" },
    { inletDirnType::_X, "_x" },
    { inletDirnType::_Y, "_y" },
    { inletDirnType::_Z, "_z" },
    { inletDirnType::X, "+x" },
    { inletDirnType::Y, "+y" },
    { inletDirnType::Z, "+z" },
    { inletDirnType::X, "x" },
    { inletDirnType::Y, "y" },
    { inletDirnType::Z, "z" },
});

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(cylinder, cyl);
addObstacleReader(cylinder, cylinder);

void Foam::PDRobstacles::cylinder::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    // Enforce complete blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;
    // if (hasPorosity(dict)) ... warn?


    dict.readEntry("point", obs.pt);
    dict.readEntry("length", obs.len());
    dict.readEntry("diameter", obs.dia());

    obs.orient = vectorComponentsNames.get("direction", dict);

    // The sortBias for later position sorting
    switch (obs.orient)
    {
        case vector::X:
            obs.sortBias = obs.len();
            break;

        default:
            obs.sortBias = 0.5*obs.dia();
            break;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(diagbeam, diag);
addObstacleReader(diagbeam, diagbeam);

void Foam::PDRobstacles::diagbeam::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    // Enforce complete blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;
    // if (hasPorosity(dict)) ... warn?


    dict.readEntry("point", obs.pt);
    dict.readEntry("length", obs.len());
    obs.dia() = Zero;
    obs.theta() = Zero;   // Fix later on

    obs.orient = vectorComponentsNames.get("direction", dict);

    // Angle (degrees) on input, limit to range [0, PI]
    scalar angle = dict.get<scalar>("angle");

    while (angle > 180) angle -= 180;
    while (angle < 0) angle += 180;

    labelPair dims;
    dict.readEntry("width", dims);

    // Swap axes when theta > PI/2
    // For 89-90 degrees it becomes -ve, which is picked up in following section
    if (angle > 89)
    {
        angle -= 90;
        dims.flip();  // Swap dimensions
    }

    obs.theta() = degToRad(angle);

    obs.wa = dims.first();
    obs.wb = dims.second();

    const scalar ctheta = cos(obs.theta());
    const scalar stheta = sin(obs.theta());

    // The sortBias for later position sorting
    switch (obs.orient)
    {
        case vector::X:
            obs.sortBias = obs.len();
            break;

        case vector::Y:
            obs.sortBias = 0.5*(obs.wa * stheta + obs.wb * ctheta);
            break;

        case vector::Z:
            obs.sortBias = 0.5*(obs.wa * ctheta + obs.wb * stheta);
            break;
    }

    // If very nearly aligned with axis, turn it into normal block,
    // to avoid 1/tan(theta) blowing up
    if (angle < 1)
    {
        Info<< "... changed diag-beam to box" << nl;

        switch (obs.orient)
        {
            case vector::X:
                obs.span = vector(obs.len(), obs.wa, obs.wb);
                break;

            case vector::Y:
                obs.span = vector(obs.wb, obs.len(), obs.wa);
                break;

            case vector::Z:
                obs.span = vector(obs.wa, obs.wb, obs.len());
                break;
        }

        // The pt was end centre (cylinder), now lower corner
        vector adjustPt = -0.5*obs.span;
        adjustPt[obs.orient] = 0;

        obs.pt -= adjustPt;

        obs.typeId = PDRobstacles::cuboid::enumTypeId;
        obs.sortBias = 0;
        obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1.0;
        obs.blowoff_type = 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(cuboid, box);

void Foam::PDRobstacles::cuboid::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    // Default - full blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;


    dict.readEntry("point", obs.pt);
    dict.readEntry("size", obs.span);

    // Optional
    obs.vbkge = getPorosity(dict);

    // Optional
    const vector blockages = getPorosities(dict);
    obs.xbkge = blockages.x();
    obs.ybkge = blockages.y();
    obs.zbkge = blockages.z();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(wallbeam, wallbeam);

void Foam::PDRobstacles::wallbeam::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    PDRobstacles::cuboid::read(obs, dict);
    obs.typeId = enumTypeId;

    // Enforce complete blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;

    // if (hasPorosity(dict)) ... warn?

    // Additional adjustment for drag etc.
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(grating, grating);
addObstacleReader(grating, grate);

void Foam::PDRobstacles::grating::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    // Initialize to full blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;

    dict.readEntry("point", obs.pt);
    dict.readEntry("size", obs.span);

    // TODO: better error handling
    // if (!equal(cmptProduct(obs.span), 0))
    // {
    //     Info<< "Type " << typeId << " has non-zero thickness.";
    //     ReportLineInfo(lineNo, inputFile);
    // }

    obs.vbkge = getPorosity(dict);

    const vector blockages = getPorosities(dict);
    obs.xbkge = blockages.x();
    obs.ybkge = blockages.y();
    obs.zbkge = blockages.z();

    // TODO: Warning if porosity was not specified?


    // TBD: Default slat width from PDR.params
    obs.slat_width = dict.getOrDefault<scalar>("slats", Zero);

    // Determine the orientation
    if (equal(obs.span.x(), 0))
    {
        obs.orient = vector::X;
    }
    else if (equal(obs.span.y(), 0))
    {
        obs.orient = vector::Y;
    }
    else
    {
        obs.orient = vector::Z;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(louver, louver);
addObstacleReader(louver, louvre);

void Foam::PDRobstacles::louver::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    // Initialize to full blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;

    dict.readEntry("point", obs.pt);
    dict.readEntry("size", obs.span);

    // TODO: better error handling
    // if (!equal(cmptProduct(obs.span), 0))
    // {
    //     Info<< "Type " << typeId << " has non-zero thickness.";
    //     ReportLineInfo(lineNo, inputFile);
    // }

    obs.vbkge = getPorosity(dict);

    const vector blockages = getPorosities(dict);
    obs.xbkge = blockages.x();
    obs.ybkge = blockages.y();
    obs.zbkge = blockages.z();

    // TODO: Warning if porosity was not specified?


    // Blowoff pressure [bar]
    const scalar blowoffPress = dict.get<scalar>("pressure");

    obs.blowoff_press = barToPa(blowoffPress);
    obs.blowoff_time = dict.getOrDefault<scalar>("time", 0);
    obs.blowoff_type = dict.getOrDefault<label>("type", 2);

    if (obs.blowoff_type == 1)
    {
        Info<< "Louver : blowoff-type 1 not yet implemented." << nl;
        // ReportLineInfo(lineNo, inputFile);

        if (obs.blowoff_time != 0)
        {
            Info<< "Louver : has blowoff time set,"
                << " not set to blow off cell-by-cell" << nl;
            // ReportLineInfo(lineNo, inputFile);
        }
    }
    else
    {
        if
        (
            (obs.blowoff_type == 1 || obs.blowoff_type == 2)
         && (blowoffPress > 0)
        )
        {
            if (blowoffPress > maxBlowoffPressure)
            {
                Info<< "Blowoff pressure (" << blowoffPress
                    << ") too high for blowoff type "
                    << obs.blowoff_type << nl;
                // ReportLineInfo(lineNo, inputFile);
            }
        }
        else
        {
            Info<< "Problem with blowoff parameters" << nl;
            // ReportLineInfo(lineNo, inputFile);
            Info<< "Pressure[bar] " << blowoffPress
                << " Blowoff type " << obs.blowoff_type
                << ", blowoff pressure " << obs.blowoff_press << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

addObstacleReader(patch, patch);

void Foam::PDRobstacles::patch::read
(
    PDRobstacle& obs,
    const dictionary& dict
)
{
    obs.PDRobstacle::readProperties(dict);
    obs.typeId = enumTypeId;

    const auto nameLen = obs.identifier.length();

    word patchName = word::validate(obs.identifier);

    if (patchName.empty())
    {
        FatalErrorInFunction
            << "RECT_PATCH without a patch name"
            << exit(FatalError);
    }
    else if (patchName.length() != nameLen)
    {
        WarningInFunction
            << "RECT_PATCH stripped invalid characters from patch name: "
            << obs.identifier
            << exit(FatalError);

        obs.identifier = std::move(patchName);
    }

    // Enforce complete blockage
    obs.xbkge = obs.ybkge = obs.zbkge = obs.vbkge = 1;

    dict.readEntry("point", obs.pt);
    dict.readEntry("size", obs.span);
    obs.inlet_dirn = inletDirnNames.get("direction", dict);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef addObstacleReader

// ************************************************************************* //
