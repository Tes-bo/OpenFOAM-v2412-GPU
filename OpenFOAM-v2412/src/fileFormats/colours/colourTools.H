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
    Foam::colourTools

Description
    Utility methods for colours and colour spaces

SourceFiles
    colourTools.C

\*---------------------------------------------------------------------------*/

#ifndef colourTools_H
#define colourTools_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace colourTools
{

/*---------------------------------------------------------------------------*\
                             Namespace colourTools
\*---------------------------------------------------------------------------*/

//- Convert RGB to HSV
void rgbToHsv(const vector& rgb, vector& hsv);

//- Convert RGB to HSV
inline vector rgbToHsv(const vector& rgb)
{
    vector hsv;
    rgbToHsv(rgb, hsv);
    return hsv;
}


//- Convert HSV to RGB
void hsvToRgb(const vector& hsv, vector& rgb);

//- Convert HSV to RGB
inline vector hsvToRgb(const vector& hsv)
{
    vector rgb;
    hsvToRgb(hsv, rgb);
    return rgb;
}


//- Convert RGB to XYZ
void rgbToXyz(const vector& rgb, vector& xyz);

//- Convert RGB to XYZ
inline vector rgbToXyz(const vector& rgb)
{
    vector xyz;
    rgbToXyz(rgb, xyz);
    return xyz;
}

//- Convert XYZ to RGB
void xyzToRgb(const vector& xyz, vector& rgb);

//- Convert XYZ to RGB
inline vector xyzToRgb(const vector& xyz)
{
    vector rgb;
    xyzToRgb(xyz, rgb);
    return rgb;
}


//- Convert LAB to XYZ
void labToXyz(const vector& lab, vector& xyz);

//- Convert LAB to XYZ
inline vector labToXyz(const vector& lab)
{
    vector xyz;
    labToXyz(lab, xyz);
    return xyz;
}


//- Convert XYZ to LAB
void xyzToLab(const vector& xyz, vector& lab);

//- Convert XYZ to LAB
inline vector xyzToLab(const vector& xyz)
{
    vector lab;
    xyzToLab(xyz, lab);
    return lab;
}


//- Convert RGB to LAB
void rgbToLab(const vector& rgb, vector& lab);

//- Convert RGB to LAB
inline vector rgbToLab(const vector& rgb)
{
    vector lab;
    rgbToLab(rgb, lab);
    return lab;
}


//- Convert LAB to RGB
void labToRgb(const vector& lab, vector& rgb);

//- Convert LAB to RGB
inline vector labToRgb(const vector& lab)
{
    vector rgb;
    labToRgb(lab, rgb);
    return rgb;
}


//- Interpolate RGB values with diverging color map
void interpolateDiverging
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2,
    vector& result
);

//- Interpolate RGB values with diverging color map
inline vector interpolateDiverging
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2
)
{
    vector result;
    interpolateDiverging(s, rgb1, rgb2, result);
    return result;
}


//- Interpolate RGB values in HSV colourspace
void interpolateHSV
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2,
    vector& result
);

//- Interpolate RGB values in HSV colourspace
inline vector interpolateHSV
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2
)
{
    vector result;
    interpolateHSV(s, rgb1, rgb2, result);
    return result;
}


//- Interpolate RGB values in RGB colourspace
inline vector interpolateRGB
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2
)
{
    return ((1-s)*rgb1 + s*rgb2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace colourTools
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
