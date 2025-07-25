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

Note
    Some implementation details from VTK vtkColorTransferFunction.cxx
    vtkMath.cxx

\*---------------------------------------------------------------------------*/

#include "colourTools.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static constexpr scalar oneThird = 1.0 / 3.0;
static constexpr scalar oneSixth = 1.0 / 6.0;
static constexpr scalar twoThird = 2.0 / 3.0;
static constexpr scalar fiveSixth = 5.0 / 6.0;


// Compute HSV from RGB
static inline void RGB_to_HSV
(
    const scalar r, const scalar g, const scalar b,
    scalar& h, scalar& s, scalar& v
)
{
    scalar cmin = r, cmax = r;

    if (g > cmax)
    {
        cmax = g;
    }
    else if (g < cmin)
    {
        cmin = g;
    }
    if (b > cmax)
    {
        cmax = b;
    }
    else if (b < cmin)
    {
        cmin = b;
    }

    v = cmax;
    s = (v > 0.0) ? ((cmax - cmin) / cmax) : 0.0;

    if (s > 0.0)
    {
        if (r == cmax)
        {
            h = oneSixth * (g - b) / (cmax - cmin);
        }
        else if (g == cmax)
        {
            h = oneThird + oneSixth * (b - r) / (cmax - cmin);
        }
        else
        {
            h = twoThird + oneSixth * (r - g) / (cmax - cmin);
        }

        if (h < 0.0)
        {
            h += 1.0;
        }
    }
    else
    {
        h = 0.0;
    }
}


// Compute HSV to RGB
static inline void HSV_to_RGB
(
    const scalar h, const scalar s, const scalar v,
    scalar& r, scalar& g, scalar& b
)
{
    if (h > oneSixth && h <= oneThird)
    {
        // green/red
        g = 1.0;
        r = (oneThird - h) / oneSixth;
        b = 0.0;
    }
    else if (h > oneThird && h <= 0.5)
    {
        // green/blue
        g = 1.0;
        b = (h - oneThird) / oneSixth;
        r = 0.0;
    }
    else if (h > 0.5 && h <= twoThird)
    {
        // blue/green
        b = 1.0;
        g = (twoThird - h) / oneSixth;
        r = 0.0;
    }
    else if (h > twoThird && h <= fiveSixth)
    {
        // blue/red
        b = 1.0;
        r = (h - twoThird) / oneSixth;
        g = 0.0;
    }
    else if (h > fiveSixth && h <= 1.0)
    {
        // red/blue
        r = 1.0;
        b = (1.0 - h) / oneSixth;
        g = 0.0;
    }
    else
    {
        // red/green
        r = 1.0;
        g = h / oneSixth;
        b = 0.0;
    }

    // Add saturation
    r = (s * r + (1.0 - s));
    g = (s * g + (1.0 - s));
    b = (s * b + (1.0 - s));

    r *= v;
    g *= v;
    b *= v;
}


// Intermediate calculation to XYZ
static inline scalar to_XYZ(scalar val)
{
    const scalar p3 = pow3(val);
    return (p3 > 0.008856 ? p3 : (val - 16.0 / 116.0) / 7.787);
}


// Intermediate calculation from XYZ
static inline scalar from_XYZ(scalar val)
{
    return (val > 0.008856) ? cbrt(val) : (7.787 * val) + (16.0 / 116.0);
}

// Observer= 2 deg Illuminant= D65
static constexpr scalar ref_X = 0.9505;
static constexpr scalar ref_Y = 1.000;
static constexpr scalar ref_Z = 1.089;

static inline void LAB_to_XYZ
(
    const scalar L, const scalar a, const scalar b,
    scalar& x, scalar& y, scalar& z
)
{
    const scalar var_Y = (L + 16.0) / 116.0;
    const scalar var_X = a / 500 + var_Y;
    const scalar var_Z = var_Y - b / 200;

    x = ref_X * to_XYZ(var_X);
    y = ref_Y * to_XYZ(var_Y);
    z = ref_Z * to_XYZ(var_Z);
}


static inline void XYZ_to_LAB
(
    const scalar x, const scalar y, const scalar z,
    scalar& L, scalar& a, scalar& b
)
{
    const scalar var_X = from_XYZ(x / ref_X);
    const scalar var_Y = from_XYZ(y / ref_Y);
    const scalar var_Z = from_XYZ(z / ref_Z);

    L = (116 * var_Y) - 16;
    a = 500 * (var_X - var_Y);
    b = 200 * (var_Y - var_Z);
}



// "Gamma correction" specified by the sRGB color space.

static inline scalar gamma_from_xyz(const scalar val)
{
    return
    (
        val > 0.0031308
      ? (1.055 * (pow(val, 1.0/2.4)) - 0.055)
      : 12.92 * val
    );
}


static inline scalar gamma_to_xyz(const scalar val)
{
    return
    (
        val > 0.04045
      ? (pow((val + 0.055) / 1.055, 2.4))
      : val / 12.92
    );
}



static inline void XYZ_to_RGB
(
    const scalar x, const scalar y, const scalar z,
    scalar& r, scalar& g, scalar& b
)
{
    r = gamma_from_xyz(x *  3.2406 + y * -1.5372 + z * -0.4986);
    g = gamma_from_xyz(x * -0.9689 + y *  1.8758 + z *  0.0415);
    b = gamma_from_xyz(x *  0.0557 + y * -0.2040 + z *  1.0570);

    // Clip colour range
    scalar cmax = r;
    if (cmax < g) cmax = g;
    if (cmax < b) cmax = b;
    if (cmax > 1.0)
    {
        r /= cmax;
        g /= cmax;
        b /= cmax;
    }

    if (r < 0) r = 0;
    if (g < 0) g = 0;
    if (b < 0) b = 0;
}


static inline void RGB_to_XYZ
(
    scalar r, scalar g, scalar b,
    scalar& x, scalar& y, scalar& z
)
{
    r = gamma_to_xyz(r);
    g = gamma_to_xyz(g);
    b = gamma_to_xyz(b);

    x = r * 0.4124 + g * 0.3576 + b * 0.1805;
    y = r * 0.2126 + g * 0.7152 + b * 0.0722;
    z = r * 0.0193 + g * 0.1192 + b * 0.9505;
}



//- Convert to special polar version of CIELAB
//  (for creating continuous diverging color maps).
inline void labToMsh(const vector& lab, vector& msh)
{
    const scalar& L = lab[0];
    const scalar& a = lab[1];
    const scalar& b = lab[2];

    msh[0] = sqrt(L*L + a*a + b*b);
    msh[1] = (msh[0] > 0.001) ? acos(L / msh[0]) : 0.0;
    msh[2] = (msh[1] > 0.001) ? atan2(b,a) : 0.0;
}


//- Convert from special polar version of CIELAB
inline void mshToLab(const vector& msh, vector& lab)
{
    lab[0] = msh[0]*cos(msh[1]);
    lab[1] = msh[0]*sin(msh[1])*cos(msh[2]);
    lab[2] = msh[0]*sin(msh[1])*sin(msh[2]);
}


// Return the smallest angle between the two
static inline scalar angleDiff(scalar angle1, scalar angle2)
{
    scalar adiff = angle1 - angle1;
    if (adiff < 0.0) adiff = -adiff;

    while (adiff >= mathematical::twoPi) adiff -= mathematical::twoPi;
    if (adiff > mathematical::pi) adiff = (mathematical::twoPi - adiff);
    return adiff;
}


// For the case when interpolating from a saturated color to an unsaturated
// color, find a hue for the unsaturated color that makes sense.
static inline scalar adjustHue(const vector& msh, scalar unsatM)
{
    if (msh[0] >= unsatM - 0.1)
    {
        // The best we can do is hold hue constant.
        return msh[2];
    }

    // This equation is designed to make the perceptual change of the
    // interpolation to be close to constant.
    const scalar hueSpin =
        msh[1]*sqrt(unsatM*unsatM - msh[0]*msh[0]) / (msh[0]*sin(msh[1]));

    // Spin hue away from 0 except in purple hues.
    if (msh[2] > -0.3*mathematical::pi)
    {
        return msh[2] + hueSpin;
    }
    else
    {
        return msh[2] - hueSpin;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::colourTools::rgbToHsv(const vector& rgb, vector& hsv)
{
    RGB_to_HSV(rgb[0], rgb[1], rgb[2], hsv[0], hsv[1], hsv[2]);
}

void Foam::colourTools::hsvToRgb(const vector& hsv, vector& rgb)
{
    HSV_to_RGB(hsv[0], hsv[1], hsv[2], rgb[0], rgb[1], rgb[2]);
}


void Foam::colourTools::rgbToXyz(const vector& rgb, vector& xyz)
{
    RGB_to_XYZ(rgb[0], rgb[1], rgb[2], xyz[0], xyz[1], xyz[2]);
}

void Foam::colourTools::xyzToRgb(const vector& xyz, vector& rgb)
{
    XYZ_to_RGB(xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
}


void Foam::colourTools::labToXyz(const vector& lab, vector& xyz)
{
    LAB_to_XYZ(lab[0], lab[1], lab[2], xyz[0], xyz[1], xyz[2]);
}


void Foam::colourTools::xyzToLab(const vector& xyz, vector& lab)
{
    XYZ_to_LAB(xyz[0], xyz[1], xyz[2], lab[0], lab[1], lab[2]);
}


void Foam::colourTools::rgbToLab(const vector& rgb, vector& lab)
{
    vector xyz;
    RGB_to_XYZ(rgb[0], rgb[1], rgb[2], xyz[0], xyz[1], xyz[2]);
    XYZ_to_LAB(xyz[0], xyz[1], xyz[2], lab[0], lab[1], lab[2]);
}


void Foam::colourTools::labToRgb(const vector& lab, vector& rgb)
{
    vector xyz;
    labToXyz(lab, xyz);
    xyzToRgb(xyz, rgb);
}


void Foam::colourTools::interpolateDiverging
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2,
    vector& result
)
{
    vector lab1, lab2;
    rgbToLab(rgb1, lab1);
    rgbToLab(rgb2, lab2);

    vector msh1, msh2;
    labToMsh(lab1, msh1);
    labToMsh(lab2, msh2);

    // If the endpoints are distinct saturated colors,
    // then place white in between them.
    if
    (
        msh1[1] > 0.05
     && msh2[1] > 0.05
     && angleDiff(msh1[2], msh2[2]) > mathematical::pi/3.0
    )
    {
        // Insert the white midpoint by setting one end to white and
        // adjusting the scalar value.

        scalar Mmid = std::max(msh1[0], msh2[0]);
        Mmid = std::max(scalar(88.0), Mmid);
        if (s < 0.5)
        {
            msh2[0] = Mmid; msh2[1] = 0; msh2[2] = 0;
            s = 2.0*s;
        }
        else
        {
            msh1[0] = Mmid; msh1[1] = 0; msh1[2] = 0;
            s = 2.0*s - 1.0;
        }
    }

    // If one color has no saturation, then its hue value is invalid.
    // In this case, we want to set it to something logical so the
    // interpolation of hue makes sense.
    if ((msh1[1] < 0.05) && (msh2[1] > 0.05))
    {
        msh1[2] = adjustHue(msh2, msh1[0]);
    }
    else if ((msh2[1] < 0.05) && (msh1[1] > 0.05))
    {
        msh2[2] = adjustHue(msh1, msh2[0]);
    }

    // Msh tmp
    vector mshTmp((1-s)*msh1 + s*msh2);

    // Convert back to RGB
    vector lab;
    mshToLab(mshTmp, lab);
    labToRgb(lab, result);
}


void Foam::colourTools::interpolateHSV
(
    scalar s,
    const vector& rgb1,
    const vector& rgb2,
    vector& result
)
{
    vector hsv1, hsv2;
    rgbToHsv(rgb1, hsv1);
    rgbToHsv(rgb2, hsv2);

    // Wrap HSV?
    if (hsv1[0] - hsv2[0] > 0.5 || hsv2[0] - hsv1[0] > 0.5)
    {
        if (hsv1[0] > hsv2[0])
        {
            hsv1[0] -= 1.0;
        }
        else
        {
            hsv2[0] -= 1.0;
        }
    }

    vector hsvTmp((1-s)*hsv1 + s*hsv2);

    if (hsvTmp[0] < 0.0)
    {
        hsvTmp[0] += 1.0;
    }

    // Convert back to RGB
    hsvToRgb(hsvTmp, result);
}


// ************************************************************************* //
