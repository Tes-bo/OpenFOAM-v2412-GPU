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

InNamespace
    Foam::stringOps

Description
    String expression evaluation.

\*---------------------------------------------------------------------------*/

#ifndef Foam_stringOpsEvaluate_H
#define Foam_stringOpsEvaluate_H

#include "labelFwd.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stringOps
{
    //- String evaluation with specified (positive, non-zero) field width
    string evaluate
    (
        label fieldWidth,
        const std::string& s,
        size_t pos = 0,
        size_t len = std::string::npos
    );

    //- A simple string evaluation that handles various basic
    //- expressions. For trivial input, use readScalar instead (faster).
    //
    //  The evaluation supports the following:
    //  - operations:  - + * / %
    //  - functions:  exp, log, log10, pow, sqrt, cbrt, sqr, mag, magSqr
    //  - trigonometric:  sin, cos, tan, asin, acos, atan, atan2, hypot
    //  - hyperbolic:  sinh, cosh, tanh
    //  - conversions:  degToRad, radToDeg
    //  - type conversion:  bool, mag
    //  - constants:  pi(), true, false
    //  - limits: neg, pos, neg0, pos0, sign, floor, ceil, round
    //  - other: rand, rand(seed)
    //  - logic: ! ? : ==  !=  <=  => < >
    //
    //  \note
    //  Unlike C/C++, the ternary and logical operations are \b not
    //  short-circuiting. So additional guards may be required.
    //
    //  \note
    //  The rand() function returns a uniform scalar on [0-1] interval
    //  and uses a constant seed.
    string evaluate
    (
        const std::string& s,
        size_t pos = 0,
        size_t len = std::string::npos
    );

} // End namespace stringOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
