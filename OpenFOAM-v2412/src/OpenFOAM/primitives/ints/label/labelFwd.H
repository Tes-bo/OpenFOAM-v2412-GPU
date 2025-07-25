/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

Description
    Typedefs for label/uLabel without requiring label.H

    The definitions of label/uLabel as an integral value depend on the
    pre-processor macro WM_LABEL_SIZE.

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_labelFwd_H
#define Foam_primitives_labelFwd_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#if WM_LABEL_SIZE == 32

typedef int32_t label;
typedef uint32_t uLabel;

#elif WM_LABEL_SIZE == 64

typedef int64_t label;
typedef uint64_t uLabel;

#else
    #error "WM_LABEL_SIZE must be set to either 32 or 64"
#endif

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
