/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Macros for accessing List elements

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListLoopM_H
#define Foam_ListLoopM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Element access looping

// Initial non-const access to list
#define List_ACCESS(type, f, fp)                \
    type* const __restrict__ fp = (f).begin()

// Initial const access to list
#define List_CONST_ACCESS(type, f, fp)          \
    const type* const __restrict__ fp = (f).begin()

// Loop over all elements
#define List_FOR_ALL(f, i)                      \
        const label _n##i = (f).size();         \
        for (label i=0; i<_n##i; ++i)

// Current element (non-const access)
#define List_ELEM(fp, i)  (fp[i])

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
