/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

Typedefs
    CGAL3DKernel

Description
    A typedef for selecting inexact (faster) or exact (slower) 3D model.

\*---------------------------------------------------------------------------*/

#ifndef Foam_CGAL3DKernel_H
#define Foam_CGAL3DKernel_H

// Silence boost bind deprecation warnings (before CGAL-5.2.1)
#include "CGAL/version.h"
#if defined(CGAL_VERSION_NR) && (CGAL_VERSION_NR < 1050211000)
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#endif
#pragma clang diagnostic ignored "-Wbitwise-instead-of-logical"

// ------------------------------------------------------------------------- //

#ifdef CGAL_INEXACT

    // Fast kernel using a double as the storage type
    #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

#else

    // Very robust but expensive kernel
    #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
