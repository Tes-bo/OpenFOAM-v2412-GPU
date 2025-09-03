/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "IOstreamOption.H"
#include "IOstreams.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Default output precision
unsigned int Foam::IOstream::precision_
(
    Foam::debug::infoSwitch("writePrecision", 6)
);


// * * * * * * * * * * * * * * Global Variables  * * * * * * * * * * * * * * //

// Global IO streams

Foam::ISstream Foam::Sin(std::cin, "Sin");
Foam::OSstream Foam::Sout(std::cout, "Sout");
Foam::OSstream Foam::Serr(std::cerr, "Serr");
Foam::OFstream Foam::Snull(nullptr);  // A "/dev/null" equivalent

Foam::prefixOSstream Foam::Pout(std::cout, "Pout");
Foam::prefixOSstream Foam::Perr(std::cerr, "Perr");


// ************************************************************************* //
