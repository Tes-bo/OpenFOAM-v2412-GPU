/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::legacyAsciiFormatter::legacyName_ = "ASCII";

const Foam::vtk::outputOptions
Foam::vtk::legacyAsciiFormatter::opts_(formatType::LEGACY_ASCII);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::legacyAsciiFormatter::legacyAsciiFormatter
(
    std::ostream& os
)
:
    asciiFormatter(os)
{}


Foam::vtk::legacyAsciiFormatter::legacyAsciiFormatter
(
    std::ostream& os,
    unsigned prec
)
:
    asciiFormatter(os, prec)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::legacyAsciiFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::legacyAsciiFormatter::name() const
{
    return legacyName_;
}


const char* Foam::vtk::legacyAsciiFormatter::encoding() const
{
    return legacyName_;
}


// ************************************************************************* //
