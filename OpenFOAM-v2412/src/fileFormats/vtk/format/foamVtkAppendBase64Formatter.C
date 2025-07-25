/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "foamVtkAppendBase64Formatter.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::appendBase64Formatter::name_ = "append";

const Foam::vtk::outputOptions
Foam::vtk::appendBase64Formatter::opts_(formatType::APPEND_BASE64);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::appendBase64Formatter::appendBase64Formatter
(
    std::ostream& os
)
:
    foamVtkBase64Layer(os),
    offset_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::appendBase64Formatter::~appendBase64Formatter()
{
    base64Layer::close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::appendBase64Formatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::appendBase64Formatter::name() const
{
    return name_;
}


uint64_t Foam::vtk::appendBase64Formatter::offset(const uint64_t numbytes)
{
    uint64_t prev = offset_;

    if (formatter::npos != numbytes)
    {
        offset_ += this->encodedLength(sizeof(uint64_t) + numbytes);
    }
    return prev;
}


// ************************************************************************* //
