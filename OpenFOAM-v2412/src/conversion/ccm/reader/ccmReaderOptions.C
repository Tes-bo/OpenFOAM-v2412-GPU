/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ccmReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ccm::reader::options::options()
:
    keepFluid_(true),
    keepPorous_(true),
    keepSolid_(true),
    mergeInterfaces_(false),
    renameInterfaces_(true),
    removeBaffles_(false),
    useNumberedNames_(false),
    mergeTol_(0.05e-3),
    undefScalar_(NAN)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ccm::reader::options::keepFluid() const
{
    return keepFluid_;
}


bool Foam::ccm::reader::options::keepPorous() const
{
    return keepPorous_;
}


bool Foam::ccm::reader::options::keepSolid() const
{
    return keepSolid_;
}


bool Foam::ccm::reader::options::keptSomeRegion() const
{
    return keepFluid_ || keepPorous_ || !keepSolid_;
}


bool Foam::ccm::reader::options::mergeInterfaces() const
{
    return mergeInterfaces_;
}


bool Foam::ccm::reader::options::renameInterfaces() const
{
    return renameInterfaces_;
}


bool Foam::ccm::reader::options::removeBaffles() const
{
    return removeBaffles_;
}


bool Foam::ccm::reader::options::useNumberedNames() const
{
    return useNumberedNames_;
}


Foam::scalar Foam::ccm::reader::options::mergeTol() const
{
    return mergeTol_;
}


Foam::scalar Foam::ccm::reader::options::undefScalar() const
{
    return undefScalar_;
}



void Foam::ccm::reader::options::keepFluid(bool b)
{
    keepFluid_ = b;
}


void Foam::ccm::reader::options::keepPorous(bool b)
{
    keepPorous_ = b;
}


void Foam::ccm::reader::options::keepSolid(bool b)
{
    keepSolid_ = b;
}


void Foam::ccm::reader::options::mergeInterfaces(bool b)
{
    mergeInterfaces_ = b;
}


void Foam::ccm::reader::options::renameInterfaces(bool b)
{
    renameInterfaces_ = b;
}


void Foam::ccm::reader::options::removeBaffles(bool b)
{
    removeBaffles_ = b;
}


void Foam::ccm::reader::options::useNumberedNames(bool b)
{
    useNumberedNames_ = b;
}


void Foam::ccm::reader::options::mergeTol(const scalar& val)
{
    mergeTol_ = val;
}


void Foam::ccm::reader::options::undefScalar(const scalar& val)
{
    undefScalar_ = val;
}


// ************************************************************************* //
