/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "adiosCoreWrite.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosFoam::adiosCoreWrite::putAttribute
(
    const fileName& name,
    const std::string& value
)
{
    writeIOPtr_->DefineAttribute<std::string>(name, value);
}


void Foam::adiosFoam::adiosCoreWrite::putIntAttribute
(
    const fileName& name,
    const int value
)
{
    writeIOPtr_->DefineAttribute<int>(name, value);
}


void Foam::adiosFoam::adiosCoreWrite::putDoubleAttribute
(
    const fileName& name,
    const double value
)
{
    writeIOPtr_->DefineAttribute<double>(name, value);
}


bool Foam::adiosFoam::adiosCoreWrite::putListAttribute
(
    const fileName& name,
    const UList<int>& list
)
{
    if (!list.empty())
    {
        writeIOPtr_->DefineAttribute<int>
        (
            name,
            list.cdata(),
            list.size()
        );
    }

    return !list.empty();
}


bool Foam::adiosFoam::adiosCoreWrite::putListAttribute
(
    const fileName& name,
    const UList<double>& list
)
{
    if (!list.empty())
    {
        writeIOPtr_->DefineAttribute<double>
        (
            name,
            list.cdata(),
            list.size()
        );
    }

    return !list.empty();
}


bool Foam::adiosFoam::adiosCoreWrite::putListAttribute
(
    const fileName& name,
    const UList<word>& list
)
{
    if (!list.empty())
    {
        writeIOPtr_->DefineAttribute<std::string>
        (
            name,
            list.cdata(),
            list.size()
        );
    }

    return !list.empty();
}


bool Foam::adiosFoam::adiosCoreWrite::putListAttribute
(
    const fileName& name,
    const UList<string>& list
)
{
    if (!list.empty())
    {
        writeIOPtr_->DefineAttribute<std::string>
        (
            name,
            list.cdata(),
            list.size()
        );
    }

    return !list.empty();
}


bool Foam::adiosFoam::adiosCoreWrite::putDimensionsAttribute
(
    const fileName& varName,
    const dimensionSet& dims
)
{
    const auto list = dims.values();

    writeIOPtr_->DefineAttribute<scalar>
    (
        varName/"dimensions",
        list.cdata(),
        list.size()
    );

    return true;
}


// ************************************************************************* //
