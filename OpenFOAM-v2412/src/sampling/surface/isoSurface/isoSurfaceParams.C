/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "isoSurfaceParams.H"
#include "dictionary.H"
#include "Switch.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::isoSurfaceParams::algorithmType
>
Foam::isoSurfaceParams::algorithmNames
({
    { algorithmType::ALGO_DEFAULT, "default" },
    { algorithmType::ALGO_CELL, "cell" },
    { algorithmType::ALGO_POINT, "point" },
    { algorithmType::ALGO_TOPO, "topo" },
});


const Foam::Enum
<
    Foam::isoSurfaceParams::filterType
>
Foam::isoSurfaceParams::filterNames
({
    { filterType::NONE, "none" },
    { filterType::PARTIAL, "partial" },
    { filterType::FULL, "full" },
    { filterType::CLEAN, "clean" },

    { filterType::CELL, "cell" },
    { filterType::DIAGCELL, "diagcell" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::isoSurfaceParams::algorithmType
Foam::isoSurfaceParams::getAlgorithmType
(
    const dictionary& dict,
    const algorithmType deflt
)
{
    word enumName;
    if
    (
        !dict.readIfPresentCompat
        (
            "isoMethod", {{"isoAlgorithm", 0}},
            enumName, keyType::LITERAL
        )
    )
    {
        return deflt;
    }

    if (!algorithmNames.found(enumName))
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << (algorithmNames) << nl
            << exit(FatalIOError);
    }

    return algorithmNames[enumName];
}


Foam::isoSurfaceParams::filterType
Foam::isoSurfaceParams::getFilterType
(
    const dictionary& dict,
    const filterType deflt
)
{
    word enumName;
    if (!dict.readIfPresent("regularise", enumName, keyType::LITERAL))
    {
        return deflt;
    }

    // Try as bool/switch
    const Switch sw = Switch::find(enumName);

    if (sw.good())
    {
        return (sw ? deflt : filterType::NONE);
    }

    // As enum
    if (!filterNames.found(enumName))
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << (filterNames) << nl
            << exit(FatalIOError);
    }

    return filterNames[enumName];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceParams::isoSurfaceParams
(
    const algorithmType algo,
    const filterType filter
) noexcept
:
    algo_(algo),
    filter_(filter),
    snap_(true),
    mergeTol_(1e-6),
    clipBounds_()
{}


Foam::isoSurfaceParams::isoSurfaceParams
(
    const dictionary& dict,
    const isoSurfaceParams& params
)
:
    isoSurfaceParams(params)
{
    algo_ = getAlgorithmType(dict, algo_);
    filter_ = getFilterType(dict, filter_);
    snap_ = dict.getOrDefault("snap", true);
    dict.readIfPresent("mergeTol", mergeTol_);
    dict.readIfPresent("bounds", clipBounds_);
}


Foam::isoSurfaceParams::isoSurfaceParams
(
    const dictionary& dict,
    const algorithmType algo,
    const filterType filter
)
:
    isoSurfaceParams(dict, isoSurfaceParams(algo, filter))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isoSurfaceParams::setClipBounds(const boundBox& bb)
{
    clipBounds_ = bb;
}


void Foam::isoSurfaceParams::print(Ostream& os) const
{
    os  << " isoMethod:" << algorithmNames[algo_]
        << " regularise:" << filterNames[filter_]
        << " snap:" << snap_;
}


// ************************************************************************* //
