/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "functionObjectBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileName
Foam::functionObjects::runTimePostPro::functionObjectBase::getFileName
(
    const word& keyword,
    const word& subDictName
) const
{
    dictionary dict;
    state_.getObjectDict(functionObjectName_, subDictName, dict);

    fileName f;
    if (dict.readIfPresent<fileName>(keyword, f))
    {
        f.expand();
    }

    return f;
}


bool Foam::functionObjects::runTimePostPro::functionObjectBase::removeFile
(
    const word& keyword,
    const word& subDictName
)
{
    // Foam::rm() ignores empty names etc.

    if (Pstream::master())
    {
        return Foam::rm(getFileName(keyword, subDictName));
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectBase::functionObjectBase
(
    const stateFunctionObject& state,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    fieldVisualisationBase(dict, colours),
    state_(state),
    functionObjectName_(dict.get<word>("functionObject")),
    liveObject_(dict.getOrDefault("liveObject", true)),
    clearObjects_(dict.getOrDefault("clearObjects", false))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::functionObjectBase::clear()
{
    return clearObjects_;
}


// ************************************************************************* //
