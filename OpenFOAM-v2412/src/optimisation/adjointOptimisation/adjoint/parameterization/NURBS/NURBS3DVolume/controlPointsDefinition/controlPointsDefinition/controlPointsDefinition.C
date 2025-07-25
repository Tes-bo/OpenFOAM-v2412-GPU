/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "controlPointsDefinition.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(controlPointsDefinition, 0);
    defineRunTimeSelectionTable(controlPointsDefinition, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::controlPointsDefinition::transformControlPoints
(
    const vector& geometryMin,
    const vector& geometryMax
)
{
    const dictionary& dict = box_.dict();
    // Translation vector
    vector position(dict.get<vector>("translation"));

    // Rotation vector
    vector rotation(dict.get<vector>("rotation"));
    const scalar deg2rad(constant::mathematical::pi/180.0);
    rotation *= deg2rad;

    // Scaling
    vector scale(dict.get<vector>("scale"));

    // Scale box
    cps_ = cmptMultiply(cps_, scale);

    // Rotation matrices
    tensor Rx
    (
        1, 0                  ,    0,
        0, ::cos(rotation.x()), -::sin(rotation.x()),
        0, ::sin(rotation.x()),  ::cos(rotation.x())
    );
    tensor Ry
    (
        ::cos(rotation.y()), 0,  ::sin(rotation.y()),
          0                , 1,  0,
       -::sin(rotation.y()), 0,  ::cos(rotation.y())
    );
    tensor Rz
    (
        ::cos(rotation.z()), -::sin(rotation.z()), 0,
        ::sin(rotation.z()),  ::cos(rotation.z()), 0,
        0,                      0                , 1
    );

    // Comined rotation matrix
    tensor R = (Rz & Rx) & Ry;

    // Rotate cps
    cps_ = R & cps_;

    // Translate cps
    cps_ += position;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controlPointsDefinition::controlPointsDefinition
(
    NURBS3DVolume& box
)
:
    box_(box),
    cps_(box.getControlPoints())
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::controlPointsDefinition> Foam::controlPointsDefinition::New
(
    NURBS3DVolume& box
)
{
    const dictionary& dict = box.dict();
    const word type(dict.get<word>("controlPointsDefinition"));

    Info<< "controlPointsDefinition type : " << type << endl;

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "controlPointsDefinition",
            type,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<controlPointsDefinition>(ctorPtr(box));
}


// ************************************************************************* //
