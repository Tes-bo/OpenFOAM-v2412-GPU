/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "pointFields.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebug(pointScalarField::Internal, 0);
defineTemplateTypeNameAndDebug(pointVectorField::Internal, 0);
defineTemplateTypeNameAndDebug(pointSphericalTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(pointSymmTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(pointTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(pointScalarField, 0);
defineTemplateTypeNameAndDebug(pointVectorField, 0);
defineTemplateTypeNameAndDebug(pointSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensorField, 0);
defineTemplateTypeNameAndDebug(pointTensorField, 0);


defineTemplateDebugSwitchWithName
(
    pointScalarField::Boundary,
    "pointScalarField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    pointVectorField::Boundary,
    "pointVectorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    pointSphericalTensorField::Boundary,
    "pointSphericalTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    pointSymmTensorField::Boundary,
    "pointSymmTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    pointTensorField::Boundary,
    "pointTensorField::Boundary",
    0
);



// Local-ops consistency enforcing

template<> int pointScalarField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "localConsistency",
    int,
    Foam::pointScalarField::Boundary::localConsistency
);

template<> int pointVectorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "pointVectorField::Boundary::localConsistency",
    int,
    Foam::pointVectorField::Boundary::localConsistency
);

template<> int pointSphericalTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "pointSphericalTensorField::Boundary::localConsistency",
    int,
    Foam::pointSphericalTensorField::Boundary::localConsistency
);

template<> int pointSymmTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "pointSymmTensorField::Boundary::localConsistency",
    int,
    Foam::pointSymmTensorField::Boundary::localConsistency
);

template<> int pointTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "pointTensorField::Boundary::localConsistency",
    int,
    Foam::pointTensorField::Boundary::localConsistency
);


} // End namespace Foam


// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Note hard-coded values are more reliable than other alternatives

const Foam::wordList Foam::fieldTypes::point
({
    "pointScalarField",
    "pointVectorField",
    "pointSphericalTensorField",
    "pointSymmTensorField",
    "pointTensorField"
});


// ************************************************************************* //
