/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "faMesh.H"
#include "areaFields.H"
#include "coupledFaPatchField.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebug(areaScalarField::Internal, 0);
defineTemplateTypeNameAndDebug(areaVectorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(areaScalarField, 0);
defineTemplateTypeNameAndDebug(areaVectorField, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField, 0);
defineTemplateTypeNameAndDebug(areaTensorField, 0);

defineTemplateDebugSwitchWithName
(
    areaScalarField::Boundary,
    "areaScalarField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaVectorField::Boundary,
    "areaVectorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaSphericalTensorField::Boundary,
    "areaSphericalTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaSymmTensorField::Boundary,
    "areaSymmTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaTensorField::Boundary,
    "areaTensorField::Boundary",
    0
);

template<> scalar areaScalarField::Boundary::tolerance
(
    debug::floatOptimisationSwitch("areaScalarField::Boundary::tolerance", 0)
);
registerOptSwitch
(
    "areaScalarField::Boundary::tolerance",
    scalar,
    Foam::areaScalarField::Boundary::tolerance
);

template<> scalar areaVectorField::Boundary::tolerance
(
    debug::floatOptimisationSwitch("areaVectorField::Boundary::tolerance", 0)
);
registerOptSwitch
(
    "areaVectorField::Boundary::tolerance",
    scalar,
    Foam::areaVectorField::Boundary::tolerance
);

template<> scalar areaSphericalTensorField::Boundary::tolerance
(
    debug::floatOptimisationSwitch
    (
        "areaSphericalTensorField::Boundary::tolerance",
        0
    )
);
registerOptSwitch
(
    "areaSphericalTensorField::Boundary::tolerance",
    scalar,
    Foam::areaSphericalTensorField::Boundary::tolerance
);

template<> scalar areaSymmTensorField::Boundary::tolerance
(
    debug::floatOptimisationSwitch
    (
        "areaSymmTensorField::Boundary::tolerance",
        0
    )
);
registerOptSwitch
(
    "areaSymmTensorField::Boundary::tolerance",
    scalar,
    Foam::areaSymmTensorField::Boundary::tolerance
);

template<> scalar areaTensorField::Boundary::tolerance
(
    debug::floatOptimisationSwitch("areaTensorField::Boundary::tolerance", 0)
);
registerOptSwitch
(
    "areaTensorField::Boundary::tolerance",
    scalar,
    Foam::areaTensorField::Boundary::tolerance
);


// Local-ops consistency enforcing

template<> int areaScalarField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "areaScalarField::Boundary::localConsistency",
    int,
    Foam::areaScalarField::Boundary::localConsistency
);

template<> int areaVectorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "areaVectorField::Boundary::localConsistency",
    int,
    Foam::areaVectorField::Boundary::localConsistency
);

template<> int areaSphericalTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "areaSphericalTensorField::Boundary::localConsistency",
    int,
    Foam::areaSphericalTensorField::Boundary::localConsistency
);

template<> int areaSymmTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "areaSymmTensorField::Boundary::localConsistency",
    int,
    Foam::areaSymmTensorField::Boundary::localConsistency
);

template<> int areaTensorField::Boundary::localConsistency
(
    debug::optimisationSwitch("localConsistency", 1)
);
registerOptSwitch
(
    "areaTensorField::Boundary::localConsistency",
    int,
    Foam::areaTensorField::Boundary::localConsistency
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations

namespace Foam
{

template<>
tmp<GeometricField<scalar, faPatchField, areaMesh>>
GeometricField<scalar, faPatchField, areaMesh>::component
(
    const direction
) const
{
    return *this;
}

template<>
void GeometricField<scalar, faPatchField, areaMesh>::replace
(
    const direction,
    const GeometricField<scalar, faPatchField, areaMesh>& gsf
)
{
    *this == gsf;
}

template<>
bool GeometricBoundaryField<scalar, faPatchField, areaMesh>::check() const
{
    return checkConsistency<coupledFaPatchField<scalar>>
    (
        Foam::areaScalarField::Boundary::tolerance,
       !(debug&4)       // make into warning if debug&4
    );
}

template<>
bool GeometricBoundaryField<vector, faPatchField, areaMesh>::check() const
{
    return checkConsistency<coupledFaPatchField<vector>>
    (
        Foam::areaScalarField::Boundary::tolerance,
       !(debug&4)       // make into warning if debug&4
    );
}

template<>
bool GeometricBoundaryField<sphericalTensor, faPatchField, areaMesh>::check
() const
{
    return checkConsistency<coupledFaPatchField<sphericalTensor>>
    (
        Foam::areaScalarField::Boundary::tolerance,
       !(debug&4)       // make into warning if debug&4
    );
}

template<>
bool GeometricBoundaryField<symmTensor, faPatchField, areaMesh>::check() const
{
    return checkConsistency<coupledFaPatchField<symmTensor>>
    (
        Foam::areaScalarField::Boundary::tolerance,
       !(debug&4)       // make into warning if debug&4
    );
}

template<>
bool GeometricBoundaryField<tensor, faPatchField, areaMesh>::check() const
{
    return checkConsistency<coupledFaPatchField<tensor>>
    (
        Foam::areaScalarField::Boundary::tolerance,
       !(debug&4)       // make into warning if debug&4
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Note hard-coded values are more reliable than other alternatives

const Foam::wordList Foam::fieldTypes::area
({
    "areaScalarField",
    "areaVectorField",
    "areaSphericalTensorField",
    "areaSymmTensorField",
    "areaTensorField"
});

const Foam::wordList Foam::fieldTypes::area_internal
({
    "areaScalarField::Internal",
    "areaVectorField::Internal",
    "areaSphericalTensorField::Internal",
    "areaSymmTensorField::Internal",
    "areaTensorField::Internal"
});


// ************************************************************************* //
