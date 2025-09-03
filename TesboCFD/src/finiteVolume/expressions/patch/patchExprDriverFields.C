/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "patchExprDriver.H"
#include "fvPatch.H"
#include "error.H"

// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::Field<bool>>
Foam::expressions::patchExpr::parseDriver::getSurfaceField<bool>
(
    const word& name
)
{
    return getVariable<bool>(name, this->size());
}


template<>
Foam::tmp<Foam::Field<bool>>
Foam::expressions::patchExpr::parseDriver::getPointField<bool>
(
    const word& name
)
{
    return getVariable<bool>(name, this->pointSize());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::field_cellSelection
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    refPtr<labelList> tselected;
    switch (setType)
    {
        case topoSetSource::sourceType::CELLZONE_SOURCE:
        case topoSetSource::sourceType::CELLSET_SOURCE:
        {
            tselected = getTopoSetLabels(name, setType);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unexpected sourceType: " << int(setType) << nl
                << exit(FatalError);
            break;
        }
    }

    // Not particularly efficient...
    labelHashSet inSelection(tselected());

    const labelList& faceCells = patch_.faceCells();
    auto tresult = tmp<boolField>::New(this->size(), false);
    auto& result = tresult.ref();

    forAll(result, facei)
    {
        if (inSelection.found(faceCells[facei]))
        {
            result[facei] = true;
        }
    }

    return tresult;
}


Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::field_faceSelection
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    refPtr<labelList> tselected;
    switch (setType)
    {
        case topoSetSource::sourceType::FACESET_SOURCE:
        case topoSetSource::sourceType::FACEZONE_SOURCE:
        {
            tselected = getTopoSetLabels(name, setType);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unexpected sourceType: " << int(setType) << nl
                << exit(FatalError);
            break;
        }
    }

    // Not particularly efficient...
    labelHashSet inSelection(tselected());

    const label patchStart = patch_.start();

    auto tresult = tmp<boolField>::New(this->size(), false);
    auto& result = tresult.ref();

    forAll(result, facei)
    {
        if (inSelection.found(facei + patchStart))
        {
            result[facei] = true;
        }
    }

    return tresult;
}


Foam::tmp<Foam::scalarField>
Foam::expressions::patchExpr::parseDriver::field_faceArea() const
{
    return patch_.magSf();
}


Foam::tmp<Foam::vectorField>
Foam::expressions::patchExpr::parseDriver::field_faceCentre() const
{
    return patch_.Cf();
}


Foam::tmp<Foam::vectorField>
Foam::expressions::patchExpr::parseDriver::field_areaNormal() const
{
    return patch_.Sf();
}


Foam::tmp<Foam::vectorField>
Foam::expressions::patchExpr::parseDriver::field_unitNormal() const
{
    return patch_.nf();
}


Foam::tmp<Foam::vectorField>
Foam::expressions::patchExpr::parseDriver::field_pointField() const
{
    return patch_.patch().localPoints();
}


Foam::tmp<Foam::scalarField>
Foam::expressions::patchExpr::parseDriver::field_rand
(
    label seed,
    bool gaussian
) const
{
    auto tfld = tmp<scalarField>::New(this->size());

    exprDriver::fill_random(tfld.ref(), seed, gaussian);

    return tfld;
}


// ************************************************************************* //
