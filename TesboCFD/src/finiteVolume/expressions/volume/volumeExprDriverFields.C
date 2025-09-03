/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "volumeExprDriver.H"
#include "fvPatch.H"
#include "error.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::expressions::volumeExpr::parseDriver::field_cellSelection
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    auto tresult = volScalarField::New
    (
        "selected",
        mesh(),
        dimensionedScalar(Zero)
    );

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
    const auto& selected = tselected();

    auto& fld = tresult.ref().primitiveFieldRef();
    UIndirectList<scalar>(fld, selected) = scalar(1);

    return tresult;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::expressions::volumeExpr::parseDriver::field_faceSelection
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    auto tresult = surfaceScalarField::New
    (
        "selected",
        mesh(),
        dimensionedScalar(Zero)
    );

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
    const auto& selected = tselected();

    const auto& bmesh = mesh().boundaryMesh();

    auto& result = tresult.ref();
    auto& fld = result.primitiveFieldRef();
    auto& bfld = result.boundaryFieldRef();

    label nErrors = 0;

    for (const label facei : selected)
    {
        if (facei < mesh().nInternalFaces())
        {
            fld[facei] = scalar(1);
        }
        else
        {
            const label patchi = bmesh.whichPatch(facei);

            if (patchi < 0)
            {
                ++nErrors;
            }
            else
            {
                bfld[patchi][facei-bmesh[patchi].start()] = scalar(1);
            }
        }
    }

    if (nErrors)
    {
        WarningInFunction
            << "The faceSet/faceZone " << name << " contained "
            << nErrors << " faces outside of the addressing range" << nl
            << nl;
    }


    return tresult;
}


Foam::tmp<Foam::pointScalarField>
Foam::expressions::volumeExpr::parseDriver::field_pointSelection
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    auto tresult = pointScalarField::New
    (
        "selected",
        pointMesh::New(mesh()),
        dimensionedScalar(Zero)
    );

    refPtr<labelList> tselected;
    switch (setType)
    {
        case topoSetSource::sourceType::POINTSET_SOURCE:
        case topoSetSource::sourceType::POINTZONE_SOURCE:
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
    const auto& selected = tselected();

    auto& fld = tresult.ref().primitiveFieldRef();
    UIndirectList<scalar>(fld, selected) = scalar(1);

    return tresult;
}



Foam::tmp<Foam::volScalarField>
Foam::expressions::volumeExpr::parseDriver::field_cellVolume() const
{
    return volScalarField::New
    (
        "vol",
        mesh(),
        dimVol,
        mesh().V()
    );
}


Foam::tmp<Foam::volVectorField>
Foam::expressions::volumeExpr::parseDriver::field_cellCentre() const
{
    return tmp<volVectorField>::New(mesh().C());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::expressions::volumeExpr::parseDriver::field_faceArea() const
{
    return surfaceScalarField::New
    (
        "face",
        mesh(),
        dimless,
        mesh().magSf()
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::expressions::volumeExpr::parseDriver::field_faceCentre() const
{
    return surfaceVectorField::New
    (
        "fpos",
        mesh(),
        dimless,
        mesh().Cf()
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::expressions::volumeExpr::parseDriver::field_areaNormal() const
{
    return surfaceVectorField::New
    (
        "face",
        mesh(),
        dimless,
        mesh().Sf()
    );
}


Foam::tmp<Foam::pointVectorField>
Foam::expressions::volumeExpr::parseDriver::field_pointField() const
{
    return pointVectorField::New
    (
        "pts",
        pointMesh::New(mesh()),
        dimless,
        mesh().points()
    );
}


Foam::tmp<Foam::volScalarField>
Foam::expressions::volumeExpr::parseDriver::field_rand
(
    label seed,
    bool gaussian
) const
{
    auto tfld = volScalarField::New("rand", mesh(), dimless);

    exprDriver::fill_random(tfld.ref().primitiveFieldRef(), seed, gaussian);

    return tfld;
}


// ************************************************************************* //
