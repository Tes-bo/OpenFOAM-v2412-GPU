/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

#include "regularisationPDE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regularisationPDE, 0);
    defineRunTimeSelectionTable(regularisationPDE, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::regularisationPDE::setValues
(
    fvScalarMatrix& bTildaEqn,
    const bool isTopoField,
    const scalar minSetValue
)
{
    const fvMesh& mesh = bTildaEqn.psi().internalField().mesh();
    DynamicList<label> cells(mesh.nCells()/100);
    DynamicList<scalar> values(mesh.nCells()/100);
    setValues(mesh, cells, values, isTopoField, minSetValue);
    bTildaEqn.setValues(cells, values);
}


void Foam::regularisationPDE::setValues
(
    const fvMesh& mesh,
    DynamicList<label>& cells,
    DynamicList<scalar>& values,
    bool isTopoField,
    const scalar minSetValue
)
{
    const labelList& IOcells = mesh.cellZones()[zones_.IOzoneID()];
    cells.append(IOcells);
    values.append(scalarField(IOcells.size(), minSetValue));

    forAll(zones_.fixedPorousZoneIDs(), zI)
    {
        const label cellZoneID = zones_.fixedPorousZoneIDs()[zI];
        const labelList& zoneCells = mesh.cellZones()[cellZoneID];
        cells.append(zoneCells);
        values.append
        (
            scalarField
            (
                zoneCells.size(),
                isTopoField ?
                zones_.fixedPorousValues()[zI] : scalar(0)
            )
        );
    }

    for (label cellZoneID : zones_.fixedZeroPorousZoneIDs())
    {
        const labelList& zoneCells = mesh.cellZones()[cellZoneID];
        cells.append(zoneCells);
        values.append(scalarField(zoneCells.size(), minSetValue));
    }
}


Foam::scalar Foam::regularisationPDE::computeRadius()
{
    scalar averageVol(gAverage(mesh_.V().field()));
    const Vector<label>& geometricD = mesh_.geometricD();
    const boundBox& bounds = mesh_.bounds();
    forAll(geometricD, iDir)
    {
        if (geometricD[iDir] == -1)
        {
            averageVol /= bounds.span()[iDir];
        }
    }
    scalar radius = pow(averageVol, scalar(1)/scalar(mesh_.nGeometricD()));
    Info<< "Computed a mean radius of " << radius << endl;
    return radius;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regularisationPDE::regularisationPDE
(
    const fvMesh& mesh,
    const dictionary& dict,
    const topOZones& zones
)
:
    mesh_(mesh),
    dict_(dict),
    zones_(zones),
    growFromWalls_(dict.getOrDefault<bool>("growFromWalls", false))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regularisationPDE> Foam::regularisationPDE::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const topOZones& zones
)
{
    const word modelType =
        dict.getOrDefault<word>("regularisationPDE", "Helmholtz");

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    Info<< "regularisationPDE type " << modelType << endl;

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "regularisationPDE",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<regularisationPDE>
    (
        ctorPtr(mesh, dict, zones)
    );
}


// ************************************************************************* //
