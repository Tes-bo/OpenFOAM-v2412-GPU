/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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

#include "Helmholtz.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "wallFvPatch.H"
#include "DynamicList.H"
#include "fvMeshSubset.H"
#include "fvm.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Helmholtz, 1);
    addToRunTimeSelectionTable(regularisationPDE, Helmholtz, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Helmholtz::solveEqn
(
    const volScalarField& aTilda,
    const scalarField& source,
    scalarField& result,
    const bool isTopoField,
    const regularisationRadius& radius,
    const scalar minSetValue,
    const bool fixATildaValues
)
{
    const fvMesh& mesh = aTilda.internalField().mesh();
    // Convergence criteria
    const label iters = dict_.getOrDefault<label>("iters", 500);
    const scalar tolerance = dict_.getOrDefault<scalar>("tolerance", 1.e-06);
    dimensionedScalar one("1", dimless, 1.);
    // Smoothed field
    volScalarField bTilda
    (
        IOobject
        (
            "bTilda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        (
            growFromWalls_ ?
            fixedValueFvPatchScalarField::typeName :
            zeroGradientFvPatchScalarField::typeName
        )
    );
    // If solution corresponds to the topology porosity field, modify boundary
    // conditions accordingly
    if (isTopoField && growFromWalls_)
    {
        // Apply a unit alphaTilda value next to all walls
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            if (isA<wallFvPatch>(patch))
            {
                bTilda.boundaryFieldRef()[patchI] == wallValue_;
            }
        }
    }
    // Source field
    DimensionedField<scalar, volMesh> sourceField
    (
        IOobject
        (
            "source",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless,
        source
    );

    for (label iter = 0; iter < iters; ++iter)
    {
        fvScalarMatrix smoothEqn
        (
            fvm::Sp(one, bTilda)
          ==
            sourceField
        );
        radius.addRegularisationTerm(smoothEqn, isTopoField);

        // Set solution in given zones
        if (fixATildaValues)
        {
            setValues(smoothEqn, isTopoField, minSetValue);
        }

        const scalar residual(mag(smoothEqn.solve().initialResidual()));

//      if (isTopoField)
//      {
//          bound(bTilda, dimensionedScalar(bTilda.dimensions(), minSetValue));
//      }

        // Print execution time
        mesh.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance)
        {
            Info<< "\n***Reached regularisation equation convergence limit, "
                   "iteration " << iter << "***\n\n";
            break;
        }
    }

    // Replace field with its regularised counterpart
    result = bTilda.primitiveField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Helmholtz::Helmholtz
(
    const fvMesh& mesh,
    const dictionary& dict,
    const topOZones& zones
)
:
    regularisationPDE(mesh, dict, zones),
    solveOnActiveCells_(dict.getOrDefault<bool>("solveOnActiveCells", false)),
    wallValue_(dict.getOrDefault<scalar>("wallValue", 1))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Helmholtz::regularise
(
    const volScalarField& aTilda,
    const scalarField& source,
    scalarField& result,
    const bool isTopoField,
    const regularisationRadius& radius,
    const scalar minSetValue,
    const bool fixATildaValues
)
{
    // Set values for all constant cells here.
    // If a subsetMesh is used, cells outside its domain will not be changed,
    // potentially leading to fixed cells not getting their correct values
    if (fixATildaValues)
    {
        DynamicList<label> cells(0);
        DynamicList<scalar> values(0);
        setValues(mesh_, cells, values, isTopoField);
        result.rmap(values, cells);
    }

    /*
    // Solve the regularisation equation on a mesh including only the active
    // cells, if needed
    if (solveOnActiveCells_)
    {
        const labelList& activeZones = zones_.adjointPorousZoneIDs();
        if (!activeZones.empty())
        {
            Info<< "Solving regularisation equation on active cells only"
                << endl;
            DynamicList<label> allActiveCells(0);
            for (const label zI : activeZones)
            {
                allActiveCells.append(mesh_.cellZones()[zI]);
            }
            fvMeshSubset::exposedPatchType = wallPolyPatch::typeName;
            fvMeshSubset subSetMesh(mesh_, allActiveCells);
            fvMesh& subMesh = subSetMesh.subMesh();

            schemesLookup& fvSchemes = static_cast<schemesLookup&>(subMesh);
            fvSchemes.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            fvSchemes.read();

            fvSolution& solution = static_cast<fvSolution&>(subMesh);
            solution.readOpt(IOobject::MUST_READ_IF_MODIFIED);
            solution.read();

            const labelList& cellMap = subSetMesh.cellMap();
            // Map input fields to subSetMesh fields
            volScalarField aTildaSub(subSetMesh.interpolate(aTilda));
            scalarField sourceSub(source, cellMap);
            scalarField resultSub(result, cellMap);
            // Solve the regularisation equation
            solveEqn
            (
                aTildaSub, sourceSub, resultSub,
                isTopoField, radius, minSetValue, fixATildaValues
            );
            // Map result back to original field
            result.rmap(resultSub, cellMap);
            Info<< "min max " << gMin(result) << " " << gMax(result) << endl;
            return;
        }
    }
    */
    solveEqn
    (
        aTilda,
        source,
        result,
        isTopoField,
        radius,
        minSetValue,
        fixATildaValues
    );
}


// ************************************************************************* //
