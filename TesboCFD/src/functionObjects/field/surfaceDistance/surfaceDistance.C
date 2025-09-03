/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "surfaceDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceDistance, 0);
    addToRunTimeSelectionTable(functionObject, surfaceDistance, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceDistance::surfaceDistance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    volScalarField* distPtr
    (
        new volScalarField
        (
            IOobject
            (
                "surfaceDistance",
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        )
    );

    regIOobject::store(distPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfaceDistance::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    doCells_ = dict.getOrDefault("calculateCells", true);

    geomPtr_.reset(nullptr);
    geomPtr_.reset
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                             // dummy name
                mesh_.time().constant(),           // directory
                "triSurface",                      // instance
                mesh_.time(),                      // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict.subDict("geometry"),
            true                // allow single-region shortcut
        )
    );

    return true;
}


bool Foam::functionObjects::surfaceDistance::execute()
{
    volScalarField& distance = mesh_.lookupObjectRef<volScalarField>
    (
        "surfaceDistance"
    );

    volScalarField::Boundary& bfld = distance.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (!polyPatch::constraintType(bfld[patchi].patch().type()))
        {
            const pointField& fc = mesh_.C().boundaryField()[patchi];

            labelList surfaces;
            List<pointIndexHit> nearestInfo;
            geomPtr_().findNearest
            (
                fc,
                scalarField(fc.size(), GREAT),
                surfaces,
                nearestInfo
            );

            scalarField dist(fc.size());
            forAll(nearestInfo, i)
            {
                dist[i] = nearestInfo[i].hitPoint().dist(fc[i]);
            }
            bfld[patchi] == dist;
        }
    }

    if (doCells_)
    {
        const pointField& cc = mesh_.C().internalField();

        labelList surfaces;
        List<pointIndexHit> nearestInfo;
        geomPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            surfaces,
            nearestInfo
        );

        forAll(nearestInfo, celli)
        {
            distance[celli] = nearestInfo[celli].hitPoint().dist(cc[celli]);
        }
    }
    distance.correctBoundaryConditions();

    return true;
}


bool Foam::functionObjects::surfaceDistance::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing distance-to-surface field" << endl;

    const volScalarField& distance =
        mesh_.lookupObject<volScalarField>("surfaceDistance");

//    volScalarField::Boundary& bfld = distance.boundaryFieldRef();
//    forAll(bfld, patchi)
//    {
//        if (!polyPatch::constraintType(bfld[patchi].patch().type()))
//        {
//            const pointField& fc = mesh_.C().boundaryField()[patchi];
//
//            labelList surfaces;
//            List<pointIndexHit> nearestInfo;
//            geomPtr_().findNearest
//            (
//                fc,
//                scalarField(fc.size(), GREAT),
//                surfaces,
//                nearestInfo
//            );
//
//            scalarField dist(fc.size());
//            forAll(nearestInfo, i)
//            {
//                dist[i] = nearestInfo[i].hitPoint().dist(fc[i]);
//            }
//            bfld[patchi] == dist;
//        }
//    }
//
//    if (doCells_)
//    {
//        const pointField& cc = mesh_.C().internalField();
//
//        labelList surfaces;
//        List<pointIndexHit> nearestInfo;
//        geomPtr_().findNearest
//        (
//            cc,
//            scalarField(cc.size(), GREAT),
//            surfaces,
//            nearestInfo
//        );
//
//        forAll(nearestInfo, celli)
//        {
//            distance[celli] = nearestInfo[celli].hitPoint().dist(cc[celli]);
//        }
//    }
//    distance.correctBoundaryConditions();
    distance.write();

    return true;
}


// ************************************************************************* //
