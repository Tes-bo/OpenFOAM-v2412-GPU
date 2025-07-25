/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "exactPatchDistMethod.H"
#include "distributedTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "DynamicField.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(exact, 0);
    addToRunTimeSelectionTable(patchDistMethod, exact, dictionary);
}
}

const Foam::distributedTriSurfaceMesh&
Foam::patchDistMethods::exact::patchSurface() const
{
    if (!patchSurfPtr_)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        Random rndGen(0);

        treeBoundBox localBb(mesh_.points());

        // Determine mesh bounding boxes:
        List<treeBoundBox> meshBb
        (
            1,
            localBb.extend(rndGen, 1E-3)
        );

        // Dummy bounds dictionary
        dictionary dict;
        dict.add("bounds", meshBb);
        dict.add
        (
            "distributionType",
            distributedTriSurfaceMesh::distributionTypeNames_
            [
                //distributedTriSurfaceMesh::FOLLOW         // use mesh bb
                //distributedTriSurfaceMesh::INDEPENDENT    // master-only
                distributedTriSurfaceMesh::DISTRIBUTED      // parallel decomp
            ]
        );
        dict.add("mergeDistance", 1e-6*localBb.mag());


        Info<< "Triangulating local patch faces" << nl << endl;

        labelList mapTriToGlobal;

        patchSurfPtr_.reset
        (
            new distributedTriSurfaceMesh
            (
                IOobject
                (
                    "wallSurface.stl",
                    mesh_.time().constant(),// directory
                    "triSurface",           // instance
                    mesh_.time(),           // registry
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                triSurfaceTools::triangulate
                (
                    pbm,
                    patchIDs_,
                    mapTriToGlobal
                ),
                dict
            )
        );

        // Do redistribution
        Info<< "Redistributing surface" << nl << endl;
        autoPtr<mapDistribute> faceMap;
        autoPtr<mapDistribute> pointMap;
        patchSurfPtr_().distribute
        (
            meshBb,
            false,           //keepNonMapped,
            faceMap,
            pointMap
        );
        faceMap.clear();
        pointMap.clear();
    }
    return patchSurfPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::exact::exact
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


Foam::patchDistMethods::exact::exact
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::exact::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::exact::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    const distributedTriSurfaceMesh& surf = patchSurface();

    List<pointIndexHit> info;
    surf.findNearest
    (
        mesh_.cellCentres(),
        scalarField(mesh_.nCells(), Foam::sqr(GREAT)),
        info
    );

    // Take over hits
    // label nHits = 0;
    forAll(info, celli)
    {
        if (info[celli].hit())
        {
            const point& cc = mesh_.cellCentres()[celli];
            y[celli] = info[celli].point().dist(cc);
            // ++nHits;
        }
        //else
        //{
        //    // Miss. Do what? Not possible with GREAT hopefully ...
        //}
    }
    y.correctBoundaryConditions();

    if (debug)
    {
        OBJstream str(mesh_.time().timePath()/"wallPoint.obj");
        Info<< type() << ": dumping nearest wall point to "
            << str.name() << endl;
        forAll(mesh_.cellCentres(), celli)
        {
            const point& cc = mesh_.cellCentres()[celli];
            str.writeLine(cc, info[celli].point());
        }
    }


    // Only calculate n if the field is defined
    if (notNull(n))
    {
        surf.getNormal(info, n.primitiveFieldRef());
        n.correctBoundaryConditions();
    }

    return true;
}


// ************************************************************************* //
