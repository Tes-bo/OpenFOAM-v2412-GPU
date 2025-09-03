/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "viewFactorModel.H"
#include "raySearchEngine.H"
#include "OBJstream.H"
#include "volFields.H"
#include "IOmapDistribute.H"
#include "scalarListIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace VF
{
    defineTypeNameAndDebug(viewFactorModel, 0);
    defineRunTimeSelectionTable(viewFactorModel, mesh);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::VF::viewFactorModel::writeRays
(
    const fileName& fName,
    const pointField& compactCf,
    const labelListList& visibleFaceFaces
)
{
    OBJstream str(fName);

    Pout<< "Writing rays to " << str.name() << endl;

    forAll(visibleFaceFaces, facei)
    {
        const labelList& visibleSlots = visibleFaceFaces[facei];

        for (const label sloti : visibleSlots)
        {
            str.write(linePointRef(compactCf[facei], compactCf[sloti]));
        }
    }

    str.flush();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VF::viewFactorModel::viewFactorModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    searchEnginePtr_(raySearchEngine::New(mesh, dict)),
    writeViewFactors_(dict.get<bool>("writeViewFactors")),
    writeRays_(dict.getOrDefault<bool>("writeRays", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VF::viewFactorModel::~viewFactorModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VF::viewFactorModel::calculate()
{
    const auto& searchEngine = *searchEnginePtr_;

    labelListList visibleFaceFaces;
    searchEngine.correct(visibleFaceFaces);
    const auto& map = searchEngine.map();

    // Construct data in compact addressing
    // Compact addressing has all local data, followed by remote contributions

    pointField compactCf;
    vectorField compactSf;
    List<List<vector>> compactFineSf;
    List<List<point>> compactFineCf;
    DynamicList<List<point>> compactPoints;
    DynamicList<label> compactPatchId;

    searchEngine.compactAddressing
    (
        map,
        compactCf,
        compactSf,
        compactFineSf,
        compactFineCf,
        compactPoints,
        compactPatchId
    );

    if (writeRays_)
    {
        // Write all rays between visible faces using .OBJ format

        writeRays
        (
            mesh_.time().path()/"allVisibleFaces.obj",
            compactCf,
            visibleFaceFaces
        );
    }

    (void)mesh_.time().cpuTimeIncrement();

    Info<< "\nCalculating view factors" << endl;

    scalarListIOList Fij
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::NO_READ
        ),
        calculate
        (
            visibleFaceFaces,
            compactCf,
            compactSf,
            compactFineSf,
            compactFineCf,
            compactPoints,
            compactPatchId
        )
    );

    const label totalPatches = mesh_.boundaryMesh().nNonProcessor();

    // Matrix sum in j(Fij) for each i
    // Note: if enclosure sum should = 1
    scalarSquareMatrix viewFactorPatch(totalPatches, Zero);

    forAll(visibleFaceFaces, startFacei)
    {
        const scalar magAi = mag(compactSf[startFacei]);

        const labelList& visibleSlots = visibleFaceFaces[startFacei];
        const label patchi = compactPatchId[startFacei];

        forAll(visibleSlots, visSloti)
        {
            const label sloti = visibleSlots[visSloti];
            const label patchj = compactPatchId[sloti];

            viewFactorPatch[patchi][patchj] += Fij[startFacei][visSloti]*magAi;
        }
    }

    reduce(viewFactorPatch, sumOp<scalarSquareMatrix>());

    const scalarList patchArea = searchEngine.patchAreas();


    if (UPstream::master())
    {
        Info<< "\nPatch view factor contributions:" << nl << endl;

        scalar vfSum = 0;

        const auto& patchIDs = searchEngine.patchIDs();
        const auto& patches = mesh_.boundaryMesh();

        for (const label patchi : patchIDs)
        {
            Info<< "    Patch " <<  patchi << ": " << patches[patchi].name()
                << endl;

            scalar vfPatch = 0;
            for (const label patchj: patchIDs)
            {
                scalar vf = viewFactorPatch[patchi][patchj]/patchArea[patchi];
                vfPatch += vf;

                Info<< "    F" << patchi << patchj << ": " << vf << endl;
            }

            Info<< "    Sum: " << vfPatch << nl << endl;

            vfSum += vfPatch;
        }

        Info<< "Sum(all patches) = " << vfSum << endl;
    }

    // Write view factors matrix in list-list form
    Info<< "\nWriting view factor matrix" << endl;
    Fij.write();

    if (writeViewFactors_)
    {
        Info<< "\nWriting view factor field" << endl;

        auto tviewFactorField =
            volScalarField::New
            (
                "viewFactorField",
                mesh_,
                dimensionedScalar(dimless, Zero)
            );

        auto& viewFactorField = tviewFactorField.ref();

        searchEngine.interpolate(viewFactorField, Fij);

        viewFactorField.write();
    }


    // Create globalFaceFaces needed to insert view factors
    // in F to the global matrix Fmatrix
    {
        labelListIOList IOglobalFaceFaces
        (
            IOobject
            (
                "globalFaceFaces",
                mesh_.facesInstance(),
                mesh_,
                IOobject::NO_READ
            ),
            visibleFaceFaces.size()
        );

        forAll(IOglobalFaceFaces, facei)
        {
            IOglobalFaceFaces[facei] = renumber
            (
                searchEngine.compactToGlobal(),
                visibleFaceFaces[facei]
            );
        }

        IOglobalFaceFaces.write();
    }

    // Write parallel map
    {
        IOmapDistribute IOmapDist
        (
            IOobject
            (
                "mapDist",
                mesh_.facesInstance(),
                mesh_,
                IOobject::NO_READ
            ),
            std::move(map)
        );

        IOmapDist.write();
    }
}


// ************************************************************************* //
