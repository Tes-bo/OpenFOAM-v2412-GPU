/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "fvGeometryScheme.H"
#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvGeometryScheme, 0);

    defineRunTimeSelectionTable(fvGeometryScheme, dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::fvGeometryScheme::setMeshPhi() const
{
    if (!mesh_.moving())
    {
        return false;
    }

    const pointField& oldPoints = mesh_.oldPoints();
    const pointField& currPoints = mesh_.points();

    if (oldPoints.size() != currPoints.size())
    {
        FatalErrorInFunction
            << "Old and current points sizes must be the same. "
            << "Old points:" << oldPoints.size()
            << " Current points:" << currPoints.size()
            << abort(FatalError);
    }

    const faceList& faces = mesh_.faces();
    const scalar rdt = 1.0/mesh_.time().deltaTValue();

    auto tmeshPhi(const_cast<fvMesh&>(mesh_).setPhi());
    if (tmeshPhi)
    {
        auto& meshPhi = tmeshPhi.ref();
        auto& meshPhii = meshPhi.primitiveFieldRef();
        forAll(meshPhii, facei)
        {
            const face& f = faces[facei];
            meshPhii[facei] = f.sweptVol(oldPoints, currPoints)*rdt;
        }

        auto& meshPhiBf = meshPhi.boundaryFieldRef();
        for (auto& meshPhip : meshPhiBf)
        {
            if (!meshPhip.size())
            {
                // Empty patches
                continue;
            }

            const auto& pp = meshPhip.patch().patch();

            forAll(pp, facei)
            {
                const face& f = pp[facei];
                meshPhip[facei] = f.sweptVol(oldPoints, currPoints)*rdt;
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvGeometryScheme> Foam::fvGeometryScheme::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& defaultScheme
)
{
    const entry* eptr = dict.findEntry("method", keyType::LITERAL);
    const word schemeName
    (
        eptr
      ? word(eptr->stream())
      : dict.getOrDefault<word>("type", defaultScheme)
    );

    DebugInFunction << "Geometry scheme = " << schemeName << endl;

    auto* ctorPtr = dictConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "fvGeometryScheme",
            schemeName,
            *dictConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvGeometryScheme::movePoints()
{
    if (mesh_.moving())
    {
        // Set the mesh motion fluxes
        setMeshPhi();

        // Clear out old geometry
        // Note: this recreates the old primitiveMesh::movePoints behaviour
        const_cast<fvMesh&>(mesh_).primitiveMesh::clearGeom();
    }
}


void Foam::fvGeometryScheme::updateMesh(const mapPolyMesh& mpm)
{}


// ************************************************************************* //
