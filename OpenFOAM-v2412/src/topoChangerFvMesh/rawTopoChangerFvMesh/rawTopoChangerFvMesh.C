/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "rawTopoChangerFvMesh.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "linear.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rawTopoChangerFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        rawTopoChangerFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        rawTopoChangerFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::rawTopoChangerFvMesh::rawTopoChangerFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    topoChangerFvMesh(io, doInit)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rawTopoChangerFvMesh::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    Info<< "rawTopoChangerFvMesh : Checking for topology changes..."
        << endl;

    // Mesh not moved/changed yet
    moving(false);
    topoChanging(false);

    // Do any topology changes. Sets topoChanging (through polyTopoChange)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true);

    const bool hasChanged = bool(topoChangeMap);

    if (hasChanged)
    {
        Info<< "rawTopoChangerFvMesh : Done topology changes..."
            << endl;

        // Temporary: fix fields on patch faces created out of nothing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Two situations:
        // - internal faces inflated out of nothing
        // - patch faces created out of previously internal faces

        // Is face mapped in any way?
        bitSet mappedFace(nFaces());

        const label nOldInternal = topoChangeMap().oldPatchStarts()[0];

        const labelList& faceMap = topoChangeMap().faceMap();
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (faceMap[facei] >= 0)
            {
                mappedFace.set(facei);
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (faceMap[facei] >= 0 && faceMap[facei] >= nOldInternal)
            {
                mappedFace.set(facei);
            }
        }

        for (const auto& map : topoChangeMap().facesFromFacesMap())
        {
            mappedFace.set(map.index());
        }

        for (const auto& map : topoChangeMap().facesFromEdgesMap())
        {
            mappedFace.set(map.index());
        }

        for (const auto& map : topoChangeMap().facesFromPointsMap())
        {
            mappedFace.set(map.index());
        }

        // Set unmapped faces to zero
        Info<< "rawTopoChangerFvMesh : zeroing unmapped boundary values." << nl;

        zeroUnmappedValues<scalar, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<vector, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<sphericalTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<symmTensor, fvPatchField, volMesh>(mappedFace);
        zeroUnmappedValues<tensor, fvPatchField, volMesh>(mappedFace);

        // Special handling for phi: set unmapped faces to recreated phi
        Info<< "rawTopoChangerFvMesh :"
            << " recreating phi for unmapped boundary values." << endl;

        const auto& U = lookupObject<volVectorField>("U");
        auto& phi = lookupObjectRef<surfaceScalarField>("phi");

        setUnmappedValues
        (
            phi,
            mappedFace,
            (linearInterpolate(U) & Sf())()
        );


        if (topoChangeMap().hasMotionPoints())
        {
            pointField newPoints = topoChangeMap().preMotionPoints();

            // Give the meshModifiers opportunity to modify points
            Info<< "rawTopoChangerFvMesh :"
                << " calling modifyMotionPoints." << endl;
            topoChanger_.modifyMotionPoints(newPoints);

            // Actually move points
            Info<< "rawTopoChangerFvMesh :"
                << " calling movePoints." << endl;

            movePoints(newPoints);
        }
    }
    else
    {
        //Pout<< "rawTopoChangerFvMesh :"
        //    << " no topology changes..." << endl;
    }

    return hasChanged;
}


// ************************************************************************* //
