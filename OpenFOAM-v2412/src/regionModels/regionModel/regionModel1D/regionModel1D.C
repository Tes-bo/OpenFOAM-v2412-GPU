/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "regionModel1D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModel1D, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionModel1D::constructMeshObjects()
{
    nMagSfPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "nMagSf",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            regionMesh(),
            dimensionedScalar(dimArea, Zero)
        )
    );
}


void Foam::regionModels::regionModel1D::initialise()
{
    if (debug)
    {
        Pout<< "regionModel1D::initialise()" << endl;
    }

    // Calculate boundaryFaceFaces and boundaryFaceCells

    DynamicList<label> faceIDs;
    DynamicList<label> cellIDs;

    label localPyrolysisFacei = 0;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        localPyrolysisFacei += ppCoupled.size();
    }

    boundaryFaceOppositeFace_.setSize(localPyrolysisFacei);
    boundaryFaceFaces_.setSize(localPyrolysisFacei);
    boundaryFaceCells_.setSize(localPyrolysisFacei);

    localPyrolysisFacei = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        forAll(ppCoupled, localFacei)
        {
            label facei = ppCoupled.start() + localFacei;
            label celli = -1;
            label nCells = 0;
            do
            {
                label ownCelli = regionMesh().faceOwner()[facei];
                if (ownCelli != celli)
                {
                    celli = ownCelli;
                }
                else
                {
                    celli = regionMesh().faceNeighbour()[facei];
                }
                nCells++;
                cellIDs.append(celli);
                const cell& cFaces = regionMesh().cells()[celli];
                faceIDs.append(facei);
                label face0 =
                    cFaces.opposingFaceLabel(facei, regionMesh().faces());
                facei = face0;
            } while (regionMesh().isInternalFace(facei));

            boundaryFaceOppositeFace_[localPyrolysisFacei] = facei;
            //faceIDs.pop_back(); //remove boundary face.

            boundaryFaceFaces_[localPyrolysisFacei].transfer(faceIDs);
            boundaryFaceCells_[localPyrolysisFacei].transfer(cellIDs);

            localPyrolysisFacei++;
            nLayers_ = nCells;
        }
    }
    faceIDs.clear();
    cellIDs.clear();

    surfaceScalarField& nMagSf = nMagSfPtr_();

    surfaceScalarField::Boundary& nMagSfBf = nMagSf.boundaryFieldRef();

    localPyrolysisFacei = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        const vectorField& pNormals = ppCoupled.faceNormals();

        nMagSfBf[patchi] = regionMesh().Sf().boundaryField()[patchi] & pNormals;

        forAll(pNormals, localFacei)
        {
            const vector n = pNormals[localFacei];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFacei++];
            forAll(faces, facei)
            {
                // facei = 0 is on boundary
                if (facei > 0)
                {
                    const label faceID = faces[facei];
                    nMagSf[faceID] = regionMesh().Sf()[faceID] & n;
                }
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::regionModel1D::read()
{
    if (regionModel::read())
    {
        return true;
    }

    return false;
}


bool Foam::regionModels::regionModel1D::read(const dictionary& dict)
{
    if (regionModel::read(dict))
    {
        moveMesh_.readIfPresent("moveMesh", coeffs_);

        return true;
    }

    return false;
}


Foam::tmp<Foam::labelField> Foam::regionModels::regionModel1D::moveMesh
(
    const scalarList& deltaV,
    const scalar minDelta
)
{
    auto tcellMoveMap = tmp<labelField>::New(regionMesh().nCells(), Zero);
    auto& cellMoveMap = tcellMoveMap.ref();

    if (!moveMesh_)
    {
        return cellMoveMap;
    }

    pointField oldPoints = regionMesh().points();
    pointField newPoints = oldPoints;

    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();

    label totalFacei = 0;
    forAll(intCoupledPatchIDs_, localPatchi)
    {
        label patchi = intCoupledPatchIDs_[localPatchi];
        const polyPatch& pp = bm[patchi];

        forAll(pp, patchFacei)
        {
            const labelList& faces = boundaryFaceFaces_[totalFacei];
            const labelList& cells = boundaryFaceCells_[totalFacei];
            const label oFace = boundaryFaceOppositeFace_[totalFacei];

            const vector n = pp.faceNormals()[patchFacei];
            const vector sf = pp.faceAreas()[patchFacei];

            List<point> oldCf(faces.size() + 1, Zero);
            List<bool> frozen(faces.size(), false);

            forAll(faces, i)
            {
                oldCf[i] = regionMesh().faceCentres()[faces[i]];
            }

            oldCf[faces.size()] = regionMesh().faceCentres()[oFace];

            forAll(faces, i)
            {
                const label celli = cells[i];

                if (mag(oldCf[i + 1] - oldCf[i]) < minDelta)
                {
                    frozen[i] = true;
                    cellMoveMap[celli] = 1;
                }
            }

            vectorField newDelta(cells.size() + 1, Zero);

            label j = 0;
            forAll(cells, i)
            {
                const label celli = cells[i];
                newDelta[j+1] = (deltaV[celli]/mag(sf))*n + newDelta[j];
                j++;
            }

            // Move the back face first
            const face of = regionMesh().faces()[oFace];
            {
                scalar omagV = mag(newDelta[newDelta.size()-1]);

                if (!frozen[cells.size()-1] && (omagV > ROOTVSMALL))
                {
                    forAll(of, pti)
                    {
                        const label pointi = of[pti];
                        newPoints[pointi] =
                            oldPoints[pointi] - newDelta[newDelta.size()-1];
                    }
                }
            }
            // Do internal faces
            for (label i=0; i < faces.size(); i++)
            {
                const label facei = faces[i];
                const face f = regionMesh().faces()[facei];

                scalar magV = mag(newDelta[i]);
                if (!frozen[i] && magV > 0)
                {
                    forAll(f, pti)
                    {
                        const label pointi = f[pti];
                        newPoints[pointi] = oldPoints[pointi] - newDelta[i];
                    }
                }
            }

            totalFacei++;
        }
    }

    // Move points
    regionMesh().movePoints(newPoints);

    return tcellMoveMap;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel(mesh, regionType),
    boundaryFaceFaces_(),
    boundaryFaceCells_(),
    boundaryFaceOppositeFace_(),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(false)
{}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(true)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();
        if (readFields)
        {
            moveMesh_.readIfPresent("moveMesh", coeffs_);
        }
    }
}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, dict, readFields),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(false)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();
        if (readFields)
        {
            moveMesh_.readIfPresent("moveMesh", coeffs_);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::~regionModel1D()
{}


// ************************************************************************* //
