/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) Wikki Ltd
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "pointVolInterpolation.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "volFields.H"
#include "emptyFvPatch.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointVolInterpolation, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointVolInterpolation::makeWeights() const
{
    if (volWeightsPtr_)
    {
        FatalErrorInFunction
            << "weighting factors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "pointVolInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const pointField& points = vMesh().points();
    const labelListList& cellPoints = vMesh().cellPoints();
    const vectorField& cellCentres = vMesh().cellCentres();

    // Allocate storage for weighting factors
    volWeightsPtr_ = new FieldField<Field, scalar>(cellCentres.size());
    FieldField<Field, scalar>& weightingFactors = *volWeightsPtr_;

    forAll(weightingFactors, pointi)
    {
        weightingFactors.set
        (
            pointi,
            new scalarField(cellPoints[pointi].size())
        );
    }


    // Calculate inverse distances between points and cell centres
    // and store in weighting factor array
    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cellPointI)
        {
            weightingFactors[cellI][cellPointI] = 1.0/
                mag
                (
                    cellCentres[cellI] - points[curCellPoints[cellPointI]]
                );
        }
    }

    scalarField pointVolSumWeights(cellCentres.size(), Zero);

    // Calculate weighting factors using inverse distance weighting
    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cellPointI)
        {
            pointVolSumWeights[cellI] += weightingFactors[cellI][cellPointI];
        }
    }

    forAll(cellCentres, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cellPointI)
        {
            weightingFactors[cellI][cellPointI] /= pointVolSumWeights[cellI];
        }
    }

    if (debug)
    {
        Info<< "pointVolInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// Do what is necessary if the mesh has changed topology
void Foam::pointVolInterpolation::clearAddressing() const
{
    deleteDemandDrivenData(patchInterpolatorsPtr_);
}


// Do what is necessary if the mesh has moved
void Foam::pointVolInterpolation::clearGeom() const
{
    deleteDemandDrivenData(volWeightsPtr_);
}


const Foam::PtrList<Foam::primitivePatchInterpolation>&
Foam::pointVolInterpolation::patchInterpolators() const
{
    if (!patchInterpolatorsPtr_)
    {
        const fvBoundaryMesh& bdry = vMesh().boundary();

        patchInterpolatorsPtr_ =
            new PtrList<primitivePatchInterpolation>(bdry.size());

        forAll(bdry, patchI)
        {
            patchInterpolatorsPtr_->set
            (
                patchI,
                new primitivePatchInterpolation(bdry[patchI].patch())
            );
        }
    }

    return *patchInterpolatorsPtr_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::pointVolInterpolation::pointVolInterpolation
(
 const pointMesh& pm,
    const fvMesh& vm
)
:
    pointMesh_(pm),
    fvMesh_(vm),
    volWeightsPtr_(nullptr),
    patchInterpolatorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::pointVolInterpolation::~pointVolInterpolation()
{
    clearAddressing();
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return point weights
const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::pointVolInterpolation::volWeights() const
{
    // If weighting factor array not present then construct
    if (!volWeightsPtr_)
    {
        makeWeights();
    }

    return *volWeightsPtr_;
}


// Do what is necessary if the mesh has moved
void Foam::pointVolInterpolation::updateTopology()
{
    clearAddressing();
    clearGeom();
}


// Do what is necessary if the mesh has moved
bool Foam::pointVolInterpolation::movePoints()
{
    clearGeom();

    return true;
}


// ************************************************************************* //
