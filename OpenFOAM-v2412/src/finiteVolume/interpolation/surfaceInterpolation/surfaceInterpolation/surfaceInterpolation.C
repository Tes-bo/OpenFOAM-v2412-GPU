/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

Description
    Cell to face interpolation scheme. Included in fvMesh.

\*---------------------------------------------------------------------------*/

#include "surfaceInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatch.H"
#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceInterpolation, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::surfaceInterpolation::clearOut()
{
    // TBD: potential to apply partial clear out only?
    // Move to fvGeometryScheme?
    weights_.clear();
    deltaCoeffs_.clear();
    nonOrthDeltaCoeffs_.clear();
    nonOrthCorrectionVectors_.clear();
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::surfaceInterpolation(const fvMesh& fvm)
:
    mesh_(fvm),
    geometryPtr_(nullptr),
    weights_(nullptr),
    deltaCoeffs_(nullptr),
    nonOrthDeltaCoeffs_(nullptr),
    nonOrthCorrectionVectors_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::~surfaceInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvGeometryScheme& Foam::surfaceInterpolation::geometry() const
{
    if (!geometryPtr_)
    {
        geometryPtr_ = fvGeometryScheme::New
        (
            mesh_,
            mesh_.schemesDict().subOrEmptyDict("geometry"),
            basicFvGeometryScheme::typeName
        );
    }

    return geometryPtr_();
}


void Foam::surfaceInterpolation::geometry(tmp<fvGeometryScheme>& schemePtr)
{
    geometryPtr_ = schemePtr;
}


const Foam::surfaceScalarField& Foam::surfaceInterpolation::weights() const
{
    if (!weights_)
    {
        weights_.reset(geometry().weights().ptr());
    }

    return weights_();
}


const Foam::surfaceScalarField& Foam::surfaceInterpolation::deltaCoeffs() const
{
    if (!deltaCoeffs_)
    {
        deltaCoeffs_.reset(geometry().deltaCoeffs().ptr());
    }

    return deltaCoeffs_();
}


const Foam::surfaceScalarField&
Foam::surfaceInterpolation::nonOrthDeltaCoeffs() const
{
    if (!nonOrthDeltaCoeffs_)
    {
        nonOrthDeltaCoeffs_.reset(geometry().nonOrthDeltaCoeffs().ptr());
    }

    return nonOrthDeltaCoeffs_();
}


const Foam::surfaceVectorField&
Foam::surfaceInterpolation::nonOrthCorrectionVectors() const
{
    if (!nonOrthCorrectionVectors_)
    {
        nonOrthCorrectionVectors_.reset
        (
            geometry().nonOrthCorrectionVectors().ptr()
        );
    }

    return nonOrthCorrectionVectors_();
}


bool Foam::surfaceInterpolation::movePoints()
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::movePoints() : "
            << "Updating geometric properties using the fvGeometryScheme"
            << endl;
    }

    // Do any primitive geometry calculation
    const_cast<fvGeometryScheme&>(geometry()).movePoints();

    clearOut();

    return true;
}


void Foam::surfaceInterpolation::updateGeom()
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::updateGeom() : "
            << "Updating geometric properties" << endl;
    }

    const_cast<fvGeometryScheme&>(geometry()).movePoints();

    clearOut();
}


void Foam::surfaceInterpolation::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::updateMesh() : "
            << "Updating geometric properties" << endl;
    }

    const_cast<fvGeometryScheme&>(geometry()).updateMesh(mpm);

    clearOut();
}


// ************************************************************************* //
