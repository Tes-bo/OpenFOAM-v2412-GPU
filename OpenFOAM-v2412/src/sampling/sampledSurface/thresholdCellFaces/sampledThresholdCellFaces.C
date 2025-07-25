/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "sampledThresholdCellFaces.H"
#include "thresholdCellFaces.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledThresholdCellFaces, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledThresholdCellFaces,
        word,
        thresholdCellFaces
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledThresholdCellFaces::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // no update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Use volField from database, or try to read it in

    const auto* cellFldPtr = fvm.findObject<volScalarField>(fieldName_);

    if (debug)
    {
        if (cellFldPtr)
        {
            InfoInFunction << "Lookup " << fieldName_ << endl;
        }
        else
        {
            InfoInFunction
                << "Reading " << fieldName_
                << " from time " << fvm.time().timeName()
                << endl;
        }
    }

    // For holding the volScalarField read in.
    autoPtr<volScalarField> fieldReadPtr;

    if (!cellFldPtr)
    {
        // Bit of a hack. Read field and store.

        fieldReadPtr = autoPtr<volScalarField>::New
        (
            IOobject
            (
                fieldName_,
                fvm.time().timeName(),
                fvm,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvm
        );
    }

    const volScalarField& cellFld =
        (fieldReadPtr ? *fieldReadPtr : *cellFldPtr);

    Mesh& mySurface = const_cast<sampledThresholdCellFaces&>(*this);

    thresholdCellFaces surf
    (
        fvm,
        cellFld.primitiveField(),
        lowerThreshold_,
        upperThreshold_,
        triangulate_
    );

    mySurface.transfer(static_cast<Mesh&>(surf));
    meshCells_.transfer(surf.meshCells());

    // Clear derived data
    sampledSurface::clearGeom();

    if (debug)
    {
        Pout<< "sampledThresholdCellFaces::updateGeometry() : constructed"
            << nl
            << "    field         : " << fieldName_ << nl
            << "    lowerLimit    : " << lowerThreshold_ << nl
            << "    upperLimit    : " << upperThreshold_ << nl
            << "    point         : " << points().size() << nl
            << "    faces         : " << faces().size() << nl
            << "    cut cells     : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledThresholdCellFaces::sampledThresholdCellFaces
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    fieldName_(dict.get<word>("field")),
    lowerThreshold_(dict.getOrDefault<scalar>("lowerLimit", -VGREAT)),
    upperThreshold_(dict.getOrDefault<scalar>("upperLimit", VGREAT)),
    triangulate_(dict.getOrDefault("triangulate", false)),
    prevTimeIndex_(-1),
    meshCells_()
{
    if (!dict.found("lowerLimit") && !dict.found("upperLimit"))
    {
        FatalErrorInFunction
            << "require at least one of 'lowerLimit' or 'upperLimit'" << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledThresholdCellFaces::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledThresholdCellFaces::expire()
{
    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledThresholdCellFaces::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField> Foam::sampledThresholdCellFaces::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledThresholdCellFaces::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledThresholdCellFaces::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledThresholdCellFaces::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledThresholdCellFaces::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledThresholdCellFaces::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledThresholdCellFaces::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledThresholdCellFaces::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledThresholdCellFaces::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledThresholdCellFaces::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledThresholdCellFaces::print(Ostream& os, int level) const
{
    os  << "sampledThresholdCellFaces: " << name() << " :"
        << " field:" << fieldName_
        << " lowerLimit:" << lowerThreshold_
        << " upperLimit:" << upperThreshold_;

    // Possibly no geom yet...
    // if (level)
    // {
    //     os  << "  faces:" << faces().size()
    //         << "  points:" << points().size();
    // }
}


// ************************************************************************* //
