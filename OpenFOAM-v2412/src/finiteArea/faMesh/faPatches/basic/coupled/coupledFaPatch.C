/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "coupledFaPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(coupledFaPatch);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::coupledFaPatch::calcTransformTensors
(
    const vector& Cf,
    const vector& Cr,
    const vector& nf,
    const vector& nr
) const
{
    if (mag(nf & nr) < 1 - SMALL)
    {
        separation_.clear();
        forwardT_.resize(1);
        reverseT_.resize(1);

        forwardT_[0] = rotationTensor(-nr, nf);
        reverseT_[0] = rotationTensor(nf, -nr);
    }
    else
    {
        forwardT_.clear();
        reverseT_.clear();

        vector separation = (nf & (Cr - Cf))*nf;

        if (mag(separation) > SMALL)
        {
            separation_.resize(1);
            separation_[0] = separation;
        }
        else
        {
            separation_.clear();
        }
    }
}


void Foam::coupledFaPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nr
) const
{
    if (sum(mag(nf & nr)) < Cf.size() - SMALL)
    {
        separation_.clear();

        forwardT_.resize_nocopy(size());
        reverseT_.resize_nocopy(size());

        forAll(forwardT_, facei)
        {
            forwardT_[facei] = rotationTensor(-nr[facei], nf[facei]);
            reverseT_[facei] = rotationTensor(nf[facei], -nr[facei]);
        }

        if (sum(mag(forwardT_ - forwardT_[0])) < SMALL)
        {
            forwardT_.resize(1);
            reverseT_.resize(1);
        }
    }
    else
    {
        forwardT_.clear();
        reverseT_.clear();

        separation_ = (nf&(Cr - Cf))*nf;

        if (sum(mag(separation_)) < SMALL)
        {
            separation_.clear();
        }
        else if (sum(mag(separation_ - separation_[0])) < SMALL)
        {
            separation_.resize(1);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::coupledFaPatch::delta() const
{
    return (edgeCentres() - edgeFaceCentres());
}


// ************************************************************************* //
