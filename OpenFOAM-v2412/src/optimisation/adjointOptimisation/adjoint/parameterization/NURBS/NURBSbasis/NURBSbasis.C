/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "NURBSbasis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(NURBSbasis, 1);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void NURBSbasis::computeKnots()
{
    // Sanity check
    if (basisDegree_ >(nCPs_ - 1))
    {
        FatalErrorInFunction
            << "B - splines basis degree can be at most equal to the "
            << "number of control points minus 1"
            << exit(FatalError);
    }

    // First zero knots
    for (label ik = 0; ik < basisDegree_ + 1; ik++)
    {
        knots_ = scalar(0);
    }

    // Intermediate knots
    label firstCPIndex(basisDegree_ + 1);
    label lastCPIndex(knots_.size() - basisDegree_ - 1);
    label size(knots_.size() - 2*basisDegree_ - 2);

    for (label ik = 0; ik < size; ik++)
    {
        knots_[ik + firstCPIndex] = scalar(ik + 1)/scalar(size + 1);
    }

    // Last unity knots
    for (label ik = 0; ik < basisDegree_ + 1; ik++)
    {
        knots_[ik + lastCPIndex] = scalar(1);
    }

    DebugInfo
        << "Using knots " << knots_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NURBSbasis::NURBSbasis
(
    const label nCPs,
    const label degree,
    const scalarField& knots
)
:
    nCPs_(nCPs),
    basisDegree_(degree),
    knots_(knots)
{}


NURBSbasis::NURBSbasis
(
    const label nCPs,
    const label degree
)
:
    nCPs_(nCPs),
    basisDegree_(degree),
    knots_((nCPs_ + basisDegree_ + 1), Zero)
{
    computeKnots();
}


NURBSbasis::NURBSbasis
(
    const dictionary& dict
)
:
    nCPs_(dict.get<label>("nCPs")),
    basisDegree_(dict.get<label>("basisDegree")),
    knots_((nCPs_ + basisDegree_ + 1), Zero)
{
    computeKnots();
}


NURBSbasis::NURBSbasis
(
    const NURBSbasis& basis
)
:
    nCPs_(basis.nCPs_),
    basisDegree_(basis.basisDegree_),
    knots_(basis.knots_)
{
    DebugInfo
        << "Copied basis function" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar NURBSbasis::basisValue
(
    const label iCP,
    const label degree,
    const scalar u
) const
{
    scalar value(0);

    if (checkRange(u, iCP, degree))
    {
        // If the degree reached zero, return step function
        if (degree == 0)
        {
            if ((u >= knots_[iCP]) && (u < knots_[iCP + 1]))
            {
                value = scalar(1);
            }
            else if ((u == 1) && (knots_[iCP + 1] == 1))
            {
                value = scalar(1);
            }
        }
        // Else, call recursively until reaching zero degree
        else
        {
            const scalar denom1(knots_[iCP + degree] - knots_[iCP]);
            const scalar denom2(knots_[iCP + degree + 1] - knots_[iCP + 1]);

            if (denom1 != 0)
            {
                value +=
                    (u - knots_[iCP])
                  * basisValue(iCP, degree - 1, u)
                  / denom1;
            }
            if (denom2 != 0)
            {
                value +=
                    (knots_[iCP + degree + 1] - u)
                  * basisValue(iCP + 1, degree - 1, u)
                  / denom2;
            }
        }
    }

    return value;
}


scalar NURBSbasis::basisDerivativeU
(
    const label iCP,
    const label degree,
    const scalar u
) const
{
    // Zero  basis function has a zero derivative,
    // irrespective of the knot span: ignore that case
    // - else, call recursively until reaching zero degree
    scalar derivative(0);

    if ((degree != 0) && checkRange(u, iCP, degree))
    {
        const scalar denom1(knots_[iCP + degree] - knots_[iCP]);
        const scalar denom2(knots_[iCP + degree + 1] - knots_[iCP + 1]);

        if (denom1 != 0)
        {
            derivative +=
                (
                    (u - knots_[iCP])
                  * basisDerivativeU(iCP, degree - 1, u)
                  + basisValue(iCP, degree - 1, u)
                ) / denom1;
        }
        if (denom2 != 0)
        {
            derivative +=
                (
                    (knots_[iCP + degree + 1] - u)
                  * basisDerivativeU(iCP + 1, degree - 1, u)
                  - basisValue(iCP + 1, degree - 1, u)
                ) / denom2;
        }
    }

    return derivative;
}


scalar NURBSbasis::basisDerivativeUU
(
    const label iCP,
    const label degree,
    const scalar u
) const
{
    // Zero  basis function has a zero derivative,
    // irrespective of the knot span: ignore that case
    // - else, call recursively until reaching zero degree
    scalar derivative(0);

    if ((degree != 0) && checkRange(u, iCP, degree))
    {
        scalar denom1 = (knots_[iCP + degree] - knots_[iCP]);
        scalar denom2 = (knots_[iCP + degree + 1] - knots_[iCP + 1]);

        if (denom1 != 0)
        {
            derivative +=
                (
                    (u - knots_[iCP])
                  * basisDerivativeUU(iCP, degree - 1, u)
                  + 2*basisDerivativeU(iCP, degree - 1, u)
                ) / denom1;
        }

        if (denom2 != 0)
        {
            derivative +=
                (
                    (knots_[iCP + degree + 1] - u)
                  * basisDerivativeUU(iCP + 1, degree - 1, u)
                  - 2*basisDerivativeU(iCP + 1, degree - 1, u)
                ) / denom2;
        }
    }

    return derivative;
}


bool NURBSbasis::checkRange
(
    const scalar u,
    const label CPI,
    const label degree
) const
{
    // Find what section of the curve given u is in.
    const scalar lowerBound(knots_[CPI]);
    const scalar upperBound(knots_[CPI + degree + 1]);

    // Check in-range requirements.
    if
    (
       ((u == scalar(1)) && (lowerBound <= u) && (u <= upperBound))
    || ((u != scalar(1)) && (lowerBound <= u) && (u < upperBound))
    )
    {
        return true;
    }

    return false;
}


label NURBSbasis::insertKnot
(
    const scalar uBar
)
{
    // Find the index of insertion, accounting for uBar=1.
    scalarList newKnots((knots_.size()+1), Zero);
    label kInsert(-1);

    for (label kI = 0; kI < (knots_.size()-1); kI++)
    {
        if (knots_[kI + 1] > uBar)
        {
            kInsert = kI;
            break;
        }
    }

    if (kInsert == -1)
    {
        kInsert = knots_.size()-1;
    }

    // Update data.
    for (label kI = 0; kI<(kInsert + 1); kI++)
    {
        newKnots[kI] = knots_[kI];
    }

    newKnots[kInsert + 1] = uBar;

    for (label kI= (kInsert + 2); kI < newKnots.size(); kI++)
    {
        newKnots[kI] = knots_[kI - 1];
    }

    // Add the new CP info and return the insertion point.
    knots_ = newKnots;
    nCPs_++;

    return kInsert;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
