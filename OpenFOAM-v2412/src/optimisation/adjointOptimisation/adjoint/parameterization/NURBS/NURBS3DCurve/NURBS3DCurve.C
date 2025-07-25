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

#include "NURBS3DCurve.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label NURBS3DCurve::sgn(const scalar val) const
{
    return val>=scalar(0) ? 1 : -1;
}


scalar NURBS3DCurve::abs(const scalar val) const
{
    return (sgn(val) == 1) ? val : -val;
}


label NURBS3DCurve::mod(const label x, const label interval) const
{
    label ratio(x%interval);
    return ratio < 0 ? ratio + interval : ratio;
}


void NURBS3DCurve::setUniformU()
{
    const vectorField& curve(*this);
    label nPts(curve.size());

    forAll(curve, ptI)
    {
        u_[ptI] = scalar(ptI)/scalar(nPts - 1);
    }
}


bool NURBS3DCurve::bound
(
    scalar& u,
    const scalar minVal = 1e-7,
    const scalar maxVal = 0.999999
) const
{
    // lower value bounding
    if (u < scalar(0))
    {
        u = minVal;
        return true;
        //Info<< "Lower bound hit." << endl;
    }

    // upper value bounding
    if (u > scalar(1))
    {
        u = maxVal;
        return true;
        //Info<< "Upper bound hit." << endl;
    }

    return false;
}


void NURBS3DCurve::setEquidistantU
(
    scalarList& U,
    const label lenAcc = 25,
    const label maxIter = 10,
    const label spacingCorrInterval=-1,
    const scalar tolerance = 1.e-5
) const
{
    const label nPts(U.size());
    const scalar xLength(length() /(nPts - 1));
    const scalar uLength(scalar(1) / scalar(nPts - 1));

    U[0] = Zero;
    U[nPts - 1] = scalar(1);

    for (label ptI = 1; ptI<(nPts - 1); ptI++)
    {
        const scalar UPrev(U[ptI - 1]);
        scalar& UCurr(U[ptI]);
        scalar direc(scalar(1));
        scalar xDiff(scalar(0));
        scalar delta(scalar(0));
        bool overReached(false);

        UCurr = UPrev + uLength;

        // Find the starting U value to ensure target is within 1 uLength.
        while (true)
        {
            bool bounded(bound(UCurr));

            delta = length(UPrev, UCurr, lenAcc);
            xDiff = xLength - delta;

            // Found the point.
            if (abs(xDiff) < tolerance)
            {
                break;
            }
            else
            {
                direc = sgn(xDiff);

                if (bounded && (direc == 1))
                {
                    // uLength addition makes U exceed 1 so it becomes bounded.
                    // However, the desired x location still exceeds how far the
                    // bounded uLength can move (~e-5 error).
                    // Must force U to be u=0.999999.
                    overReached = true;
                    break;
                }
                else if (direc == scalar(1))
                {
                    UCurr += uLength;
                }
                else
                {
                    break;
                }
            }
        }

        if (!overReached)
        {
            label iter(0);

            while (iter < maxIter)
            {
                // Set the new search length to ensure convergence and next v.
                direc /= scalar(2);
                UCurr += direc * uLength;

                bound(UCurr);

                // Can employ an occasional tolerance check from beg of curve.
                if (
                        (spacingCorrInterval != -1)
                     && (mod(ptI, spacingCorrInterval) == 0)
                   )
                {
                    delta = length(scalar(0), UCurr, ptI*lenAcc);
                    xDiff = (ptI * xLength) - delta;
                }
                else
                {
                    delta = length(UPrev, UCurr, lenAcc);
                    xDiff = xLength - delta;
                }

                // Break if found point or find the new search direction.
                if (abs(xDiff) < tolerance)
                {
                    break;
                }
                else
                {
                    const scalar oldDirec(direc);
                    direc = sgn(xDiff) * abs(oldDirec);
                }

                iter++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NURBS3DCurve::NURBS3DCurve
(
    const NURBSbasis& basis,
    const List<vector>& CPs,
    const List<scalar>& weights,
    const scalarField& u,
    const label nPts,
    const word name
)
:
    vectorField(nPts, Zero),
    CPs_(CPs),
    weights_(weights),
    u_(u),
    name_(name),
    basis_(basis),

    givenInitNrm_(Zero),

    nrmOrientation_(ALIGNED)

{
    buildCurve();
}


NURBS3DCurve::NURBS3DCurve
(
    const NURBSbasis& basis,
    const List<vector>& CPs,
    const List<scalar>& weights,
    const label nPts,
    const word name
)
:
    vectorField(nPts, Zero),
    CPs_(CPs),
    weights_(weights),
    u_(nPts, Zero),
    name_(name),
    basis_(basis),

    givenInitNrm_(Zero),

    nrmOrientation_(ALIGNED)

{
    setUniformU();
    buildCurve();
}


NURBS3DCurve::NURBS3DCurve
(
    const NURBSbasis& basis,
    const List<vector>& CPs,
    const label nPts,
    const word name
)
:
    vectorField(nPts, Zero),
    CPs_(CPs),
    weights_(CPs.size(), scalar(1)),
    u_(nPts, Zero),
    name_(name),
    basis_(basis),

    givenInitNrm_(Zero),

    nrmOrientation_(1)

{
    setUniformU();
    buildCurve();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Set Functions * * * * * * * * * * * * * * //

void NURBS3DCurve::setNrm3DOrientation
(
    const vector& givenNrm,
    const vector& givenTan
)
{
    givenInitNrm_ = givenNrm;

    const vector tan(curveDerivativeU(Zero));
    vector curveNrm(tan ^ givenTan);

    if ((givenNrm & curveNrm) >= scalar(0))
    {
        nrmOrientation_ = ALIGNED;
    }
    else
    {
        nrmOrientation_ = OPPOSED;
    }

    Info<< "Initial nrmOrientation after comparison to NURBS u = 0 nrm : "
        << nrmOrientation_
        << endl;
}


void NURBS3DCurve::setNrm2DOrientation
(
    const vector& givenNrm,
    const scalar zVal
)
{
    givenInitNrm_ = givenNrm;

    const vector tan(curveDerivativeU(Zero));
    vector curveNrm(Zero);

    curveNrm.x() = -tan.y();
    curveNrm.y() =  tan.x();
    curveNrm.z() =  zVal;

    if ((givenNrm & curveNrm) >= scalar(0))
    {
        nrmOrientation_ = ALIGNED;
    }
    else
    {
        nrmOrientation_ = OPPOSED;
    }

    Info<< "Initial nrmOrientation after comparison to NURBS u = 0 nrm : "
        << nrmOrientation_
        << endl;
}


void NURBS3DCurve::flipNrmOrientation()
{
    if (nrmOrientation_ == ALIGNED)
    {
        nrmOrientation_ = OPPOSED;
    }
    else
    {
        nrmOrientation_ = ALIGNED;
    }
}


void NURBS3DCurve::setCPs(const List<vector>& CPs)
{
    CPs_ = CPs;
}


void NURBS3DCurve::setWeights(const List<scalar>& weights)
{
    weights_ = weights;
}


void NURBS3DCurve::setName(const word& name)
{
    name_ = name;
}


void NURBS3DCurve::buildCurve()
{
    const label degree(basis_.degree());

    // Iterate over the CPs of this curve.
    forAll(*this, ptI)
    {
        this->operator[](ptI) = vector::zero;

        const scalar u(u_[ptI]);

        // Compute denominator.
        scalar denom(Zero);

        forAll(CPs_, CPJ)
        {
            denom += basis_.basisValue(CPJ, degree, u)*weights_[CPJ];
        }

        // Compute curve point.
        forAll(CPs_, CPI)
        {
            this->operator[](ptI)
                +=   CPs_[CPI]
                   * basis_.basisValue(CPI, degree, u)
                   * weights_[CPI]/denom;
        }
    }
}


void NURBS3DCurve::invert()
{
    Info<< "Inverting NURBS curve " << name_ << endl;

    const label nCPs(CPs_.size());
    List<vector> invertedCPs(nCPs, Zero);
    List<scalar> invertedWeights(nCPs, Zero);

    for (label CPI = 0; CPI<nCPs; CPI++)
    {
        invertedCPs[CPI]     = CPs_[nCPs - 1 - CPI];
        invertedWeights[CPI] = weights_[nCPs - 1 - CPI];
    }

    CPs_ = invertedCPs;
    weights_ = invertedWeights;

    buildCurve();
}


void NURBS3DCurve::makeEquidistant
(
    const label lenAcc,
    const label maxIter,
    const label spacingCorrInterval,
    const scalar tolerance
)
{
/*
    Info<< "Making points equidistant is physical space on curve "
        << name_
        << endl;
*/
    setEquidistantU
    (
        u_,
        lenAcc,
        maxIter,
        spacingCorrInterval,
        tolerance
    );

    buildCurve();
}


// * * * * * * * * * * * * *  Point Calc Functions * * * * * * * * * * * * * //

vector NURBS3DCurve::curvePoint(const scalar u) const
{
    // Compute denominator.
    const label degree(basis_.degree());
    scalar NW(Zero);

    forAll(CPs_, CPI)
    {
        NW += basis_.basisValue(CPI, degree, u)*weights_[CPI];
    }

    // Compute curve point.
    vector point(Zero);

    forAll(CPs_, CPI)
    {
        point +=
            CPs_[CPI]
           *basis_.basisValue(CPI, degree, u)
           *weights_[CPI]/NW;
    }

    return point;
}


scalar NURBS3DCurve::findClosestCurvePoint
(
    const vector& targetPoint,
    const label maxIter,
    const scalar tolerance
)
{
    // Loop through curve points to find the closest one for initialization.
    const vectorField& curve(*this);
    scalar dist(GREAT);
    label closeI(-1);

    forAll(curve, ptI)
    {
        const scalar distLoc(mag(targetPoint - curve[ptI]));

        if (distLoc < dist)
        {
            dist = distLoc;
            closeI = ptI;
        }
    }

    label iter(0);
    scalar u(scalar(closeI)/scalar(this->size()-1));
    vector xu(curvePoint(u));
    scalar res(GREAT);

    do
    {
        vector dxdu(curveDerivativeU(u));
        vector d2xdu2(curveDerivativeUU(u));
        scalar lhs((dxdu&dxdu) + ((xu - targetPoint) & d2xdu2));
        scalar rhs(-((xu - targetPoint) & dxdu));

        // Update parametric coordinate and compute new point and
        // bound param coordinates if needed.
        // Compute residual.
        u += rhs/lhs;

        bound(u);

        xu = curvePoint(u);
        res = mag((xu - targetPoint) & curveDerivativeU(u));

        //Info<< "dxdu" << dxdu << endl;
        //Info<< "d2xdu2" << d2xdu2 << endl;
        //Info<< "Iteration " << iter << ", Residual " << res << endl;
    }
    while ((iter++< maxIter) && (res > tolerance));
/*
    // debug info
    Info<< "Residual after " << iter << " iterations : " << res
        << "\nparametric coordinate " <<  u
        << "\ncurve point closest to target " <<  xu
        << "\ntarget being " <<  targetPoint
        <<  endl;
*/
    // warning if method did not reach threshold
    if (iter > maxIter)
    {
        WarningInFunction
            << "Finding curve point closest to " << targetPoint << " failed."
            << endl;
    }

    return u;
}


scalar NURBS3DCurve::findClosestCurvePoint
(
    const vector& targetPoint,
    const scalar initGuess,
    const label maxIter,
    const scalar tolerance
)
{
    // Loop through curve points to find the closest one for initialization.
    label iter(0);
    scalar u(initGuess);
    vector xu(curvePoint(u));
    scalar res(GREAT);

    do
    {
        vector dxdu(curveDerivativeU(u));
        vector d2xdu2(curveDerivativeUU(u));
        scalar lhs((dxdu&dxdu) + ((xu - targetPoint) & d2xdu2));
        scalar rhs(-((xu - targetPoint) & dxdu));

        // Update parametric coordinate and compute new point and
        // bound param coordinates if needed.
        // Compute residual.
        u += rhs/lhs;

        bound(u);

        xu = curvePoint(u);
        res = mag((xu - targetPoint) & curveDerivativeU(u));

        //Info<< "dxdu" << dxdu << endl;
        //Info<< "d2xdu2" << d2xdu2 << endl;
        //Info<< "bound u " << u << endl;
        //Info<< "u " << u << endl;
        //Info<< "Iteration " << iter << ", Residual " << res << endl;
    }
    while ((iter++< maxIter) && (res > tolerance));
/*
    // debug info
    Info<< "Residual after " << iter << " iterations : " << res
        << "\nparametric coordinate "         <<  u
        << "\ncurve point closest to target " <<  xu
        << "\ntarget being "                   <<  targetPoint
        <<  endl;
*/
    // warning if method did not reach threshold
    if (iter > maxIter)
    {
        WarningInFunction
            << "Finding curve point closest to " << targetPoint
            << " failed."
            << endl;
    }

    return u;
}


const vector NURBS3DCurve::nrm3D(const vector& refTan, const scalar u) const
{
    vector curveNrm(Zero);

    if (nrmOrientation_ == ALIGNED)
    {
        curveNrm = curveDerivativeU(u) ^ refTan;
    }
    else
    {
        curveNrm = refTan ^ curveDerivativeU(u);
    }

    curveNrm.normalise();

    return curveNrm;
}


const vector NURBS3DCurve::nrm2D(const scalar zVal, const scalar u) const
{
    const vector tan(curveDerivativeU(u));
    vector curveNrm(Zero);

    curveNrm.x()  = -nrmOrientation_*tan.y();
    curveNrm.y()  =  nrmOrientation_*tan.x();
    curveNrm.z()  =  zVal;
    curveNrm /= mag(curveNrm);

    return curveNrm;
}


void NURBS3DCurve::insertKnot
(
    const scalarField& oldKnots,
    const scalar uBar,
    const label kInsert
)
{
    // Get the req ref info.
    // Insertion into curve of non-uniform weight is not currently supported.
    const label degree(basis_.degree());
    const label nCPs(basis_.nCPs());
    List<vector> newCPs(nCPs, Zero);
    List<scalar> newWeights(nCPs, scalar(1));

    // Compute the new CPs and their weights.
    for (label CPI = 0; CPI < (kInsert - degree + 1); CPI++)
    {
        newCPs[CPI] = CPs_[CPI];
    }

    for (label CPI = (kInsert - degree + 1); CPI < (kInsert + 1); CPI++)
    {
        const scalar uIOld(oldKnots[CPI]);
        const scalar uIDOld(oldKnots[CPI + degree]);
        const scalar ratio((uBar - uIOld) /(uIDOld - uIOld));

        newCPs[CPI] = (ratio*CPs_[CPI] + (1 - ratio)*CPs_[CPI - 1]);
    }

    for (label CPI= (kInsert + 1); CPI<newCPs.size(); CPI++)
    {
        newCPs[CPI] = CPs_[CPI - 1];
    }

    // Reset the CPs and weights and recompute the curve.
    CPs_ = newCPs;
    weights_ = newWeights;

    buildCurve();
}


scalarList NURBS3DCurve::genEquidistant
(
    const label nPts,
    const label lenAcc,
    const label maxIter,
    const label spacingCorrInterval,
    const scalar tolerance
)
{
/*
    Info<< "Generating points equidistant in physical space on curve "
        << name_
        << endl;
*/
    scalarList U(nPts, Zero);

    setEquidistantU
    (
        U,
        lenAcc,
        maxIter,
        spacingCorrInterval,
        tolerance
    );

    return U;
}


// * * * * * * * * * * * * * *  Location Functions * * * * * * * * * * * * * //

bool NURBS3DCurve::checkRange
(
    const scalar u,
    const label CPI,
    const label degree
) const
{
    return basis_.checkRange(u, CPI, degree);
}


bool NURBS3DCurve::checkRange(const scalar u, const label CPI) const
{
    const label degree(basis_.degree());
    return basis_.checkRange(u, CPI, degree);
}


scalar NURBS3DCurve::length(const label uIStart, const label uIEnd) const
{
    // Compute derivatives wrt u for the given u interval.
    const label lenSize(uIEnd - uIStart + 1);
    vectorField dxdu(lenSize, Zero);
    scalar length(Zero);

    forAll(dxdu, uI)
    {
        dxdu[uI] = curveDerivativeU(u_[uIStart + uI]);
    }

    // Integrate.
    for (label uI = 0; uI < (lenSize - 1); uI++)
    {
        length +=
            0.5
           *(mag(dxdu[uI + 1]) + mag(dxdu[uI]))
           *(u_[uIStart + uI + 1]-u_[uIStart + uI]);
    }

    return length;
}


scalar NURBS3DCurve::length
(
    const scalar uStart,
    const scalar uEnd,
    const label nPts
) const
{
    // Compute derivatives wrt u for the given u interval.
    scalar length(Zero);
    scalarField localU(nPts, Zero);
    vectorField dxdu(nPts, Zero);

    forAll(localU, uI)
    {
        localU[uI] = uStart + scalar(uI)/scalar(nPts - 1)*(uEnd - uStart);
        dxdu[uI]   = curveDerivativeU(localU[uI]);
    }

    // Integrate.
    for (label uI = 0; uI < (nPts - 1); uI++)
    {
        length +=
            0.5
           *(mag(dxdu[uI + 1]) + mag(dxdu[uI]))
           *(localU[uI + 1]-localU[uI]);
    }

    return length;
}


scalar NURBS3DCurve::length() const
{
    return length(scalar(0), (u_.size() - 1));
}


// * * * * * * * * * * * * *  Derivative Functions * * * * * * * * * * * * * //

vector NURBS3DCurve::curveDerivativeU(const scalar u) const
{
    const label degree(basis_.degree());
    vector NWP(Zero);
    vector dNduWP(Zero);
    scalar NW(Zero);
    scalar dNduW(Zero);

    forAll(CPs_, CPI)
    {
        const scalar basisValue(basis_.basisValue(CPI, degree, u));
        const scalar basisDeriv(basis_.basisDerivativeU(CPI, degree, u));

        NWP += basisValue * weights_[CPI] * CPs_[CPI];
        dNduWP += basisDeriv * weights_[CPI] * CPs_[CPI];
        NW += basisValue * weights_[CPI];
        dNduW += basisDeriv * weights_[CPI];
    }

    const vector uDerivative((dNduWP - NWP*dNduW/NW)/NW);

    return uDerivative;
}


vector NURBS3DCurve::curveDerivativeUU(const scalar u) const
{
    const label degree(basis_.degree());
    vector NWP(Zero);
    vector dNduWP(Zero);
    vector d2Ndu2WP(Zero);
    scalar NW(Zero);
    scalar dNduW(Zero);
    scalar d2Ndu2W(Zero);

    forAll(CPs_, CPI)
    {
        const scalar basisValue(basis_.basisValue(CPI, degree, u));
        const scalar basisDeriv(basis_.basisDerivativeU(CPI, degree, u));
        const scalar basis2Deriv(basis_.basisDerivativeUU(CPI, degree, u));

        NWP += basisValue * weights_[CPI] * CPs_[CPI];
        dNduWP += basisDeriv * weights_[CPI] * CPs_[CPI];
        d2Ndu2WP += basis2Deriv * weights_[CPI] * CPs_[CPI];
        NW += basisValue * weights_[CPI];
        dNduW += basisDeriv * weights_[CPI];
        d2Ndu2W += basis2Deriv * weights_[CPI];
    }

    const vector uuDerivative
    (
        (
            d2Ndu2WP
          - scalar(2)*dNduWP*dNduW/NW
          - NWP*d2Ndu2W/NW
          + scalar(2)*NWP*dNduW*dNduW/NW/NW
        ) / NW
    );

    return uuDerivative;
}


scalar NURBS3DCurve::curveDerivativeCP
(
    const scalar u,
    const label CPI
)
{
    // compute denominator.
    const label degree(basis_.degree());
    scalar NW(Zero);

    forAll(CPs_, CPJ)
    {
        NW += basis_.basisValue(CPJ, degree, u) * weights_[CPJ];
    }

    const scalar basisValueI(basis_.basisValue(CPI, degree, u));
    const scalar CPDerivative(basisValueI * weights_[CPI] / NW);

    return CPDerivative;
}


vector NURBS3DCurve::curveDerivativeWeight
(
    const scalar u,
    const label CPI
)
{
    // Compute nominator, denominator.
    const label degree(basis_.degree());
    vector NWP(Zero);
    scalar NW(Zero);

    forAll(CPs_, CPJ)
    {
        const scalar basisValue(basis_.basisValue(CPJ, degree, u));
        NWP += basisValue * weights_[CPJ] * CPs_[CPJ];
        NW += basisValue * weights_[CPJ];
    }

    // Compute derivative.
    const scalar basisValueI(basis_.basisValue(CPI, degree, u));
    const vector WDerivative(basisValueI/NW * (CPs_[CPI] - NWP/NW));

    return WDerivative;
}


scalar NURBS3DCurve::lengthDerivativeU
(
    const scalar uStart,
    const scalar uEnd,
    const label nPts
)
{
    // compute derivatives wrt u for the given u interval
    vectorField dxdu(nPts, Zero);
    vectorField d2xdu2(nPts, Zero);
    scalarField localU(nPts, Zero);
    scalar lDerivative(Zero);

    forAll(localU, uI)
    {
        scalar& uLocal(localU[uI]);
        uLocal = uStart + scalar(uI)/scalar(nPts - 1)*(uEnd - uStart);
        dxdu[uI]   = curveDerivativeU(uLocal);
        d2xdu2[uI] = curveDerivativeUU(uLocal);
    }

    // Integrate.
    for (label uI = 0; uI<(nPts - 1); uI++)
    {
        lDerivative +=
            0.5
          * (
                (dxdu[uI + 1]&d2xdu2[uI + 1])/mag(dxdu[uI + 1])
              + (dxdu[uI]&d2xdu2[uI])/mag(dxdu[uI])
            )
          * (localU[uI + 1]-localU[uI]);
    }

    return lDerivative;
}


// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //

//                             -- Inlined in H--

// * * * * * * * * * * * * * * * Write Functions  * * * * * * * * * * * * * //

void NURBS3DCurve::write()
{
    write(name_);
}


void NURBS3DCurve::write
(
    const word fileName
)
{
    if (Pstream::master())
    {
        OFstream curveFile(fileName);
        OFstream curveFileCPs(fileName + "CPs");
        OFstream curveFileBases(fileName + "Bases");
        label degree(basis_.degree());

        vectorField& field(*this);

        forAll(*this, uI)
        {
            curveFile << field[uI].x() << " "
                      << field[uI].y() << " "
                      << field[uI].z()
                      << endl;
        }

        forAll(CPs_, CPI)
        {
            curveFileCPs << CPs_[CPI].x() << " "
                         << CPs_[CPI].y() << " "
                         << CPs_[CPI].z()
                         << endl;
        }

        forAll(*this, uI)
        {
            const scalar u(u_[uI]);
            scalarField basesValues(CPs_.size());

            curveFileBases << u << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
    }
}


void NURBS3DCurve::write
(
    const fileName dirName,
    const fileName fileName
)
{
    if (Pstream::master())
    {
        OFstream curveFile(dirName/fileName);
        OFstream curveFileCPs(dirName/fileName + "CPs");
        OFstream curveFileBases(dirName/fileName + "Bases");
        label degree(basis_.degree());

        vectorField& field(*this);

        forAll(*this, uI)
        {
            curveFile << field[uI].x() << " "
                      << field[uI].y() << " "
                      << field[uI].z()
                      << endl;
        }

        forAll(CPs_, CPI)
        {
            curveFileCPs << CPs_[CPI].x() << " "
                         << CPs_[CPI].y() << " "
                         << CPs_[CPI].z()
                         << endl;
        }

        forAll(*this, uI)
        {
            const scalar u(u_[uI]);
            scalarField basesValues(CPs_.size());

            curveFileBases << u << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
    }
}


void NURBS3DCurve::writeWParses()
{
    writeWParses(name_);
}


void NURBS3DCurve::writeWParses
(
    const word fileName
)
{
    if (Pstream::master())
    {
        OFstream curveFile(fileName);
        OFstream curveFileCPs(fileName + "CPs");
        OFstream curveFileBases(fileName + "Bases");
        label degree(basis_.degree());

        vectorField& field(*this);

        forAll(*this, uI)
        {
            curveFile << "("
                      << field[uI].x() << " "
                      << field[uI].y() << " "
                      << field[uI].z() << ")"
                      << endl;
        }

        forAll(CPs_, CPI)
        {
            curveFileCPs << "("
                         << CPs_[CPI].x() << " "
                         << CPs_[CPI].y() << " "
                         << CPs_[CPI].z() << ")"
                         << endl;
        }

        forAll(*this, uI)
        {
            const scalar u(u_[uI]);
            scalarField basesValues(CPs_.size());

            curveFileBases << u << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
    }
}


void NURBS3DCurve::writeWParses
(
    const fileName dirName,
    const fileName fileName
)
{
    if (Pstream::master())
    {
        OFstream curveFile(dirName/fileName);
        OFstream curveFileCPs(dirName/fileName + "CPs");
        OFstream curveFileBases(dirName/fileName + "Bases");
        label degree(basis_.degree());

        vectorField& field(*this);

        forAll(*this, uI)
        {
            curveFile << "("
                      << field[uI].x() << " "
                      << field[uI].y() << " "
                      << field[uI].z() << ")"
                      << endl;
        }

        forAll(CPs_, CPI)
        {
            curveFileCPs << "("
                         << CPs_[CPI].x() << " "
                         << CPs_[CPI].y() << " "
                         << CPs_[CPI].z() << ")"
                         << endl;
        }

        forAll(*this, uI)
        {
            const scalar u(u_[uI]);
            scalarField basesValues(CPs_.size());

            curveFileBases << u << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
