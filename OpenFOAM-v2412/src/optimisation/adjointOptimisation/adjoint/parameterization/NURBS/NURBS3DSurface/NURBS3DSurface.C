/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "NURBS3DSurface.H"
#include "vtkSurfaceWriter.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label NURBS3DSurface::sgn(const scalar val) const
{
    return val >= scalar(0) ? 1 : -1;
}


scalar NURBS3DSurface::abs(const scalar val) const
{
    return (sgn(val) == 1)? val: -val;
}


label NURBS3DSurface::mod(const label x, const label interval) const
{
    label ratio(x%interval);
    return ratio < 0 ? ratio+interval : ratio;
}


void NURBS3DSurface::setCPUVLinking()
{
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());

    CPsUCPIs_.setSize(uNCPs*vNCPs, -1);
    CPsVCPIs_.setSize(uNCPs*vNCPs, -1);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);
            CPsUCPIs_[CPI] = uCPI;
            CPsVCPIs_[CPI] = vCPI;
        }
    }
}


void NURBS3DSurface::setUniformUV
(
    scalarList& u,
    scalarList& v,
    const label nUPts,
    const label nVPts
) const
{
    u.setSize(nUPts*nVPts, Zero);
    v.setSize(nUPts*nVPts, Zero);

    for (label uI = 0; uI<nUPts; uI++)
    {
        scalar uValue = scalar(uI)/scalar(nUPts - 1);
        for (label vI = 0; vI<nVPts; vI++)
        {
            const label ptI(uI*nVPts + vI);
            u[ptI] = uValue;
            v[ptI] = scalar(vI)/scalar(nVPts - 1);
        }
    }
}


void NURBS3DSurface::setUniformUV()
{
    setUniformUV(u_, v_, nUPts_, nVPts_);
}


bool NURBS3DSurface::bound
(
    scalar& u,
    scalar& v,
    const scalar minVal,
    const scalar maxVal
) const
{
    bool boundPoint =
        boundDirection(u, minVal, maxVal)
     || boundDirection(v, minVal, maxVal);

    return boundPoint;
}


bool NURBS3DSurface::boundDirection
(
    scalar& u,
    const scalar minVal,
    const scalar maxVal
) const
{
    bool boundPoint(false);

    if (u < scalar(0))
    {
        u = minVal;
        boundPoint = true;
    }
    else if (u > scalar(1))
    {
        u = maxVal;
        boundPoint = true;
    }

    return boundPoint;
}


void NURBS3DSurface::setEquidistantR
(
    scalarList& R,
    const scalar SHeld,
    const label paramR,
    const label lenAcc = 25,
    const label maxIter = 10,
    const label spacingCorrInterval = -1,
    const scalar tolerance = 1.e-5
) const
{
    const label  nPts(R.size());
    scalar SNull(SHeld);
    scalar xLength(Zero);
    const scalar rLength(scalar(1) / scalar(nPts - 1));

    if (paramR == PARAMU)
    {
        xLength = lengthU(SHeld) / (nPts - 1);
    }
    else
    {
        xLength = lengthV(SHeld) / (nPts - 1);
    }

    R[0] = Zero;
    R[nPts-1] = scalar(1);

    for (label ptI=1; ptI<(nPts - 1); ptI++)
    {
        const scalar& RPrev(R[ptI - 1]);
        scalar& RCurr(R[ptI]);
        scalar direc(1);
        scalar xDiff(0);
        scalar delta(0);
        bool overReached(false);

        RCurr = RPrev + rLength;

        // Find the starting U value to ensure target is within 1 rLength.
        while (true)
        {
            bool bounded(false);

            if (paramR == PARAMU)
            {
                bounded = bound(RCurr, SNull);
                delta = lengthU(SHeld, RPrev, RCurr, lenAcc);
            }
            else
            {
                bounded = bound(SNull, RCurr);
                delta = lengthV(SHeld, RPrev, RCurr, lenAcc);
            }

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
                    // rLength addition makes U exceed 1 so it becomes bounded.
                    // However, the desired x location still exceeds how far the
                    // bounded rLength can move (~e-5 error).
                    // Must force U to be u=0.999999.
                    overReached = true;
                    break;
                }
                else if (direc == scalar(1))
                {
                    RCurr += rLength;
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
                RCurr += direc * rLength;

                if (paramR == PARAMU)
                {
                    bound(RCurr, SNull);
                }
                else
                {
                    bound(SNull, RCurr);
                }

                // Can employ an occasional tolerance check from beg of curve.
                if
                (
                    (spacingCorrInterval           != -1)
                 && (mod(ptI, spacingCorrInterval) ==  0)
                )
                {
                    if (paramR == PARAMU)
                    {
                        delta = lengthU(SHeld, Zero, RCurr, ptI*lenAcc);
                    }
                    else
                    {
                        delta = lengthV(SHeld, Zero, RCurr, ptI*lenAcc);
                    }

                    xDiff = (ptI * xLength) - delta;
                }
                else
                {
                    if (paramR == PARAMU)
                    {
                        delta = lengthU(SHeld, RPrev, RCurr, lenAcc);
                    }
                    else
                    {
                        delta = lengthV(SHeld, RPrev, RCurr, lenAcc);
                    }

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

                ++iter;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const label nUPts,
    const label nVPts,
    const NURBSbasis& uBasis,
    const NURBSbasis& vBasis,
    const word name
)
:
    vectorField (nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(CPs.size(), scalar(1)),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_(name),
    uBasis_(uBasis),
    vBasis_(vBasis),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}


NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const List<scalar>& weights,
    const label nUPts,
    const label nVPts,
    const NURBSbasis& uBasis,
    const NURBSbasis& vBasis,
    const word name
)
:
    vectorField(nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(weights),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_ (name),
    uBasis_(uBasis),
    vBasis_(vBasis),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}


NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const label nUPts,
    const label nVPts,
    const label uDegree,
    const label vDegree,
    const label nCPsU,
    const label nCPsV,
    const word name
)
:
    vectorField(nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(CPs.size(), scalar(1)),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_(name),
    uBasis_(nCPsU, uDegree),
    vBasis_(nCPsV, vDegree),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    // Sanity checks
    if (nCPsU*nCPsV != CPs_.size())
    {
        FatalErrorInFunction
            << "nCPsU*nCPsV " << nCPsU*nCPsV
            << " not equal to size of CPs " << CPs_.size()
            << exit(FatalError);
    }
    // Initialize surface
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}


NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const label nUPts,
    const label nVPts,
    const label uDegree,
    const label vDegree,
    const label nCPsU,
    const label nCPsV,
    const scalarField& knotsU,
    const scalarField& knotsV,
    const word name
)
:
    vectorField(nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(CPs.size(), scalar(1)),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_(name),
    uBasis_(nCPsU, uDegree, knotsU),
    vBasis_(nCPsV, vDegree, knotsV),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    // Sanity checks
    if (nCPsU*nCPsV != CPs_.size())
    {
        FatalErrorInFunction
            << "nCPsU*nCPsV " << nCPsU*nCPsV
            << " not equal to size of CPs " << CPs_.size()
            << exit(FatalError);
    }
    // initialize surface
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}


NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const List<scalar>& weights,
    const label nUPts,
    const label nVPts,
    const label uDegree,
    const label vDegree,
    const label nCPsU,
    const label nCPsV,
    const word name
)
:
    vectorField(nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(weights),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_(name),
    uBasis_(nCPsU, uDegree),
    vBasis_(nCPsV, vDegree),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    // Sanity checks
    if (nCPsU*nCPsV != CPs_.size())
    {
        FatalErrorInFunction
            << "nCPsU*nCPsV " << nCPsU*nCPsV
            << " not equal to size of CPs " << CPs_.size()
            << exit(FatalError);
    }

    // Initialize surface
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}


NURBS3DSurface::NURBS3DSurface
(
    const List<vector>& CPs,
    const List<scalar>& weights,
    const label nUPts,
    const label nVPts,
    const label uDegree,
    const label vDegree,
    const label nCPsU,
    const label nCPsV,
    const scalarField& knotsU,
    const scalarField& knotsV,
    const word name
)
:
    vectorField(nUPts*nVPts, Zero),

    CPs_(CPs),
    u_(nUPts*nVPts, Zero),
    v_(nUPts*nVPts, Zero),
    weights_(weights),
    nUPts_(nUPts),
    nVPts_(nVPts),
    name_(name),
    uBasis_(nCPsU, uDegree, knotsU),
    vBasis_(nCPsV, vDegree, knotsV),

    givenInitNrm_(Zero),

    CPsUCPIs_(0),
    CPsVCPIs_(0),

    nrmOrientation_(ALIGNED),

    boundaryCPIDs_(nullptr)
{
    // Sanity checks
    if (nCPsU*nCPsV != CPs_.size())
    {
        FatalErrorInFunction
            << "nCPsU*nCPsV " << nCPsU*nCPsV
            << " not equal to size of CPs " << CPs_.size()
            << exit(FatalError);
    }
    // Initialize surface
    setUniformUV();
    buildSurface();
    setCPUVLinking();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Set Functions * * * * * * * * * * * * * * //

void NURBS3DSurface::setNrmOrientation
(
    const vector& givenNrm,
    const scalar u,
    const scalar v
)
{
    vector surfaceNrm(surfaceDerivativeU(u, v) ^ surfaceDerivativeV(u, v));

    givenInitNrm_ = givenNrm;
    surfaceNrm /= mag(surfaceNrm);

    const scalar relation(givenNrm & surfaceNrm);

    if (relation >= 0)
    {
        nrmOrientation_ = ALIGNED;
    }
    else
    {
        nrmOrientation_ = OPPOSED;
    }

    Info<< "Initial nrmOrientation after comparison to NURBS u="
        << u << ",v=" << v << " nrm: " << nrmOrientation_
        << endl;
}


void NURBS3DSurface::flipNrmOrientation()
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


void NURBS3DSurface::setCPs(const List<vector>& CPs)
{
    CPs_ = CPs;
}


void NURBS3DSurface::setWeights(const scalarList& weights)
{
    weights_ = weights;
}


void NURBS3DSurface::setName(const word& name)
{
    name_ = name;
}


void NURBS3DSurface::buildSurface()
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
/*
    Info<< "\nuDegree: " << uDegree << "\nvDegree: " << vDegree
         << "\nnUPts: "   << nUPts_  << "\nnVPts: "   << nVPts_
         << "\nuNCPs: "   << uNCPs   << "\nvNCPs: "   << vNCPs
         << "\nNURBSSurface:\nCPs: " << CPs
         << endl;
*/
    vectorField& field = *this;
    field = vector::zero;

    for (label uI = 0; uI<nUPts_; uI++)
    {
        for (label vI = 0; vI<nVPts_; vI++)
        {
            const label ptI(uI*nVPts_ + vI);
            const scalar& u(u_[ptI]);
            const scalar& v(v_[ptI]);
            scalar NMW(Zero);

            // Compute denominator.
            for (label vCPI = 0; vCPI<vNCPs; vCPI++)
            {
                for (label uCPI = 0; uCPI<uNCPs; uCPI++)
                {
                    const label CPI(vCPI*uNCPs + uCPI);

                    NMW +=
                        uBasis_.basisValue(uCPI, uDegree, u)
                      * vBasis_.basisValue(vCPI, vDegree, v)
                      * weights_[CPI];
                }
            }

            // Compute the points.
            for (label vCPI = 0; vCPI<vNCPs; vCPI++)
            {
                for (label uCPI = 0; uCPI<uNCPs; uCPI++)
                {
                    const label CPI(vCPI*uNCPs + uCPI);

                    this->operator[](ptI) +=
                             CPs_[CPI]
                           * uBasis_.basisValue(uCPI, uDegree, u)
                           * vBasis_.basisValue(vCPI, vDegree, v)
                           * weights_[CPI]/NMW;
                }
            }
        }
    }
}



void NURBS3DSurface::invertU()
{
    Info<< "Inverting NURBS surface " << name_ << " in u." << endl;

    const label  uNCPs(uBasis_.nCPs());
    const label  vNCPs(vBasis_.nCPs());
    List<vector> invertedCPs(CPs_.size(), Zero);
    List<scalar> invertedWeights(CPs_.size(), Zero);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);
            const label invUCPI(uNCPs-1-uCPI);
            const label uInvCPI(vCPI*uNCPs + invUCPI);

            invertedCPs[CPI] = CPs_[uInvCPI];
            invertedWeights[CPI] = weights_[uInvCPI];
        }
    }

    CPs_ = invertedCPs;
    weights_ = invertedWeights;

    buildSurface();
}


void NURBS3DSurface::invertV()
{
    Info<< "Inverting NURBS surface " << name_ << " in v." << endl;

    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    List<vector> invertedCPs(CPs_.size(), Zero);
    List<scalar> invertedWeights(CPs_.size(), Zero);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);
            const label invVCPI(vNCPs-1-vCPI);
            const label vInvCPI(invVCPI*uNCPs + uCPI);

            invertedCPs[CPI] = CPs_[vInvCPI];
            invertedWeights[CPI] = weights_[vInvCPI];
        }
    }

    CPs_ = invertedCPs;
    weights_ = invertedWeights;

    buildSurface();
}


void NURBS3DSurface::invertUV()
{
    Info<< "Inverting NURBS surface " << name_ << " in u and v." << endl;

    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    List<vector> invertedCPs(CPs_.size(), Zero);
    List<scalar> invertedWeights(CPs_.size(), Zero);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);
            const label invUCPI(uNCPs - 1 - uCPI);
            const label invVCPI(vNCPs - 1 - vCPI);
            const label uvInvCPI(invVCPI*uNCPs + invUCPI);

            invertedCPs[CPI] = CPs_[uvInvCPI];
            invertedWeights[CPI] = weights_[uvInvCPI];
        }
    }

    CPs_ = invertedCPs;
    weights_ = invertedWeights;

    buildSurface();
}


void NURBS3DSurface::makeEquidistant
(
    const label lenAcc,
    const label maxIter,
    const label spacingCorrInterval,
    const scalar tolerance
)
{
/*
    Info<< "Making points equidistant is physical space on surface "
         << name_
         << endl;
*/
    // Equidistant spacing in u along v isoLines.
    for (label vI = 0; vI<nVPts_; vI++)
    {
        scalarList UI(nUPts_, Zero);
        const scalar VHeld(v_[vI]);
        labelList uAddressing(nUPts_, -1);

        // Set the point uAddressing to re-assign correct u_ values later.
        forAll(uAddressing, uI)
        {
            const label ptI(uI*nVPts_ + vI);
            uAddressing[uI] = ptI;
        }
        // Set equidistant u values.
        setEquidistantR
        (
            UI,
            VHeld,
            PARAMU,
            lenAcc,
            maxIter,
            spacingCorrInterval,
            tolerance
        );

        // Re-assign new equidistant u values.
        forAll(UI, uI)
        {
            const label& uAddress(uAddressing[uI]);
            u_[uAddress] = UI[uI];
        }
    }

    // Equidistant spacing in v along u isoLines.
    for (label uI = 0; uI<nUPts_; uI++)
    {
        scalarList VI(nVPts_, Zero);
        const scalar UHeld(u_[uI*nVPts_]);
        labelList vAddressing(nUPts_, -1);

        // Set the point vAddressing to re-assign correct u_ values later.
        forAll(vAddressing, vI)
        {
            const label ptI(uI*nVPts_ + vI);
            vAddressing[vI] = ptI;
        }

        // Set equidistant u values.
        setEquidistantR
        (
            VI,
            UHeld,
            PARAMV,
            lenAcc,
            maxIter,
            spacingCorrInterval,
            tolerance
        );

        // Re-assign new equidistant u values.
        forAll(VI, vI)
        {
            const label& vAddress(vAddressing[vI]);
            v_[vAddress] = VI[vI];
        }
    }

    buildSurface();
}


// * * * * * * * * * * * * *  Point Calc Functions * * * * * * * * * * * * * //

vector NURBS3DSurface::surfacePoint
(
    const scalar& u,
    const scalar& v
)
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    scalar NMW(Zero);

    // Compute denominator.
    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);

            NMW +=
                uBasis_.basisValue(uCPI, uDegree, u)
              * vBasis_.basisValue(vCPI, vDegree, v)
              * weights_[CPI];
        }
    }

    // Compute the points.
    vector point(Zero);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);

            point +=
                CPs_[CPI]
              * uBasis_.basisValue(uCPI, uDegree, u)
              * vBasis_.basisValue(vCPI, vDegree, v)
              * weights_[CPI]/NMW;
        }
    }

    return point;
}


scalarList NURBS3DSurface::findClosestSurfacePoint
(
    const vector& targetPoint,
    const label maxIter,
    const scalar tolerance
)
{
    // Loop through surface points to find the closest one for initialization.
    const vectorField& surface(*this);
    scalar dist(GREAT);
    label closePtI(-1);

    forAll(surface, ptI)
    {
        const scalar distLoc(mag(targetPoint-surface[ptI]));

        if (distLoc < dist)
        {
            dist = distLoc;
            closePtI = ptI;
        }
    }

    label iter(0);
    scalar u(u_[closePtI]);
    scalar v(v_[closePtI]);
    vector xuv(surfacePoint(u, v));
    scalar res(GREAT);
    scalar resOld(GREAT);
    scalar resDeriv(GREAT);
    label nBoundsU(0);
    label nBoundsV(0);

    do
    {
        /*
        const vector dxdu(surfaceDerivativeU(u, v));
        const vector dxdv(surfaceDerivativeV(u, v));
        const vector d2xdu2(surfaceDerivativeUU(u, v));
        const vector d2xdv2(surfaceDerivativeVV(u, v));
        const scalar uLHS((dxdu&dxdu) + ((xuv-targetPoint) & d2xdu2));
        const scalar uRHS(-((xuv-targetPoint) & dxdu));
        const scalar vLHS((dxdv&dxdv) + ((xuv-targetPoint) & d2xdv2));
        const scalar vRHS(-((xuv-targetPoint) & dxdv));

        // Update parametric coordinate and compute new point and
        // bound param coordinates if needed.
        // Compute residual.
        u += uRHS/(uLHS+SMALL);
        v += vRHS/(vLHS+SMALL);

        bound(u, v);

        xuv = surfacePoint(u, v);
        res =   mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
              + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
        */

        const vector dxdu(surfaceDerivativeU(u, v));
        const vector dxdv(surfaceDerivativeV(u, v));
        const vector d2xdu2(surfaceDerivativeUU(u, v));
        const vector d2xdv2(surfaceDerivativeVV(u, v));
        const vector d2xduv(surfaceDerivativeUV(u, v));
        const scalar a((dxdu&dxdu) + ((xuv-targetPoint) & d2xdu2));
        const scalar b((dxdu&dxdv) + ((xuv-targetPoint) & d2xduv));
        const scalar c=b;
        const scalar d((dxdv&dxdv) + ((xuv-targetPoint) & d2xdv2));
        const scalar invDenom = 1./(a*d-b*c);

        const scalar uRHS(-((xuv-targetPoint) & dxdu));
        const scalar vRHS(-((xuv-targetPoint) & dxdv));

        // Update parametric coordinate and compute new point and
        // bound param coordinates if needed.
        // Compute residual.
        u += ( d*uRHS-b*vRHS)*invDenom;
        v += (-c*uRHS+a*vRHS)*invDenom;

        if (boundDirection(u))
        {
            nBoundsU++;
        }
        if (boundDirection(v))
        {
            nBoundsV++;
        }

        xuv = surfacePoint(u, v);
        // If manual assignment in u is required, deal only with the v eqn
        if (nBoundsU >= 5)
        {
            res = mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
            resDeriv = mag(res-resOld)/resOld;
            resOld = res;
//          Info<< targetPoint << " " << res << endl;
        }
        // If manual assignment in v is required, deal only with the u eqn
        else if (nBoundsV >= 5)
        {
            res = mag((xuv-targetPoint) & surfaceDerivativeU(u, v));
            resDeriv = mag(res-resOld)/resOld;
            resOld = res;
//          Info<< targetPoint << " " << res << endl;
        }
        else if (nBoundsU <= 5 && nBoundsV <= 5)
        {
            res = mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
                + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
            resDeriv = mag(res-resOld)/resOld;
            resOld = res;
        }
        else
        {
            WarningInFunction
                << "More than 5 bounds in both the u and v directions!"
                << "Something seems weird" << nBoundsU << " " << nBoundsV
                << endl;
            res  = mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
                 + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
            resDeriv = mag(res-resOld)/resOld;
            resOld = res;
        }
    }
    while ((iter++ < maxIter) && (res > tolerance));

    scalarList closestParameters(2, u);
    closestParameters[1] = v;
    // warning if method did not reach threshold
    if (iter > maxIter)
    {
        WarningInFunction
            << "Finding surface point closest to " << targetPoint
            << " for surface " << name_ << " failed \n"
            << " Number of bounding operations in u,v "
            << nBoundsU << " " << nBoundsV << endl
            << " Residual value and derivative " << res << " " << resDeriv
            << endl << endl;
        closestParameters = -1;
    }

    return closestParameters;
}


tmp<vector2DField> NURBS3DSurface::findClosestSurfacePoint
(
    const vectorField& targetPoints,
    const label maxIter,
    const scalar tolerance
)
{
    auto tparamCoors = tmp<vector2DField>::New(targetPoints.size(), Zero);
    auto& paramCoors = tparamCoors.ref();

    const vectorField& surface(*this);
    label nBoundedPoints(0);
    scalar maxResidual(0);
    scalar maxResidualDeriv(0);
    forAll(paramCoors, pI)
    {
        const vector& targetPoint(targetPoints[pI]);

        // Loop through surface points to find the closest one for
        // initialization.  The initialization could possibly be done with a
        // geodesical Laplace. Potentially faster?
        scalar dist(GREAT);
        label closePtI(-1);

        forAll(surface, ptI)
        {
            const scalar distLoc(mag(targetPoint - surface[ptI]));

            if (distLoc < dist)
            {
                dist     = distLoc;
                closePtI = ptI;
            }
        }

        label iter(0);
        scalar u(u_[closePtI]);
        scalar v(v_[closePtI]);
        vector xuv(surfacePoint(u, v));
        scalar res(GREAT);
        scalar resOld(GREAT);
        scalar resDeriv(GREAT);
        label nBoundsU(0);
        label nBoundsV(0);

        do
        {
            const vector dxdu(surfaceDerivativeU(u, v));
            const vector dxdv(surfaceDerivativeV(u, v));
            const vector d2xdu2(surfaceDerivativeUU(u, v));
            const vector d2xdv2(surfaceDerivativeVV(u, v));
            const vector d2xduv(surfaceDerivativeUV(u, v));
            const scalar a((dxdu&dxdu) + ((xuv-targetPoint) & d2xdu2));
            const scalar b((dxdu&dxdv) + ((xuv-targetPoint) & d2xduv));
            const scalar c=b;
            const scalar d((dxdv&dxdv) + ((xuv-targetPoint) & d2xdv2));
            const scalar invDenom = 1./(a*d-b*c);

            const scalar uRHS(-((xuv-targetPoint) & dxdu));
            const scalar vRHS(-((xuv-targetPoint) & dxdv));

            // Update parametric coordinate and compute new point and
            // bound param coordinates if needed.
            // Compute residual.
            u += ( d*uRHS-b*vRHS)*invDenom;
            v += (-c*uRHS+a*vRHS)*invDenom;

            if (boundDirection(u))
            {
                nBoundsU++;
            }
            if (boundDirection(v))
            {
                nBoundsV++;
            }

            xuv = surfacePoint(u, v);
            // If manual assignment in u is required, deal only with the v eqn
            if (nBoundsU >= 5)
            {
                res = mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
                resDeriv = mag(res-resOld)/resOld;
                resOld = res;
            }
            // If manual assignment in b is required, deal only with the u eqn
            else if (nBoundsV >= 5)
            {
                res = mag((xuv-targetPoint) & surfaceDerivativeU(u, v));
                resDeriv = mag(res-resOld)/resOld;
                resOld = res;
            }
            else if (nBoundsU<=5 && nBoundsV<=5)
            {
                res = mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
                         + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
                resDeriv = mag(res-resOld)/resOld;
                resOld = res;
            }
            else
            {
                WarningInFunction
                    << "More than 5 bounding operations in both the u and v directions!"
                    << "u direction " << nBoundsU << endl
                    << "v direction " << nBoundsV << endl
                    << endl;
                res =
                    mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
                  + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
                resDeriv = mag(res-resOld)/resOld;
                resOld = res;
            }
        }
        while ((iter++ < maxIter) && (res > tolerance));

        // warning if method did not reach threshold
        if (iter > maxIter)
        {
            nBoundedPoints++;
            maxResidual = max(maxResidual, res);
            maxResidualDeriv = max(maxResidualDeriv, resDeriv);
        }

        paramCoors[pI].x() = u;
        paramCoors[pI].y() = v;
    }

    Info<< "Points that couldn't reach the residual limit::             "
        << returnReduce(nBoundedPoints, sumOp<label>()) << nl
        << "Max residual after reaching iterations limit::              "
        << returnReduce(maxResidual, maxOp<scalar>()) << nl
        << "Max residual derivative after reaching iterations limit::   "
        << returnReduce(maxResidualDeriv, maxOp<scalar>()) << nl
        << endl;

    return tparamCoors;
}


scalarList NURBS3DSurface::findClosestSurfacePoint
(
    const vector& targetPoint,
    const scalar& uInitGuess,
    const scalar& vInitGuess,
    const label maxIter,
    const scalar tolerance
)
{
    // Loop through surface points to find the closest one for initialization.
    label iter(0);
    scalar u(uInitGuess);
    scalar v(vInitGuess);
    vector xuv(surfacePoint(u, v));
    scalar res(GREAT);

    do
    {
        const vector dxdu(surfaceDerivativeU(u, v));
        const vector dxdv(surfaceDerivativeV(u, v));
        const vector d2xdu2(surfaceDerivativeUU(u, v));
        const vector d2xdv2(surfaceDerivativeVV(u, v));
        const scalar uLHS((dxdu&dxdu) + ((xuv-targetPoint) & d2xdu2));
        const scalar uRHS(-((xuv-targetPoint) & dxdu));
        const scalar vLHS((dxdv&dxdv) + ((xuv-targetPoint) & d2xdv2));
        const scalar vRHS(-((xuv-targetPoint) & dxdv));

        // Update parametric coordinate and compute new point and
        // bound param coordinates if needed.
        // Compute residual.
        u += uRHS/(uLHS+SMALL);
        v += vRHS/(vLHS+SMALL);

        bound(u, v);

        xuv = surfacePoint(u, v);
        res =
            mag((xuv-targetPoint) & surfaceDerivativeU(u, v))
          + mag((xuv-targetPoint) & surfaceDerivativeV(u, v));
    }
    while ((iter++ < maxIter) && (res > tolerance));

    // warning if method did not reach threshold
    if (iter > maxIter)
    {
        WarningInFunction
            << "Finding surface point closest to " << targetPoint << " failed."
            << endl;
    }

    scalarList closestParameters(2, u);
    closestParameters[1] = v;

    return closestParameters;
}


const vector NURBS3DSurface::nrm(scalar u, scalar v)
{
    vector surfaceNrm(Zero);

    if (nrmOrientation_ == ALIGNED)
    {
        surfaceNrm = surfaceDerivativeU(u, v) ^ surfaceDerivativeV(u, v);
    }
    else
    {
        surfaceNrm = surfaceDerivativeV(u, v) ^ surfaceDerivativeU(u, v);
    }

    surfaceNrm /= mag(surfaceNrm);

    return surfaceNrm;
}


List<scalarList> NURBS3DSurface::genEquidistant
(
    const label  nUPts,
    const label  nVPts,
    const label  lenAcc,
    const label  maxIter,
    const label  spacingCorrInterval,
    const scalar tolerance
)
{
/*
    Info<< "Generating points equidistant in physical space on surface "
        << name_
        << endl;
*/
    // Preset U and V with uniform values.
    List<scalarList> UV(NPARAMS, scalarList(0));

    scalarList& U(UV[PARAMU]);
    scalarList& V(UV[PARAMV]);

    setUniformUV(U, V, nUPts, nVPts);

    // Equidistant spacing in u along v isoLines.
    for (label vI = 0; vI<nVPts; vI++)
    {
        scalarList UI(nUPts, Zero);
        const scalar VHeld(V[vI]);
        labelList uAddressing(nUPts, -1);

        // Set the point uAddressing to re-assign correct U values later.
        forAll(uAddressing, uI)
        {
            const label ptI(uI*nVPts + vI);
            uAddressing[uI] = ptI;
        }
        // Set equidistant u values.
        setEquidistantR
        (
            UI,
            VHeld,
            PARAMU,
            lenAcc,
            maxIter,
            spacingCorrInterval,
            tolerance
        );

        // Re-assign new equidistant u values.
        forAll(UI, uI)
        {
            const label& uAddress(uAddressing[uI]);
            U[uAddress] = UI[uI];
        }
    }

    // Equidistant spacing in v along u isoLines.
    for (label uI = 0; uI<nUPts; uI++)
    {
        scalarList VI(nVPts, Zero);
        const scalar UHeld(U[uI*nVPts]);
        labelList vAddressing(nUPts, -1);

        // Set the point vAddressing to re-assign correct V values later.
        forAll(vAddressing, vI)
        {
            const label ptI(uI*nVPts + vI);
            vAddressing[vI] = ptI;
        }

        // Set equidistant u values.
        setEquidistantR
        (
            VI,
            UHeld,
            PARAMV,
            lenAcc,
            maxIter,
            spacingCorrInterval,
            tolerance
        );

        // Re-assign new equidistant u values.
        forAll(VI, vI)
        {
            const label& vAddress(vAddressing[vI]);
            V[vAddress] = VI[vI];
        }
    }

    return UV;
}


// * * * * * * * * * * * * * *  Location Functions * * * * * * * * * * * * * //

bool NURBS3DSurface::checkRangeU
(
    const scalar u,
    const label CPI,
    const label uDegree
) const
{
    const label uCPI(CPsUCPIs_[CPI]);

    return uBasis_.checkRange(u, uCPI, uDegree);
}


bool NURBS3DSurface::checkRangeU
(
    const scalar u,
    const label CPI
) const
{
    const label uDegree(uBasis_.degree());

    return checkRangeU(u, CPI, uDegree);
}


bool NURBS3DSurface::checkRangeV
(
    const scalar v,
    const label CPI,
    const label vDegree
) const
{
    const label vCPI(CPsVCPIs_[CPI]);

    return vBasis_.checkRange(v, vCPI, vDegree);
}


bool NURBS3DSurface::checkRangeV
(
    const scalar v,
    const label CPI
) const
{
    const label vDegree(vBasis_.degree());

    return checkRangeV(v, CPI, vDegree);
}


bool NURBS3DSurface::checkRangeUV
(
    const scalar v,
    const scalar u,
    const label CPI,
    const label uDegree,
    const label vDegree
) const
{
    if (checkRangeU(u, CPI, uDegree) && checkRangeV(v, CPI, vDegree))
    {
        return true;
    }

    return false;
}


bool NURBS3DSurface::checkRangeUV
(
    const scalar v,
    const scalar u,
    const label CPI
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());

    return checkRangeUV(u, v, CPI, uDegree, vDegree);
}


scalar NURBS3DSurface::lengthU
(
    const label vIConst,
    const label uIStart,
    const label uIEnd
) const
{
    // Compute derivatives wrt u for the given u interval.
    const label uLenSize(uIEnd - uIStart + 1);
    vectorField dxdu(uLenSize, Zero);
    scalar uLength(Zero);

    forAll(dxdu, uI)
    {
        const label ptI((uIStart+uI)*nVPts_ + vIConst);
        const label& u(u_[ptI]);
        const label& v(v_[ptI]);

        dxdu[uI] = surfaceDerivativeU(u, v);
    }

    // Integrate.
    for (label uI = 0; uI<(uLenSize - 1); uI++)
    {
        const label ptI((uIStart+uI)*nVPts_ + vIConst);

        uLength +=
            0.5*(mag(dxdu[uI + 1]) + mag(dxdu[uI]))*(u_[ptI + 1]-u_[ptI]);
    }

    return uLength;
}


scalar NURBS3DSurface::lengthU
(
    const scalar vConst,
    const scalar uStart,
    const scalar uEnd,
    const label nPts
) const
{
    // Compute derivatives wrt u for the given u interval.
    vectorField dxdu(nPts, Zero);
    scalarField localU(nPts, Zero);
    scalar uLength(Zero);

    forAll(localU, uI)
    {
        scalar& uLocal(localU[uI]);
        uLocal = uStart + scalar(uI)/scalar(nPts-1)*(uEnd-uStart);
        dxdu[uI] = surfaceDerivativeU(uLocal, vConst);
    }

    // Integrate.
    for (label uI = 0; uI<(nPts - 1); uI++)
    {
        uLength +=
            0.5*(mag(dxdu[uI + 1]) + mag(dxdu[uI]))*(localU[uI + 1]-localU[uI]);
    }

    return uLength;
}


scalar NURBS3DSurface::lengthU(const label vIConst) const
{
    return lengthU(vIConst, 0, (nUPts_ - 1));
}


scalar NURBS3DSurface::lengthU(const scalar vConst) const
{
    return lengthU(vConst, scalar(0), scalar(1), 100);
}


scalar NURBS3DSurface::lengthV
(
    const label uIConst,
    const label vIStart,
    const label vIEnd
) const
{
    // Compute derivatives wrt v for the given v interval.
    const label vLenSize(vIEnd - vIStart + 1);
    vectorField dxdv(vLenSize, Zero);
    scalar vLength(Zero);

    forAll(dxdv, vI)
    {
        const label ptI((uIConst)*nVPts_ + (vIStart+vI));
        const label& u(u_[ptI]);
        const label& v(v_[ptI]);

        dxdv[vI] = surfaceDerivativeV(u, v);
    }

    // Integrate.
    for (label vI = 0; vI<(vLenSize - 1); vI++)
    {
        const label ptI((uIConst)*nVPts_ + (vIStart + vI));

        vLength +=
            0.5*(mag(dxdv[vI + 1]) + mag(dxdv[vI]))*(v_[ptI + 1] - v_[ptI]);
    }

    return vLength;
}


scalar NURBS3DSurface::lengthV
(
    const scalar uConst,
    const scalar vStart,
    const scalar vEnd,
    const label nPts
) const
{
    // Compute derivatives wrt v for the given v interval.
    vectorField dxdv(nPts, Zero);
    scalarField localV(nPts, Zero);
    scalar vLength(Zero);

    forAll(localV, vI)
    {
        scalar& vLocal(localV[vI]);
        vLocal = vStart + scalar(vI)/scalar(nPts - 1)*(vEnd - vStart);
        dxdv[vI] = surfaceDerivativeV(uConst, vLocal);
    }

    // Integrate.
    for (label vI = 0; vI<(nPts - 1); vI++)
    {
        vLength +=
            0.5*(mag(dxdv[vI + 1]) + mag(dxdv[vI]))*(localV[vI + 1]-localV[vI]);
    }

    return vLength;
}


scalar NURBS3DSurface::lengthV(const label uIConst) const
{
    return lengthV(uIConst, 0, (nVPts_ - 1));
}


scalar NURBS3DSurface::lengthV(const scalar uConst) const
{
    return lengthV(uConst, scalar(0), scalar(1), 100);
}


// * * * * * * * * * * * * *  Derivative Functions * * * * * * * * * * * * * //

vector NURBS3DSurface::surfaceDerivativeU
(
    const scalar uIn,
    const scalar vIn
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    vector NMWP(Zero);
    vector dNduMWP(Zero);
    scalar NMW(Zero);
    scalar dNduMW(Zero);

    scalar u = uIn;
    scalar v = vIn;
    bound(u, v);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label CPI(vCPI*uNCPs + uCPI);
            const scalar uBasisValue(uBasis_.basisValue(uCPI, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPI, vDegree, v));
            const scalar uBasisDeriv
                (uBasis_.basisDerivativeU(uCPI, uDegree, u));

            NMWP    += uBasisValue * vBasisValue * weights_[CPI] * CPs_[CPI];
            dNduMWP += uBasisDeriv * vBasisValue * weights_[CPI] * CPs_[CPI];
            NMW     += uBasisValue * vBasisValue * weights_[CPI];
            dNduMW  += uBasisDeriv * vBasisValue * weights_[CPI];
        }
    }

    const vector uDerivative((dNduMWP - dNduMW*NMWP/(NMW+SMALL))/(NMW+SMALL));

    return uDerivative;
}


vector NURBS3DSurface::surfaceDerivativeV
(
    const scalar uIn,
    const scalar vIn
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    vector NMWP(Zero);
    vector dMdvNWP(Zero);
    scalar NMW(Zero);
    scalar dMdvNW(Zero);

    scalar u = uIn;
    scalar v = vIn;
    bound(u, v);

    for (label vCPI = 0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI = 0; uCPI<uNCPs; uCPI++)
        {
            const label  CPI(vCPI*uNCPs + uCPI);
            const scalar uBasisValue(uBasis_.basisValue(uCPI, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPI, vDegree, v));
            const scalar vBasisDeriv
                (vBasis_.basisDerivativeU(vCPI, vDegree, v));

            NMWP    += uBasisValue * vBasisValue * weights_[CPI] * CPs_[CPI];
            dMdvNWP += vBasisDeriv * uBasisValue * weights_[CPI] * CPs_[CPI];
            NMW     += uBasisValue * vBasisValue * weights_[CPI];
            dMdvNW  += vBasisDeriv * uBasisValue * weights_[CPI];
        }
    }

    const vector vDerivative((dMdvNWP - dMdvNW*NMWP/(NMW+SMALL))/(NMW+SMALL));

    return vDerivative;
}


vector NURBS3DSurface::surfaceDerivativeUV
(
    const scalar uIn,
    const scalar vIn
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    vector NMWP(Zero);
    vector dNduMWP(Zero);
    vector dMdvNWP(Zero);
    vector dNMduvWP(Zero);
    scalar NMW(Zero);
    scalar dNduMW(Zero);
    scalar dMdvNW(Zero);
    scalar dNMduvW(Zero);

    scalar u = uIn;
    scalar v = vIn;
    bound(u, v);

    for (label vCPI=0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI=0; uCPI<uNCPs; uCPI++)
        {
            const label  CPI(vCPI*uNCPs + uCPI);
            const scalar uBasisValue(uBasis_.basisValue(uCPI, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPI, vDegree, v));
            const scalar uBasisDeriv
                (uBasis_.basisDerivativeU(uCPI, uDegree, u));
            const scalar vBasisDeriv
                (vBasis_.basisDerivativeU(vCPI, vDegree, v));
            //Info<< "vCPI=" << vCPI << ",uCPI=" << uCPI
            //    << "  N=" << uBasisValue << ",N'=" << uBasisDeriv
            //    << "  M=" << vBasisValue << ",M'=" << vBasisDeriv
            //    << endl;
            NMWP     += vBasisValue * uBasisValue * weights_[CPI] * CPs_[CPI];
            dNduMWP  += uBasisDeriv * vBasisValue * weights_[CPI] * CPs_[CPI];
            dMdvNWP  += vBasisDeriv * uBasisValue * weights_[CPI] * CPs_[CPI];
            dNMduvWP += uBasisDeriv * vBasisDeriv * weights_[CPI] * CPs_[CPI];
            NMW      += vBasisValue * uBasisValue * weights_[CPI];
            dNduMW   += uBasisDeriv * vBasisValue * weights_[CPI];
            dMdvNW   += vBasisDeriv * uBasisValue * weights_[CPI];
            dNMduvW  += uBasisDeriv * vBasisDeriv * weights_[CPI];
        }
    }

    const vector uvDerivative
    (
        (
            dNMduvWP
          - (dNMduvW*NMWP + dMdvNW*dNduMWP + dNduMW*dMdvNWP)/(NMW+SMALL)
          + scalar(2)*dNduMW*dMdvNW*NMWP/(NMW+SMALL)/(NMW+SMALL)
        ) / (NMW+SMALL)
    );

    return uvDerivative;
}


vector NURBS3DSurface::surfaceDerivativeUU
(
    const scalar uIn,
    const scalar vIn
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    vector NMWP(Zero);
    vector dNduMWP(Zero);
    vector d2Ndu2MWP(Zero);
    scalar NMW(Zero);
    scalar dNduMW(Zero);
    scalar d2Ndu2MW(Zero);

    scalar u = uIn;
    scalar v = vIn;
    bound(u, v);

    for (label vCPI=0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI=0; uCPI<uNCPs; uCPI++)
        {
            const label  CPI(vCPI*uNCPs + uCPI);
            const scalar uBasisValue(uBasis_.basisValue(uCPI, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPI, vDegree, v));
            const scalar uBasisDeriv
                (uBasis_.basisDerivativeU(uCPI, uDegree, u));
            const scalar uBasis2Deriv
                (uBasis_.basisDerivativeUU(uCPI, uDegree, u));

            NMWP      += uBasisValue  * vBasisValue * weights_[CPI] * CPs_[CPI];
            dNduMWP   += uBasisDeriv  * vBasisValue * weights_[CPI] * CPs_[CPI];
            d2Ndu2MWP += uBasis2Deriv * vBasisValue * weights_[CPI] * CPs_[CPI];
            NMW       += uBasisValue  * vBasisValue * weights_[CPI];
            dNduMW    += uBasisDeriv  * vBasisValue * weights_[CPI];
            d2Ndu2MW  += uBasis2Deriv * vBasisValue * weights_[CPI];
        }
    }

    const vector uuDerivative
    (
        (
            d2Ndu2MWP
          - scalar(2)*dNduMW*dNduMWP/(NMW+SMALL)
          - d2Ndu2MW*NMWP/(NMW+SMALL)
          + scalar(2)*dNduMW*dNduMW*NMWP/(NMW+SMALL)/(NMW+SMALL)
        ) / (NMW+SMALL)
    );

    return uuDerivative;
}


vector NURBS3DSurface::surfaceDerivativeVV
(
    const scalar uIn,
    const scalar vIn
) const
{
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    vector NMWP(Zero);
    vector dMdvNWP(Zero);
    vector d2Mdv2NWP(Zero);
    scalar NMW(Zero);
    scalar dMdvNW(Zero);
    scalar d2Mdv2NW(Zero);

    scalar u = uIn;
    scalar v = vIn;
    bound(u, v);

    for (label vCPI=0; vCPI<vNCPs; vCPI++)
    {
        for (label uCPI=0; uCPI<uNCPs; uCPI++)
        {
            const label  CPI(vCPI*uNCPs + uCPI);
            const scalar uBasisValue(uBasis_.basisValue(uCPI, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPI, vDegree, v));
            const scalar vBasisDeriv
                (vBasis_.basisDerivativeU(vCPI, vDegree, v));
            const scalar vBasis2Deriv
                (vBasis_.basisDerivativeUU(vCPI, vDegree, v));

            NMWP      += vBasisValue  * uBasisValue * weights_[CPI] * CPs_[CPI];
            dMdvNWP   += vBasisDeriv  * uBasisValue * weights_[CPI] * CPs_[CPI];
            d2Mdv2NWP += vBasis2Deriv * uBasisValue * weights_[CPI] * CPs_[CPI];
            NMW       += vBasisValue  * uBasisValue * weights_[CPI];
            dMdvNW    += vBasisDeriv  * uBasisValue * weights_[CPI];
            d2Mdv2NW  += vBasis2Deriv * uBasisValue * weights_[CPI];
        }
    }

    const vector vvDerivative
    (
        (
            d2Mdv2NWP
          - scalar(2)*dMdvNW*dMdvNWP/(NMW+SMALL)
          - d2Mdv2NW*NMWP/(NMW+SMALL)
          + scalar(2)*dMdvNW*dMdvNW*NMWP/(NMW+SMALL)/(NMW+SMALL)
        ) / (NMW+SMALL)
    );

    return vvDerivative;
}


scalar NURBS3DSurface::surfaceDerivativeCP
(
    const scalar u,
    const scalar v,
    const label CPI
) const
{
    //Info<< "u,v,cpI " << u << " " << v << " " << CPI << endl;
    // compute denominator.
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    const label uCPI(CPsUCPIs_[CPI]);
    const label vCPI(CPsVCPIs_[CPI]);
    scalar NMW(Zero);

    for (label vCPJ=0; vCPJ<vNCPs; vCPJ++)
    {
        for (label uCPJ=0; uCPJ<uNCPs; uCPJ++)
        {
            const label CPJ(vCPJ*uNCPs + uCPJ);
            const scalar uBasisValue(uBasis_.basisValue(uCPJ, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPJ, vDegree, v));

            NMW += vBasisValue * uBasisValue * weights_[CPJ];
        }
    }
    //Info<< "denom " << NMW << endl;

    const scalar CPDerivative
    (
          uBasis_.basisValue(uCPI, uDegree, u)
        * vBasis_.basisValue(vCPI, vDegree, v)
        * weights_[CPI]
        / (NMW+SMALL)
    );

    return CPDerivative;
}


vector NURBS3DSurface::surfaceDerivativeW
(
    const scalar u,
    const scalar v,
    const label CPI
) const
{
    // Compute nominator, denominator.
    const label uDegree(uBasis_.degree());
    const label vDegree(vBasis_.degree());
    const label uNCPs(uBasis_.nCPs());
    const label vNCPs(vBasis_.nCPs());
    const label uCPI(CPsUCPIs_[CPI]);
    const label vCPI(CPsVCPIs_[CPI]);
    vector NMWP(Zero);
    scalar NMW(Zero);

    for (label vCPJ=0; vCPJ<vNCPs; vCPJ++)
    {
        for (label uCPJ=0; uCPJ<uNCPs; uCPJ++)
        {
            const label  CPJ(vCPJ*uNCPs + uCPJ);
            const scalar uBasisValue(uBasis_.basisValue(uCPJ, uDegree, u));
            const scalar vBasisValue(vBasis_.basisValue(vCPJ, vDegree, v));

            NMWP += uBasisValue * vBasisValue * weights_[CPJ] * CPs_[CPJ];
            NMW  += uBasisValue * vBasisValue * weights_[CPJ];
        }
    }

    // Compute derivative.
    const vector WDerivative
    (
          uBasis_.basisValue(uCPI, uDegree, u)
        * vBasis_.basisValue(vCPI, vDegree, v)
        * (CPs_[CPI] - (NMWP / (NMW+SMALL)))
        / (NMW+SMALL)
    );

    return WDerivative;
}


scalar NURBS3DSurface::lengthDerivativeU
(
    const scalar vConst,
    const scalar uStart,
    const scalar uEnd,
    const label nPts
) const
{
    // compute derivatives wrt u for the given u interval
    vectorField dxdu(nPts, Zero);
    vectorField d2xdu2(nPts, Zero);
    scalarField localU(nPts, Zero);
    scalar ulDerivative(Zero);

    forAll(localU, uI)
    {
        scalar& uLocal(localU[uI]);
        uLocal = uStart + scalar(uI)/scalar(nPts-1)*(uEnd-uStart);
        dxdu[uI] = surfaceDerivativeU(uLocal, vConst);
        d2xdu2[uI] = surfaceDerivativeUU(uLocal, vConst);
    }

    // Integrate.
    for (label uI=0; uI<(nPts-1); uI++)
    {
        ulDerivative +=
            0.5
         * (
               (dxdu[uI+1]&d2xdu2[uI+1])/(mag(dxdu[uI+1])+SMALL)
             + (dxdu[uI  ]&d2xdu2[uI  ])/(mag(dxdu[uI  ])+SMALL)
           )
         * (localU[uI+1]-localU[uI]);
    }

    return ulDerivative;
}


scalar NURBS3DSurface::lengthDerivativeV
(
    const scalar uConst,
    const scalar vStart,
    const scalar vEnd,
    const label nPts
) const
{
    // Compute derivatives wrt v for the given v interval.
    vectorField dxdv(nPts, Zero);
    vectorField d2xdv2(nPts, Zero);
    scalarField localV(nPts, Zero);
    scalar vlDerivative(Zero);

    forAll(localV, vI)
    {
        scalar& vLocal(localV[vI]);
        vLocal = vStart + scalar(vI)/scalar(nPts-1)*(vEnd-vStart);
        dxdv[vI] = surfaceDerivativeV(uConst, vLocal);
        d2xdv2[vI] = surfaceDerivativeVV(uConst, vLocal);
    }

    // Integrate.
    for (label vI=0; vI<(nPts-1); vI++)
    {
        vlDerivative +=
            0.5
         * (
               (dxdv[vI+1]&d2xdv2[vI+1])/(mag(dxdv[vI+1])+SMALL)
             + (dxdv[vI  ]&d2xdv2[vI  ])/(mag(dxdv[vI  ])+SMALL)
           )
         * (localV[vI+1]-localV[vI]);
    }

    return vlDerivative;
}


// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //

const labelList& NURBS3DSurface::getBoundaryCPIDs()
{
    if (!boundaryCPIDs_)
    {
        const label uNCPs(uBasis_.nCPs());
        const label vNCPs(vBasis_.nCPs());
        const label nBoundCPs(2*uNCPs+2*vNCPs-4);
        boundaryCPIDs_.reset(new labelList(nBoundCPs, -1));
        whichBoundaryCPID_.reset(new labelList(uNCPs*vNCPs, -1));

        // v-constant cps
        label bID(0);
        for (label vI=0; vI<vNCPs; vI+=vNCPs-1)
        {
            for (label uI=0; uI<uNCPs; uI++)
            {
                const label CPI(vI*uNCPs + uI);
                whichBoundaryCPID_()[CPI] = bID;
                boundaryCPIDs_()[bID++]   = CPI;
            }
        }
        // u-constant cps
        for (label uI=0; uI<uNCPs; uI+=uNCPs-1)
        {
            // corner CPS already accounted for
            for (label vI=1; vI<vNCPs-1; vI++)
            {
                const label CPI(vI*uNCPs + uI);
                whichBoundaryCPID_()[CPI] = bID;
                boundaryCPIDs_()[bID++]   = CPI;
            }
        }
    }

    return boundaryCPIDs_();
}


const labelList& NURBS3DSurface::getBoundaryCPIs()
{
    return getBoundaryCPIDs();
}


const label& NURBS3DSurface::whichBoundaryCPI(const label& globalCPI)
{
    if (!whichBoundaryCPID_)
    {
        getBoundaryCPIDs();
    }

    return whichBoundaryCPID_()[globalCPI];
}


// * * * * * * * * * * * * * * * Write Functions  * * * * * * * * * * * * * //

void NURBS3DSurface::write()
{
    write(name_);
}


void NURBS3DSurface::write(const word fileName)
{
    if (Pstream::master())
    {
        OFstream surfaceFile(fileName);
        OFstream surfaceFileCPs(fileName+"CPs");
        vectorField& surface(*this);

        forAll(*this, ptI)
        {
            surfaceFile
                << surface[ptI].x() << " "
                << surface[ptI].y() << " "
                << surface[ptI].z()
                << endl;
        }

        forAll(CPs_, CPI)
        {
            surfaceFileCPs
                << CPs_[CPI].x() << " "
                << CPs_[CPI].y() << " "
                << CPs_[CPI].z()
                << endl;
        }
/*
        const label  uDegree(uBasis_.degree());
        const label  vDegree(vBasis_.degree());
        const label  uNCPs(uBasis_.nCPs());
        const label  vNCPs(vBasis_.nCPs());

        OFstream     surfaceFileUBases(fileName+"UBases");
        OFstream     surfaceFileVBases(fileName+"VBases");

        forAll(*this, ptI)
        {
            const scalar& u(u_[ptI]);
            const scalar& v(v_[ptI]);
            scalarField   uBasesValues(uNCPs);
            scalarField   vBasesValues(vNCPs);

            surfaceFileUBases << u << " ";
            surfaceFileVBases << v << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
*/
    }
}


void NURBS3DSurface::write(const fileName dirName, const fileName fileName)
{
    if (Pstream::master())
    {
        OFstream surfaceFile(dirName/fileName);
        OFstream surfaceFileCPs(dirName/fileName+"CPs");
        vectorField& surface(*this);

        forAll(*this, ptI)
        {
            surfaceFile
                << surface[ptI].x() << " "
                << surface[ptI].y() << " "
                << surface[ptI].z()
                << endl;
        }

        forAll(CPs_, CPI)
        {
            surfaceFileCPs
                << CPs_[CPI].x() << " "
                << CPs_[CPI].y() << " "
                << CPs_[CPI].z()
                << endl;
        }
/*
        const label  uDegree(uBasis_.degree());
        const label  vDegree(vBasis_.degree());
        const label  uNCPs(uBasis_.nCPs());
        const label  vNCPs(vBasis_.nCPs());

        OFstream     surfaceFileUBases(dirName/fileName+"UBases");
        OFstream     surfaceFileVBases(dirName/fileName+"VBases");

        forAll(*this, ptI)
        {
            const scalar& u(u_[ptI]);
            const scalar& v(v_[ptI]);
            scalarField   uBasesValues(uNCPs);
            scalarField   vBasesValues(vNCPs);

            surfaceFileUBases << u << " ";
            surfaceFileVBases << v << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
*/
    }
}


void NURBS3DSurface::writeWParses()
{
    writeWParses(name_);
}


void NURBS3DSurface::writeWParses(const word fileName)
{
    if (Pstream::master())
    {
        OFstream surfaceFile(fileName);
        OFstream surfaceFileCPs(fileName+"CPs");
        vectorField& surface(*this);

        forAll(*this, ptI)
        {
            surfaceFile
                << "("
                << surface[ptI].x() << " "
                << surface[ptI].y() << " "
                << surface[ptI].z() << ")"
                << endl;
        }

        forAll(CPs_, CPI)
        {
            surfaceFileCPs
                << "("
                << CPs_[CPI].x() << " "
                << CPs_[CPI].y() << " "
                << CPs_[CPI].z() << ")"
                << endl;
        }
/*
        const label  uDegree(uBasis_.degree());
        const label  vDegree(vBasis_.degree());
        const label  uNCPs(uBasis_.nCPs());
        const label  vNCPs(vBasis_.nCPs());

        OFstream     surfaceFileUBases(fileName+"UBases");
        OFstream     surfaceFileVBases(fileName+"VBases");

        forAll(*this, ptI)
        {
            const scalar& u(u_[ptI]);
            const scalar& v(v_[ptI]);
            scalarField   uBasesValues(uNCPs);
            scalarField   vBasesValues(vNCPs);

            surfaceFileUBases << u << " ";
            surfaceFileVBases << v << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
*/
    }
}


void NURBS3DSurface::writeWParses
(
    const fileName dirName,
    const fileName fileName
)
{
    if (Pstream::master())
    {
        OFstream surfaceFile(dirName/fileName);
        OFstream surfaceFileCPs(dirName/fileName+"CPs");
        vectorField& surface(*this);

        forAll(*this, ptI)
        {
            surfaceFile
                << "("
                << surface[ptI].x() << " "
                << surface[ptI].y() << " "
                << surface[ptI].z() << ")"
                << endl;
        }

        forAll(CPs_, CPI)
        {
            surfaceFileCPs
                << "("
                << CPs_[CPI].x() << " "
                << CPs_[CPI].y() << " "
                << CPs_[CPI].z() << ")"
                << endl;
        }
/*
        const label  uDegree(uBasis_.degree());
        const label  vDegree(vBasis_.degree());
        const label  uNCPs(uBasis_.nCPs());
        const label  vNCPs(vBasis_.nCPs());

        OFstream     surfaceFileUBases(dirName/fileName+"UBases");
        OFstream     surfaceFileVBases(dirName/fileName+"VBases");

        forAll(*this, ptI)
        {
            const scalar& u(u_[ptI]);
            const scalar& v(v_[ptI]);
            scalarField   uBasesValues(uNCPs);
            scalarField   vBasesValues(vNCPs);

            surfaceFileUBases << u << " ";
            surfaceFileVBases << v << " ";

            forAll(CPs_, CPI)
            {
                basesValues[CPI] = basis_.basisValue(CPI, degree, u);
                curveFileBases <<  basesValues[CPI] << " ";
            }

            curveFileBases << endl;
        }
*/
    }
}


void NURBS3DSurface::writeVTK
(
    const fileName vtkDirName,
    const fileName vtkFileName
)
{
    if (Pstream::master())
    {
        if (vtkFileName.has_ext())
        {
            FatalErrorInFunction
                << "Do not supply a file extension."
                << exit(FatalError);
        }

        // Build the surface.
        buildSurface();

        // Write the surface as vtk.
        OFstream surfaceFile(vtkFileName);
        pointField& surfacePoints(*this);
        faceList surfaceFaces((nUPts_ - 1)*(nUPts_ - 1), face(4));

        for (label fuI = 0; fuI < (nUPts_ - 1); fuI++)
        {
            for (label fvI = 0; fvI < (nVPts_ - 1); fvI++)
            {
                const label fI(fuI*(nUPts_ - 1) + fvI);
                face& surfaceFace(surfaceFaces[fI]);

                surfaceFace[0] = (fuI)*nVPts_ + (fvI);
                surfaceFace[1] = (fuI + 1)*nVPts_ + (fvI);
                surfaceFace[2] = (fuI + 1)*nVPts_ + (fvI + 1);
                surfaceFace[3] = (fuI)*nVPts_ + (fvI + 1);
            }
        }

        surfaceWriters::vtkWriter writer;

        writer.open
        (
            surfacePoints,
            surfaceFaces,
            vtkDirName/vtkFileName,
            false
        );

        writer.close();

        // Write the control surface as vtk.
        fileName vtkCPFileName(vtkFileName+"CPs");
        pointField surfaceCPPoints(CPs_);
        const label uNCPs(uBasis_.nCPs());
        const label vNCPs(vBasis_.nCPs());
        faceList surfaceCPFaces((uNCPs-1)*(vNCPs-1), face(4));

        for (label fvCPI=0; fvCPI<(vNCPs-1); fvCPI++)
        {
            for (label fuCPI=0; fuCPI<(uNCPs-1); fuCPI++)
            {
                const label fCPI(fvCPI*(uNCPs-1) + fuCPI);
                face& surfaceCPFace(surfaceCPFaces[fCPI]);

                surfaceCPFace[0] = (fvCPI)*uNCPs + (fuCPI);
                surfaceCPFace[1] = (fvCPI + 1)*uNCPs + (fuCPI);
                surfaceCPFace[2] = (fvCPI + 1)*uNCPs + (fuCPI + 1);
                surfaceCPFace[3] = (fvCPI)*uNCPs + (fuCPI + 1);
            }
        }

        writer.open
        (
            surfaceCPPoints,
            surfaceCPFaces,
            vtkDirName/vtkCPFileName,
            false
        );

        writer.close();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
