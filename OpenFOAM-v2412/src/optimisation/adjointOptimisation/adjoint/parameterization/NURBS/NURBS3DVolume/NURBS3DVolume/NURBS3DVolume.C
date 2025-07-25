/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "NURBS3DVolume.H"

#include "OFstream.H"
#include "Time.H"
#include "deltaBoundary.H"
#include "coupledFvPatch.H"
#include "controlPointsDefinition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NURBS3DVolume, 0);
    defineRunTimeSelectionTable(NURBS3DVolume, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::NURBS3DVolume::findPointsInBox(const vectorField& meshPoints)
{
    // It is considered an error to recompute points in the control boxes
    if (mapPtr_ || reverseMapPtr_)
    {
        FatalErrorInFunction
           << "Attempting to recompute points residing within control boxes"
           << exit(FatalError);
    }

    mapPtr_.reset(new labelList(meshPoints.size(), -1));
    reverseMapPtr_.reset(new labelList(meshPoints.size(), -1));
    labelList& map = mapPtr_();
    labelList& reverseMap = reverseMapPtr_();

    // Identify points inside morphing boxes
    scalar lowerX = min(cps_.component(0));
    scalar upperX = max(cps_.component(0));
    scalar lowerY = min(cps_.component(1));
    scalar upperY = max(cps_.component(1));
    scalar lowerZ = min(cps_.component(2));
    scalar upperZ = max(cps_.component(2));

    Info<< "Control Points bounds \n"
        << "\tX1 : (" << lowerX << " " << upperX <<  ")\n"
        << "\tX2 : (" << lowerY << " " << upperY <<  ")\n"
        << "\tX3 : (" << lowerZ << " " << upperZ <<  ")\n" << endl;

    label count(0);
    forAll(meshPoints, pI)
    {
        const vector& pointI = meshPoints[pI];
        if
        (
            pointI.x() >= lowerX && pointI.x() <= upperX
         && pointI.y() >= lowerY && pointI.y() <= upperY
         && pointI.z() >= lowerZ && pointI.z() <= upperZ
        )
        {
            map[count] = pI;
            reverseMap[pI] = count;
            ++count;
        }
    }

    // Resize lists
    map.setSize(count);

    reduce(count, sumOp<label>());
    Info<< "Initially found " << count << " points inside control boxes"
        << endl;
}


void Foam::NURBS3DVolume::computeParametricCoordinates
(
    const vectorField& points
)
{
    scalar timeBef = mesh_.time().elapsedCpuTime();

    if (parametricCoordinatesPtr_)
    {
        FatalErrorInFunction
           << "Attempting to recompute parametric coordinates"
           << exit(FatalError);
    }

    labelList& map = mapPtr_();
    labelList& reverseMap = reverseMapPtr_();

    parametricCoordinatesPtr_.reset
    (
        new pointVectorField
        (
            IOobject
            (
                "parametricCoordinates" + name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh_),
            dimensionedVector(dimless, Zero)
        )
    );
    vectorField& paramCoors = parametricCoordinatesPtr_().primitiveFieldRef();

    // If already present, read from file
    if
    (
        parametricCoordinatesPtr_().typeHeaderOk<pointVectorField>(true)
     && readStoredData_
    )
    {
        // These are the parametric coordinates that have been computed
        // correctly (hopefully!) in a previous run.
        // findPointsInBox has possibly found more points and lists need
        // resizing
        Info<< "Reading parametric coordinates from file" << endl;
        IOobject& header = parametricCoordinatesPtr_().ref();
        parametricCoordinatesPtr_() =
            pointVectorField
            (
                header,
                pointMesh::New(mesh_)
            );

        // Initialize intermediate fields with sizes from findPointInBox
        labelList actualMap(map.size());

        // Read and store
        label curIndex(0);
        forAll(map, pI)
        {
            const label globalPointIndex = map[pI];
            if (paramCoors[globalPointIndex] != vector::zero)
            {
                actualMap[curIndex] = map[pI];
                reverseMap[globalPointIndex] = curIndex;
                ++curIndex;
            }
            else
            {
                reverseMap[globalPointIndex] = -1;
            }
        }

        // Resize intermediate
        actualMap.setSize(curIndex);

        reduce(curIndex, sumOp<label>());
        Info<< "Read non-zero parametric coordinates for " << curIndex
            << " points" << endl;

        // Update lists with the appropriate entries
        map = actualMap;
    }
    // Else, compute parametric coordinates iteratively
    else
    {
        // Initialize parametric coordinates based on min/max of control points
        scalar minX1 = min(cps_.component(0));
        scalar maxX1 = max(cps_.component(0));
        scalar minX2 = min(cps_.component(1));
        scalar maxX2 = max(cps_.component(1));
        scalar minX3 = min(cps_.component(2));
        scalar maxX3 = max(cps_.component(2));

        scalar oneOverDenomX(1./(maxX1 - minX1));
        scalar oneOverDenomY(1./(maxX2 - minX2));
        scalar oneOverDenomZ(1./(maxX3 - minX3));

        forAll(map, pI)
        {
            const label globalPI = map[pI];
            paramCoors[globalPI].x() = (points[pI].x() - minX1)*oneOverDenomX;
            paramCoors[globalPI].y() = (points[pI].y() - minX2)*oneOverDenomY;
            paramCoors[globalPI].z() = (points[pI].z() - minX3)*oneOverDenomZ;
        }

        // Indices of points that failed to converge
        // (i.e. are bounded for nMaxBound iters)
        boolList dropOffPoints(map.size(), false);
        label nDropedPoints(0);

        // Initial cartesian coordinates
        tmp<vectorField> tsplinesBasedCoors(coordinates(paramCoors));
        vectorField& splinesBasedCoors = tsplinesBasedCoors.ref();

        // Newton-Raphson loop to compute parametric coordinates
        // based on cartesian coordinates and the known control points
        Info<< "Mapping of mesh points to parametric space for box " << name_
            << " ..." << endl;
        // Do loop on a point-basis and check residual of each point equation
        label maxIterNeeded(0);
        forAll(points, pI)
        {
            label iter(0);
            label nBoundIters(0);
            vector res(GREAT, GREAT, GREAT);
            do
            {
                const label globalPI = map[pI];
                vector& uVec = paramCoors[globalPI];
                vector& coorPointI = splinesBasedCoors[pI];
                uVec += ((inv(JacobianUVW(uVec))) & (points[pI] - coorPointI));
                // Bounding might be needed for the first iterations
                // If multiple bounds happen, point is outside of the control
                // boxes and should be discarded
                if (bound(uVec))
                {
                    ++nBoundIters;
                }
                if (nBoundIters > nMaxBound_)
                {
                    dropOffPoints[pI] = true;
                    ++nDropedPoints;
                    break;
                }
                // Update current cartesian coordinates based on parametric ones
                coorPointI = coordinates(uVec);
                // Compute residual
                res = cmptMag(points[pI] - coorPointI);
            }
            while
            (
                (iter++ < maxIter_)
             && (
                       res.component(0) > tolerance_
                    || res.component(1) > tolerance_
                    || res.component(2) > tolerance_
                )
            );
            if (iter > maxIter_)
            {
                WarningInFunction
                    << "Mapping to parametric space for point " << pI
                    << " failed." << endl
                    << "Residual after " << maxIter_ + 1 << " iterations : "
                    << res << endl
                    << "parametric coordinates " <<  paramCoors[map[pI]]
                    <<  endl
                    << "Local system coordinates " <<  points[pI] <<  endl
                    << "Threshold residual per direction : " << tolerance_
                    << endl;
            }
            maxIterNeeded = max(maxIterNeeded, iter);
        }
        reduce(maxIterNeeded, maxOp<label>());

        label nParameterizedPoints = map.size() - nDropedPoints;

        // Resize mapping lists and parametric coordinates after dropping some
        // points
        labelList mapOld(map);

        map.setSize(nParameterizedPoints);

        label curIndex(0);
        forAll(dropOffPoints, pI)
        {
            if (!dropOffPoints[pI])
            {
                map[curIndex] = mapOld[pI];
                reverseMap[mapOld[pI]] = curIndex;
                ++curIndex;
            }
            else
            {
                paramCoors[mapOld[pI]] = vector::zero;
                reverseMap[mapOld[pI]] = -1;
            }
        }

        reduce(nDropedPoints, sumOp<label>());
        reduce(nParameterizedPoints, sumOp<label>());
        Info<< "Found " << nDropedPoints
            << " to discard from morphing boxes" << endl;
        Info<< "Keeping " << nParameterizedPoints
            << " parameterized points in boxes"  << endl;

        splinesBasedCoors = coordinates(paramCoors)();
        scalar maxDiff(-GREAT);
        forAll(splinesBasedCoors, pI)
        {
            scalar diff =
                mag(splinesBasedCoors[pI] - localSystemCoordinates_[map[pI]]);
            if (diff > maxDiff)
            {
                maxDiff = diff;
            }
        }
        reduce(maxDiff, maxOp<scalar>());
        scalar timeAft = mesh_.time().elapsedCpuTime();
        Info<< "\tMapping completed in " << timeAft - timeBef << " seconds"
            << endl;
        Info<< "\tMax iterations per point needed to compute parametric "
            << "coordinates : "
            << maxIterNeeded << endl;
        Info<< "\tMax difference between original mesh points and "
            << "parameterized ones "
            << maxDiff << endl;
    }
}


void Foam::NURBS3DVolume::computeParametricCoordinates
(
    tmp<vectorField> tPoints
)
{
    const vectorField& points = tPoints();
    computeParametricCoordinates(points);
}


bool Foam::NURBS3DVolume::bound
(
    vector& vec,
    scalar minValue,
    scalar maxValue
) const
{
    bool boundPoint(false);
    // Lower value bounding
    if (vec.x() < scalar(0))
    {
        vec.x() = minValue;
        boundPoint = true;
    }
    if (vec.y() < scalar(0))
    {
        vec.y() = minValue;
        boundPoint = true;
    }
    if (vec.z() < scalar(0))
    {
        vec.z() = minValue;
        boundPoint = true;
    }
    // Upper value bounding
    if (vec.x() > 1)
    {
        vec.x() = maxValue;
        boundPoint = true;
    }
    if (vec.y() > 1)
    {
         vec.y() = maxValue;
        boundPoint = true;
    }
    if (vec.z() > 1)
    {
        vec.z() = maxValue;
        boundPoint = true;
    }
    return boundPoint;
}


void Foam::NURBS3DVolume::makeFolders()
{
    if (Pstream::master())
    {
        mkDir(mesh_.time().globalPath()/"optimisation"/cpsFolder_);
    }
}


void Foam::NURBS3DVolume::determineActiveDesignVariablesAndPoints()
{
    label nCPs = cps_.size();
    activeControlPoints_ = boolList(nCPs, true);
    activeDesignVariables_ = boolList(3*nCPs, true);

    // Check whether all boundary control points should be confined
    confineBoundaryControlPoints();

    // Apply confinement to maintain continuity
    continuityRealatedConfinement();

    // Confine user-specified directions
    confineControlPointsDirections();

    // Determine active control points. A control point is considered active
    // if at least one of its components is free to move
    forAll(activeControlPoints_, cpI)
    {
        if
        (
            !activeDesignVariables_[3*cpI]
         && !activeDesignVariables_[3*cpI + 1]
         && !activeDesignVariables_[3*cpI + 2]
        )
        {
            activeControlPoints_[cpI] = false;
        }
    }

    // Deactivate design varibles if they affect no mesh point
    confineInertControlPoints();
}


void Foam::NURBS3DVolume::confineInertControlPoints()
{
    // Loop over all active control points and check whether they parameterize
    // at least one mesh point of the ones laying within the morphing box
    const labelList& map = getMap();
    const pointVectorField& parametricCoors = getParametricCoordinates();

    const scalarField& knotsU = basisU_.knots();
    const scalarField& knotsV = basisV_.knots();
    const scalarField& knotsW = basisW_.knots();

    const label pU = basisU_.degree();
    const label pV = basisV_.degree();
    const label pW = basisW_.degree();
    forAll(activeControlPoints_, cpI)
    {
        if (activeControlPoints_[cpI])
        {
            // Get i, j, k corresponding to this control point
            label i(-1),  j(-1), k(-1);
            getIJK(i, j, k, cpI);

            const scalar uMin(knotsU[i]);
            const scalar uMax(knotsU[i + pU + 1]);
            const scalar vMin(knotsV[j]);
            const scalar vMax(knotsV[j + pV + 1]);
            const scalar wMin(knotsW[k]);
            const scalar wMax(knotsW[k + pW + 1]);
            bool foundParamPt(false);
            for (const label paramI : map)
            {
                const vector& paramCoors = parametricCoors()[paramI];
                if
                (
                    paramCoors.x() >= uMin && paramCoors.x() < uMax
                 && paramCoors.y() >= vMin && paramCoors.y() < vMax
                 && paramCoors.z() >= wMin && paramCoors.z() < wMax
                )
                {
                    foundParamPt = true;
                    break;
                }
            }
            reduce(foundParamPt, orOp<bool>());
            if (!foundParamPt)
            {
                activeControlPoints_[cpI] = false;
                activeDesignVariables_[3*cpI] = false;
                activeDesignVariables_[3*cpI + 1] = false;
                activeDesignVariables_[3*cpI + 2] = false;
                Info<< "Disabling control " << cpI
                    << " and variables "
                    << "(" << 3*cpI << ", " << 3*cpI+1 << ", " << 3*cpI+2 << ")"
                    << " since they does not parameterize any mesh point"
                    << endl;
            }
        }
    }
}


void Foam::NURBS3DVolume::confineBoundaryControlPoints()
{
    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    // Zero movement of the boundary control points. Active by default
    if (confineBoundaryControlPoints_)
    {
        // Side patches
        for (label iCPw = 0; iCPw < nCPsW; iCPw += nCPsW - 1)
        {
            for (label iCPv = 0; iCPv < nCPsV; iCPv++)
            {
                for (label iCPu = 0; iCPu < nCPsU; iCPu++)
                {
                    confineControlPoint(getCPID(iCPu, iCPv, iCPw));
                }
            }
        }
        // Front-back patches
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPv = 0; iCPv < nCPsV; iCPv++)
            {
                for (label iCPu = 0; iCPu < nCPsU; iCPu += nCPsU - 1)
                {
                    confineControlPoint(getCPID(iCPu, iCPv, iCPw));
                }
            }
        }
        // Top-bottom patches
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPv = 0; iCPv < nCPsV; iCPv += nCPsV - 1)
            {
                for (label iCPu = 0; iCPu < nCPsU; iCPu++)
                {
                    confineControlPoint(getCPID(iCPu, iCPv, iCPw));
                }
            }
        }
    }
}


void Foam::NURBS3DVolume::continuityRealatedConfinement()
{
    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    // Zero movement to a number of x-constant slices of cps in order to
    // preserve continuity at the boundary of the parameterized space
    forAll(confineUMinCPs_, iCPu)
    {
        const boolVector& confineSlice = confineUMinCPs_[iCPu];
        // Control points at the start of the parameterized space
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPv = 0; iCPv < nCPsV; iCPv++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }

    forAll(confineUMaxCPs_, sliceI)
    {
        const boolVector& confineSlice = confineUMaxCPs_[sliceI];
        label iCPu = nCPsU - 1 - sliceI;
        // Control points at the end of the parameterized space
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPv = 0; iCPv < nCPsV; iCPv++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }

    // Zero movement to a number of y-constant slices of cps in order to
    // preserve continuity at the boundary of the parameterized space
    forAll(confineVMinCPs_, iCPv)
    {
        const boolVector& confineSlice = confineVMinCPs_[iCPv];
        // Control points at the start of the parameterized space
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }

    forAll(confineVMaxCPs_, sliceI)
    {
        const boolVector& confineSlice = confineVMaxCPs_[sliceI];
        label iCPv = nCPsV - 1 - sliceI;
        // Control points at the end of the parameterized space
        for (label iCPw = 0; iCPw < nCPsW; iCPw++)
        {
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }

    // Zero movement to a number of w-constant slices of cps in order to
    // preserve continuity at the boundary of the parameterized space
    forAll(confineWMinCPs_, iCPw)
    {
        const boolVector& confineSlice = confineWMinCPs_[iCPw];
        // Control points at the start of the parameterized space
        for (label iCPv = 0; iCPv < nCPsV; iCPv++)
        {
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }

    forAll(confineWMaxCPs_, sliceI)
    {
        const boolVector& confineSlice = confineWMaxCPs_[sliceI];
        label iCPw = nCPsW - 1 - sliceI;
        // Control points at the end of the parameterized space
        for (label iCPv = 0; iCPv < nCPsV; iCPv++)
        {
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                confineControlPoint(getCPID(iCPu, iCPv, iCPw), confineSlice);
            }
        }
    }
}


void Foam::NURBS3DVolume::confineControlPointsDirections()
{
    for (label cpI = 0; cpI < cps_.size(); ++cpI)
    {
        if (confineUMovement_) activeDesignVariables_[3*cpI] = false;
        if (confineVMovement_) activeDesignVariables_[3*cpI + 1] = false;
        if (confineWMovement_) activeDesignVariables_[3*cpI + 2] = false;
    }
}


void Foam::NURBS3DVolume::confineControlPoint(const label cpI)
{
    if (cpI < 0 || cpI > cps_.size() -1)
    {
        FatalErrorInFunction
           << "Attempted to confine control point movement for a control point "
           << " ID which is out of bounds"
           << exit(FatalError);
    }
    else
    {
        activeDesignVariables_[3*cpI] = false;
        activeDesignVariables_[3*cpI + 1] = false;
        activeDesignVariables_[3*cpI + 2] = false;
    }
}


void Foam::NURBS3DVolume::confineControlPoint
(
    const label cpI,
    const boolVector& confineDirections
)
{
    if (cpI < 0 || cpI > cps_.size() -1)
    {
        FatalErrorInFunction
           << "Attempted to confine control point movement for a control point "
           << " ID which is out of bounds"
           << exit(FatalError);
    }
    else
    {
        if (confineDirections.x()) activeDesignVariables_[3*cpI] = false;
        if (confineDirections.y()) activeDesignVariables_[3*cpI + 1] = false;
        if (confineDirections.z()) activeDesignVariables_[3*cpI + 2] = false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NURBS3DVolume::NURBS3DVolume
(
    const dictionary& dict,
    const fvMesh& mesh,
    bool computeParamCoors
)
:
    localIOdictionary
    (
        IOobject
        (
            dict.dictName() + "cpsBsplines",
            mesh.time().timeName(),
            fileName("uniform")/fileName("volumetricBSplines"),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    mesh_(mesh),
    dict_(dict),
    name_(dict.dictName()),
    basisU_(dict.get<label>("nCPsU"), dict.get<label>("degreeU")),
    basisV_(dict.get<label>("nCPsV"), dict.get<label>("degreeV")),
    basisW_(dict.get<label>("nCPsW"), dict.get<label>("degreeW")),
    maxIter_(dict.getOrDefault<label>("maxIterations", 10)),
    tolerance_(dict.getOrDefault<scalar>("tolerance", 1.e-10)),
    nMaxBound_(dict.getOrDefault<scalar>("nMaxBoundIterations", 4)),
    cps_(),
    mapPtr_(nullptr),
    reverseMapPtr_(nullptr),
    parametricCoordinatesPtr_(nullptr),
    localSystemCoordinates_(mesh_.nPoints(), Zero),
    confineUMovement_
    (
        dict.getOrDefaultCompat<bool>
        (
            "confineUMovement", {{"confineX1movement", 1912}}, false
        )
    ),
    confineVMovement_
    (
        dict.getOrDefaultCompat<bool>
        (
            "confineVMovement", {{"confineX2movement", 1912}}, false
        )
    ),
    confineWMovement_
    (
        dict.getOrDefaultCompat<bool>
        (
            "confineWMovement", {{"confineX3movement", 1912}}, false
        )
    ),
    confineBoundaryControlPoints_
    (
        dict.getOrDefault<bool>("confineBoundaryControlPoints", true)
    ),
    confineUMinCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineUMinCPs", {{"boundUMinCPs", 1912}}, boolVectorList()
        )
    ),
    confineUMaxCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineUMaxCPs", {{"boundUMaxCPs", 1912}}, boolVectorList()
        )
    ),
    confineVMinCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineVMinCPs", {{"boundVMinCPs", 1912}}, boolVectorList()
        )
    ),
    confineVMaxCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineVMaxCPs", {{"boundVMaxCPs", 1912}}, boolVectorList()
        )
    ),
    confineWMinCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineWMinCPs", {{"boundWMinCPs", 1912}}, boolVectorList()
        )
    ),
    confineWMaxCPs_
    (
        dict.getOrDefaultCompat<boolVectorList>
        (
            "confineWMaxCPs", {{"boundWMaxCPs", 1912}}, boolVectorList()
        )
    ),
    activeControlPoints_(0), //zero here, execute sanity checks first
    activeDesignVariables_(0), //zero here, execute sanity checks first
    cpsFolder_("controlPoints"),
    readStoredData_(dict.getOrDefault<bool>("readStoredData", true))
{
    // Create folders
    makeFolders();

    // Sanity checks
    if
    (
        (confineUMinCPs_.size() + confineUMaxCPs_.size() >= basisU_.nCPs())
     || (confineVMinCPs_.size() + confineVMaxCPs_.size() >= basisV_.nCPs())
     || (confineWMinCPs_.size() + confineWMaxCPs_.size() >= basisW_.nCPs())
    )
    {
        FatalErrorInFunction
           << "Number of control point slices to be kept frozen at "
           << "the boundaries is invalid \n"
           << "Number of control points in u " << basisU_.nCPs() << "\n"
           << "Number of control points in v " << basisV_.nCPs() << "\n"
           << "Number of control points in w " << basisW_.nCPs() << "\n"
           << exit(FatalError);
    }

    // Construct control points, if not already read from file
    if (found("controlPoints"))
    {
        cps_ =
            vectorField
            (
                "controlPoints",
                *this,
                basisU_.nCPs()*basisV_.nCPs()*basisW_.nCPs()
            );
    }
    else
    {
        controlPointsDefinition::New(*this);
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::NURBS3DVolume> Foam::NURBS3DVolume::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    bool computeParamCoors
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "NURBS3DVolume type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "type",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<NURBS3DVolume>(ctorPtr(dict, mesh, computeParamCoors));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::vectorField> Foam::NURBS3DVolume::computeParametricCoordinates
(
    const vectorField& points
) const
{
    tmp<vectorField> tparamCoors(tmp<vectorField>::New(points.size(), Zero));
    vectorField& paramCoors = tparamCoors.ref();
    // Initialize parametric coordinates based on min/max of control points
    scalar minX1 = min(cps_.component(0));
    scalar maxX1 = max(cps_.component(0));
    scalar minX2 = min(cps_.component(1));
    scalar maxX2 = max(cps_.component(1));
    scalar minX3 = min(cps_.component(2));
    scalar maxX3 = max(cps_.component(2));

    scalar oneOverDenomX(1./(maxX1 - minX1));
    scalar oneOverDenomY(1./(maxX2 - minX2));
    scalar oneOverDenomZ(1./(maxX3 - minX3));

    forAll(points, pI)
    {
        paramCoors[pI].x() = (points[pI].x() - minX1)*oneOverDenomX;
        paramCoors[pI].y() = (points[pI].y() - minX2)*oneOverDenomY;
        paramCoors[pI].z() = (points[pI].z() - minX3)*oneOverDenomZ;
    }

    // Indices of points that failed to converge
    // (i.e. are bounded for nMaxBound iters)
    boolList dropOffPoints(points.size(), false);
    label nDropedPoints(0);

    // Initial cartesian coordinates
    tmp<vectorField> tsplinesBasedCoors(coordinates(paramCoors));
    vectorField& splinesBasedCoors = tsplinesBasedCoors.ref();

    // Newton-Raphson loop to compute parametric coordinates
    // based on cartesian coordinates and the known control points
    Info<< "Mapping of mesh points to parametric space for box " << name_
        << " ..." << endl;
    // Do loop on a point-basis and check residual of each point equation
    label maxIterNeeded(0);
    forAll(points, pI)
    {
        label iter(0);
        label nBoundIters(0);
        vector res(GREAT, GREAT, GREAT);
        do
        {
            vector& uVec = paramCoors[pI];
            vector& coorPointI = splinesBasedCoors[pI];
            uVec += ((inv(JacobianUVW(uVec))) & (points[pI] - coorPointI));
            // Bounding might be needed for the first iterations
            // If multiple bounds happen, point is outside of the control
            // boxes and should be discarded
            if (bound(uVec))
            {
                ++nBoundIters;
            }
            if (nBoundIters > nMaxBound_)
            {
                dropOffPoints[pI] = true;
                ++nDropedPoints;
                break;
            }
            // Update current cartesian coordinates based on parametric ones
            coorPointI = coordinates(uVec);
            // Compute residual
            res = cmptMag(points[pI] - coorPointI);
        }
        while
        (
            (iter++ < maxIter_)
         && (
                   res.component(0) > tolerance_
                || res.component(1) > tolerance_
                || res.component(2) > tolerance_
            )
        );
        if (iter > maxIter_)
        {
            WarningInFunction
                << "Mapping to parametric space for point " << pI
                << " failed." << endl
                << "Residual after " << maxIter_ + 1 << " iterations : "
                << res << endl
                << "parametric coordinates " <<  paramCoors[pI]
                <<  endl
                << "Local system coordinates " <<  points[pI] <<  endl
                << "Threshold residual per direction : " << tolerance_
                << endl;
        }
        maxIterNeeded = max(maxIterNeeded, iter);
    }
    reduce(maxIterNeeded, maxOp<label>());

    label nParameterizedPoints = points.size() - nDropedPoints;

    reduce(nDropedPoints, sumOp<label>());
    reduce(nParameterizedPoints, sumOp<label>());
    Info<< "Found " << nDropedPoints
        << " to discard from morphing boxes" << endl;
    Info<< "Keeping " << nParameterizedPoints
        << " parameterized points in boxes"  << endl;
    return tparamCoors;
}


Foam::vector Foam::NURBS3DVolume::volumeDerivativeU
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    const label degreeU = basisU_.degree();
    const label degreeV = basisV_.degree();
    const label degreeW = basisW_.degree();

    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    vector derivative(Zero);

    for (label iCPw = 0; iCPw < nCPsW; ++iCPw)
    {
        const scalar basisW(basisW_.basisValue(iCPw, degreeW, w));
        for (label iCPv = 0; iCPv < nCPsV; ++iCPv)
        {
            const scalar basisVW = basisW*basisV_.basisValue(iCPv, degreeV, v);
            for (label iCPu = 0; iCPu < nCPsU; ++iCPu)
            {
                derivative +=
                    cps_[getCPID(iCPu, iCPv, iCPw)]
                   *basisU_.basisDerivativeU(iCPu, degreeU, u)
                   *basisVW;
            }
        }
    }

    return derivative;
}


Foam::vector Foam::NURBS3DVolume::volumeDerivativeV
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    const label degreeU = basisU_.degree();
    const label degreeV = basisV_.degree();
    const label degreeW = basisW_.degree();

    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    vector derivative(Zero);

    for (label iCPw = 0; iCPw < nCPsW; ++iCPw)
    {
        const scalar basisW(basisW_.basisValue(iCPw, degreeW, w));
        for (label iCPv = 0; iCPv < nCPsV; ++iCPv)
        {
            const scalar basisWDeriV =
                basisW*basisV_.basisDerivativeU(iCPv, degreeV, v);
            for (label iCPu = 0; iCPu < nCPsU; ++iCPu)
            {
                derivative +=
                    cps_[getCPID(iCPu, iCPv, iCPw)]
                   *basisU_.basisValue(iCPu, degreeU, u)
                   *basisWDeriV;
            }
        }
    }

    return derivative;
}


Foam::vector Foam::NURBS3DVolume::volumeDerivativeW
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    const label degreeU = basisU_.degree();
    const label degreeV = basisV_.degree();
    const label degreeW = basisW_.degree();

    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    vector derivative(Zero);

    for (label iCPw = 0; iCPw < nCPsW; iCPw++)
    {
        const scalar derivW(basisW_.basisDerivativeU(iCPw, degreeW, w));
        for (label iCPv = 0; iCPv < nCPsV; iCPv++)
        {
            const scalar derivWBasisV =
                derivW*basisV_.basisValue(iCPv, degreeV, v);
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                derivative +=
                    cps_[getCPID(iCPu, iCPv, iCPw)]
                   *basisU_.basisValue(iCPu, degreeU, u)
                   *derivWBasisV;
            }
        }
    }

    return derivative;
}


Foam::tensor Foam::NURBS3DVolume::JacobianUVW
(
    const vector& uVector
) const
{
    const scalar u = uVector.x();
    const scalar v = uVector.y();
    const scalar w = uVector.z();

    vector uDeriv = volumeDerivativeU(u, v, w);
    vector vDeriv = volumeDerivativeV(u, v, w);
    vector wDeriv = volumeDerivativeW(u, v, w);

    tensor Jacobian(uDeriv, vDeriv, wDeriv, true);

    return Jacobian;
}


Foam::scalar Foam::NURBS3DVolume::volumeDerivativeCP
(
    const vector& uVector,
    const label cpI
) const
{
    const scalar u = uVector.x();
    const scalar v = uVector.y();
    const scalar w = uVector.z();

    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();

    const label degreeU = basisU_.degree();
    const label degreeV = basisV_.degree();
    const label degreeW = basisW_.degree();

    label iCPw = cpI/label(nCPsU*nCPsV);
    label iCPv = (cpI - iCPw*nCPsU*nCPsV)/nCPsU;
    label iCPu = (cpI - iCPw*nCPsU*nCPsV - iCPv*nCPsU);

    // Normally, this should be a tensor, however the parameterization is
    // isotropic. Hence the tensor degenerates to a diagonal tensor with all
    // diagonal elements being equal. This returns the (unique) diag element
    scalar derivative =
        basisU_.basisValue(iCPu, degreeU, u)
       *basisV_.basisValue(iCPv, degreeV, v)
       *basisW_.basisValue(iCPw, degreeW, w);

    return derivative;
}


Foam::vectorField Foam::NURBS3DVolume::computeControlPointSensitivities
(
    const pointVectorField& pointSens,
    const labelList& sensitivityPatchIDs
)
{
    vectorField controlPointDerivs(cps_.size(), Zero);

    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    forAll(controlPointDerivs, cpI)
    {
        forAll(sensitivityPatchIDs, pI)
        {
            const label patchI = sensitivityPatchIDs[pI];
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            const labelList& meshPoints = patch.meshPoints();

            forAll(meshPoints, mpI)
            {
                const label globalIndex = meshPoints[mpI];
                const label whichPointInBox = reverseMapPtr_()[globalIndex];

                // If point resides within control points box,
                // add contribution to cp derivative
                if (whichPointInBox != -1)
                {
                    controlPointDerivs[cpI] +=
                    (
                        pointSens[globalIndex]
                      & transformationTensorDxDb(globalIndex)
                    )
                   *volumeDerivativeCP
                    (
                        parametricCoordinates[globalIndex],
                        cpI
                    );
                }
            }
        }
    }

    // Sum contributions from all processors
    Pstream::listCombineReduce(controlPointDerivs, plusEqOp<vector>());

    return controlPointDerivs;
}


Foam::vectorField Foam::NURBS3DVolume::computeControlPointSensitivities
(
    const volVectorField& faceSens,
    const labelList& sensitivityPatchIDs
)
{
    return
        computeControlPointSensitivities
        (
            faceSens.boundaryField(),
            sensitivityPatchIDs
        );
}


Foam::vectorField Foam::NURBS3DVolume::computeControlPointSensitivities
(
    const boundaryVectorField& faceSens,
    const labelList& sensitivityPatchIDs
)
{
    // Return field
    vectorField controlPointDerivs(cps_.size(), Zero);

    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Auxiliary quantities
    deltaBoundary deltaBoundary(mesh_);
    const labelList& reverseMap = reverseMapPtr_();

    forAll(controlPointDerivs, cpI)
    {
        forAll(sensitivityPatchIDs, pI)
        {
            const label patchI = sensitivityPatchIDs[pI];
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            const label patchStart = patch.start();
            const fvPatchVectorField& patchSens = faceSens[patchI];

            // loop over patch faces
            forAll(patch, fI)
            {
                const face& fGlobal = mesh_.faces()[fI + patchStart];
                const pointField facePoints = fGlobal.points(mesh_.points());
                // loop over face points
                tensorField facePointDerivs(facePoints.size(), Zero);
                forAll(fGlobal, pI)
                {
                    const label globalIndex = fGlobal[pI]; //global point index
                    const label whichPointInBox = reverseMap[globalIndex];
                    // if point resides within control points box,
                    // add contribution to d( facePoints )/db
                    if (whichPointInBox != -1)
                    {
                        // TENSOR-BASED
                        //~~~~~~~~~~~~~
                        facePointDerivs[pI] =
                            transformationTensorDxDb(globalIndex)
                          * volumeDerivativeCP
                            (
                                parametricCoordinates[globalIndex],
                                cpI
                            );

                    }
                }

                tensor fCtrs_d =
                    deltaBoundary.makeFaceCentresAndAreas_d
                    (
                        facePoints,
                        facePointDerivs
                    )[0];
                controlPointDerivs[cpI] += patchSens[fI] & fCtrs_d;
            }
        }
    }
    // Sum contributions from all processors
    Pstream::listCombineReduce(controlPointDerivs, plusEqOp<vector>());

    return controlPointDerivs;
}


Foam::vector Foam::NURBS3DVolume::computeControlPointSensitivities
(
    const vectorField& faceSens,
    const label patchI,
    const label cpI
)
{
    // Return vector
    vector cpSens(Zero);
    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Auxiliary quantities
    deltaBoundary deltaBoundary(mesh_);
    const labelList& reverseMap = reverseMapPtr_();

    const polyPatch& patch = mesh_.boundaryMesh()[patchI];
    const label patchStart = patch.start();
    // Loop over patch faces
    forAll(patch, fI)
    {
        const face& fGlobal = mesh_.faces()[fI + patchStart];
        const pointField facePoints = fGlobal.points(mesh_.points());
        // Loop over face points
        tensorField facePointDerivs(facePoints.size(), Zero);
        forAll(fGlobal, pI)
        {
            const label globalIndex = fGlobal[pI];  //global point index
            const label whichPointInBox = reverseMap[globalIndex];
            // If point resides within control points box,
            // add contribution to d( facePoints )/db
            if (whichPointInBox != -1)
            {
                // TENSOR-BASED
                //~~~~~~~~~~~~~
                facePointDerivs[pI] =
                    transformationTensorDxDb(globalIndex)
                   *volumeDerivativeCP
                    (
                        parametricCoordinates[globalIndex],
                        cpI
                    );
            }
        }

        tensor fCtrs_d =
            deltaBoundary.makeFaceCentresAndAreas_d
            (
                facePoints,
                facePointDerivs
            )[0];
        cpSens += faceSens[fI] & fCtrs_d;
    }
    // Sum contributions from all processors
    reduce(cpSens, sumOp<vector>());

    return cpSens;
}


Foam::tmp<Foam::tensorField> Foam::NURBS3DVolume::dndbBasedSensitivities
(
    const label patchI,
    const label cpI,
    bool DimensionedNormalSens
)
{
    const fvPatch& patch = mesh_.boundary()[patchI];
    const polyPatch& ppatch = patch.patch();
    // Return field
    auto tdndbSens = tmp<tensorField>::New(patch.size(), Zero);
    auto& dndbSens = tdndbSens.ref();

    // Auxiliary quantities
    deltaBoundary deltaBoundary(mesh_);
    const label patchStart = ppatch.start();
    const labelList& reverseMap = reverseMapPtr_();

    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Loop over patch faces
    forAll(patch, fI)
    {
        const face& fGlobal = mesh_.faces()[fI + patchStart];
        const pointField facePoints = fGlobal.points(mesh_.points());
        // Loop over face points
        tensorField facePointDerivs(facePoints.size(), Zero);
        forAll(fGlobal, pI)
        {
            const label globalIndex = fGlobal[pI];  //global point index
            const label whichPointInBox = reverseMap[globalIndex];
            // If point resides within control points box,
            // add contribution to d( facePoints )/db
            if (whichPointInBox != -1)
            {
                // TENSOR-BASED
                //~~~~~~~~~~~~~
                facePointDerivs[pI] =
                    transformationTensorDxDb(globalIndex)
                   *volumeDerivativeCP
                    (
                        parametricCoordinates[globalIndex],
                        cpI
                    );
            }
        }

        // Determine whether to return variance of dimensioned or unit normal
        tensorField dNdbSens =
            deltaBoundary.makeFaceCentresAndAreas_d
            (
                facePoints,
                facePointDerivs
            );

        if (DimensionedNormalSens)
        {
            dndbSens[fI] = dNdbSens[1];
        }
        else
        {
            dndbSens[fI] = dNdbSens[2];
        }
    }

    return tdndbSens;
}


Foam::tmp<Foam::tensorField> Foam::NURBS3DVolume::patchDxDb
(
    const label patchI,
    const label cpI
)
{
    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Patch data
    const polyPatch& patch = mesh_.boundaryMesh()[patchI];
    const labelList& meshPoints = patch.meshPoints();

    // Return field
    auto tdxdb = tmp<tensorField>::New(patch.nPoints(), Zero);
    auto& dxdb = tdxdb.ref();

    forAll(meshPoints, pI)
    {
        const label globalIndex = meshPoints[pI];  //global point index
        const label whichPointInBox = reverseMapPtr_()[globalIndex];

        // If point resides within control points box, find dxdb
        if (whichPointInBox != -1)
        {
            dxdb[pI] =
                transformationTensorDxDb(globalIndex)
               *volumeDerivativeCP
                (
                    parametricCoordinates[globalIndex],
                    cpI
                );
        }
    }

    return tdxdb;
}


Foam::tmp<Foam::tensorField> Foam::NURBS3DVolume::patchDxDbFace
(
    const label patchI,
    const label cpI
)
{
    // get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Patch data
    const polyPatch& patch = mesh_.boundaryMesh()[patchI];
    const label patchStart = patch.start();

    // Return field
    auto tdxdb = tmp<tensorField>::New(patch.size(), Zero);
    auto& dxdb = tdxdb.ref();

    // Mesh differentiation engine
    deltaBoundary deltaBound(mesh_);

    forAll(patch, fI)
    {
        const face& fGlobal = mesh_.faces()[fI + patchStart];
        const pointField facePoints = fGlobal.points(mesh_.points());
        // Loop over face points
        tensorField facePointDerivs(facePoints.size(), Zero);
        forAll(fGlobal, pI)
        {
            const label globalIndex = fGlobal[pI];  //global point index
            const label whichPointInBox = reverseMapPtr_()[globalIndex];
            // If point resides within control points box,
            // add contribution to d( facePoints )/db
            if (whichPointInBox != -1)
            {
                // TENSOR-BASED
                //~~~~~~~~~~~~~
                facePointDerivs[pI] =
                    transformationTensorDxDb(globalIndex)
                   *volumeDerivativeCP
                    (
                        parametricCoordinates[globalIndex],
                        cpI
                    );

            }
        }
        dxdb[fI] =
            deltaBound.makeFaceCentresAndAreas_d
            (
                facePoints,
                facePointDerivs
            )[0];
    }

    return tdxdb;
}


Foam::tmp<Foam::vectorField> Foam::NURBS3DVolume::coordinates
(
    const vectorField& uVector
) const
{
    const label nPoints = mapPtr_().size();
    auto tpoints = tmp<vectorField>::New(nPoints, Zero);
    auto& points = tpoints.ref();

    forAll(points, pI)
    {
        const label globalPI = mapPtr_()[pI];
        points[pI] = coordinates(uVector[globalPI]);
    }

    return tpoints;
}


Foam::vector Foam::NURBS3DVolume::coordinates
(
    const vector& uVector
) const
{
    const label degreeU = basisU_.degree();
    const label degreeV = basisV_.degree();
    const label degreeW = basisW_.degree();

    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    const label nCPsW = basisW_.nCPs();

    const scalar u = uVector.x();
    const scalar v = uVector.y();
    const scalar w = uVector.z();

    vector point(Zero);
    for (label iCPw = 0; iCPw < nCPsW; iCPw++)
    {
        const scalar basisW(basisW_.basisValue(iCPw, degreeW, w));
        for (label iCPv = 0; iCPv < nCPsV; iCPv++)
        {
            const scalar basisVW =
                basisW*basisV_.basisValue(iCPv, degreeV, v);
            for (label iCPu = 0; iCPu < nCPsU; iCPu++)
            {
                point +=
                   cps_[getCPID(iCPu, iCPv, iCPw)]
                  *basisU_.basisValue(iCPu, degreeU, u)
                  *basisVW;
            }
        }
    }

    return point;
}


Foam::tmp<Foam::vectorField> Foam::NURBS3DVolume::computeNewPoints
(
    const vectorField& controlPointsMovement
)
{
    // Get parametric coordinates and map
    const vectorField& paramCoors = getParametricCoordinates();
    const labelList& map = mapPtr_();

    // Update control points position
    cps_ += controlPointsMovement;
    writeCps("cpsBsplines"+mesh_.time().timeName());

    // Compute new mesh points based on updated control points
    tmp<vectorField> tparameterizedPoints = coordinates(paramCoors);
    const vectorField& parameterizedPoints = tparameterizedPoints();

    // Return field. Initialized with current mesh points
    auto tnewPoints = tmp<vectorField>::New(mesh_.points());
    auto& newPoints = tnewPoints.ref();

    // Update position of parameterized points
    forAll(parameterizedPoints, pI)
    {
        newPoints[map[pI]] = transformPointToCartesian(parameterizedPoints[pI]);
    }

    // Update coordinates in the local system based on the cartesian points
    updateLocalCoordinateSystem(newPoints);
    DebugInfo
        << "Max mesh movement equal to "
        << gMax(mag(newPoints - mesh_.points())) << endl;

    return tnewPoints;
}


Foam::tmp<Foam::vectorField> Foam::NURBS3DVolume::computeNewBoundaryPoints
(
    const vectorField& controlPointsMovement,
    const labelList& patchesToBeMoved,
    const bool updateCPs
)
{
    // Get parametric coordinates
    const vectorField& paramCoors = getParametricCoordinates();

    // Update control points position
    cps_ += controlPointsMovement;

    if (updateCPs)
    {
        writeCps("cpsBsplines"+mesh_.time().timeName());
    }

    // Return field. Initialized with current mesh points
    auto tnewPoints = tmp<vectorField>::New(mesh_.points());
    auto& newPoints = tnewPoints.ref();

    // Update position of parameterized boundary points
    for (const label patchI : patchesToBeMoved)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        const labelList& meshPoints = patch.meshPoints();

        for (const label globalIndex : meshPoints)
        {
            const label whichPointInBox = reverseMapPtr_()[globalIndex];
            // If point resides within control points box,
            // compute new cartesian coordinates
            if (whichPointInBox != -1)
            {
                newPoints[globalIndex] =
                    transformPointToCartesian
                    (
                        coordinates
                        (
                            paramCoors[globalIndex]
                        )
                    );
            }
        }
    }

    if (updateCPs)
    {
        // Update coordinates in the local system based on the cartesian points
        updateLocalCoordinateSystem(newPoints);
    }
    else
    {
        // Move control points to their initial position
        cps_ -= controlPointsMovement;
    }

    DebugInfo
        << "Max mesh movement equal to "
        << gMax(mag(newPoints - mesh_.points())) << endl;

    return tnewPoints;
}


Foam::label Foam::NURBS3DVolume::getCPID
(
    const label i,
    const label j,
    const label k
) const
{
    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();

    return k*nCPsU*nCPsV + j*nCPsU + i;
}


void Foam::NURBS3DVolume::getIJK
(
    label&  i,
    label&  j,
    label&  k,
    const label cpID
) const
{
    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();
    k = cpID/(nCPsU*nCPsV);
    const label inKplaneID = (cpID - k*(nCPsU*nCPsV));
    j = inKplaneID/nCPsU;
    i = inKplaneID - j*nCPsU;
}


void Foam::NURBS3DVolume::setControlPoints(const vectorField& newCps)
{
    if (cps_.size() != newCps.size())
    {
        FatalErrorInFunction
           << "Attempting to replace control points with a set of "
           << "different size"
           << exit(FatalError);
    }
    cps_ = newCps;
}


void Foam::NURBS3DVolume::boundControlPointMovement
(
    vectorField& controlPointsMovement
) const
{
    forAll(controlPointsMovement, cpI)
    {
        if (!activeDesignVariables_[3*cpI])
        {
            controlPointsMovement[cpI].x() = Zero;
        }
        if (!activeDesignVariables_[3*cpI + 1])
        {
            controlPointsMovement[cpI].y() = Zero;
        }
        if (!activeDesignVariables_[3*cpI + 2])
        {
            controlPointsMovement[cpI].z() = Zero;
        }
    }
}


Foam::scalar Foam::NURBS3DVolume::computeMaxBoundaryDisplacement
(
    const vectorField& controlPointsMovement,
    const labelList& patchesToBeMoved
)
{
    // Backup old cps
    vectorField oldCPs = cps_;
    // Get parametric coordinates
    const vectorField& paramCoors = getParametricCoordinates();
    // Update control points position
    cps_ += controlPointsMovement;
    // Update position of parameterized boundary points
    scalar maxDisplacement(Zero);
    for (const label patchI : patchesToBeMoved)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        const labelList& meshPoints = patch.meshPoints();

        for (const label globalIndex : meshPoints)
        {
            const label whichPointInBox = reverseMapPtr_()[globalIndex];
            // If point resides within control points box,
            // compute new cartesian coordinates
            if (whichPointInBox != -1)
            {
                vector newPoint =
                    transformPointToCartesian
                    (
                       coordinates
                       (
                           paramCoors[globalIndex]
                       )
                    );
                maxDisplacement =
                    max
                    (
                        maxDisplacement,
                        mag(newPoint - mesh_.points()[globalIndex])
                    );
            }
        }
    }
    reduce(maxDisplacement, maxOp<scalar>());
    cps_ = oldCPs;

    return maxDisplacement;
}


Foam::tmp<Foam::vectorField> Foam::NURBS3DVolume::getPointsInBox()
{
    if (!mapPtr_)
    {
        findPointsInBox(localSystemCoordinates_);
    }
    tmp<vectorField> pointsInBox
    (
        new vectorField(localSystemCoordinates_, mapPtr_())
    );

    return pointsInBox;
}


const Foam::labelList& Foam::NURBS3DVolume::getMap()
{
    if (!mapPtr_)
    {
        findPointsInBox(localSystemCoordinates_);
    }

    return mapPtr_();
}


const Foam::labelList& Foam::NURBS3DVolume::getReverseMap()
{
    if (!reverseMapPtr_)
    {
        findPointsInBox(localSystemCoordinates_);
    }

    return reverseMapPtr_();
}


const Foam::pointVectorField& Foam::NURBS3DVolume::getParametricCoordinates()
{
    // If not computed yet, compute parametric coordinates
    if (!parametricCoordinatesPtr_)
    {
        // Find mesh points inside control points box
        // if they have been identified yet
        if (!mapPtr_)
        {
            findPointsInBox(localSystemCoordinates_);
        }
        computeParametricCoordinates(getPointsInBox()());
    }

    return parametricCoordinatesPtr_();
}


Foam::tmp<Foam::pointTensorField> Foam::NURBS3DVolume::getDxDb(const label cpI)
{
    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Set return field to zero
    tmp<pointTensorField> tDxDb
    (
        new pointTensorField
        (
            IOobject
            (
                "DxDb",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pointMesh::New(mesh_),
            dimensionedTensor(dimless, Zero)
        )
    );

    pointTensorField& DxDb = tDxDb.ref();

    // All points outside the control box remain unmoved.
    // Loop over only points within the control box
    const labelList& map = mapPtr_();
    for (const label globalIndex : map)
    {
        DxDb[globalIndex] =
            transformationTensorDxDb(globalIndex)
           *volumeDerivativeCP
            (
                parametricCoordinates[globalIndex],
                cpI
            );
    }

    return tDxDb;
}


Foam::tmp<Foam::volTensorField> Foam::NURBS3DVolume::getDxCellsDb
(
    const label cpI
)
{
    // Get parametric coordinates
    const vectorField& parametricCoordinates = getParametricCoordinates();

    // Set return field to zero
    auto tDxDb = volTensorField::New
    (
        "DxDb",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedTensor(dimless, Zero)
    );

    volTensorField& DxDb = tDxDb.ref();
    deltaBoundary deltaBound(mesh_);
    const labelListList& pointCells = mesh_.pointCells();

    // All points outside the control box remain unmoved.
    // Loop over only points within the control box
    const labelList& map = mapPtr_();
    for (const label globalIndex : map)
    {
        tensor pointDxDb =
            transformationTensorDxDb(globalIndex)
           *volumeDerivativeCP
            (
                parametricCoordinates[globalIndex],
                cpI
            );
        const labelList& pointCellsI = pointCells[globalIndex];
        tmp<tensorField> tC_d = deltaBound.cellCenters_d(globalIndex);
        const tensorField& C_d = tC_d();
        forAll(pointCellsI, cI)
        {
            const label cellI = pointCellsI[cI];
            DxDb[cellI] += C_d[cI] & pointDxDb;
        }
    }

    // Assign boundary values since the grad of this field is often needed
    forAll(mesh_.boundary(), pI)
    {
        const fvPatch& patch = mesh_.boundary()[pI];
        if (!isA<coupledFvPatch>(patch))
        {
            DxDb.boundaryFieldRef()[pI] = patchDxDbFace(pI, cpI);
        }
    }

    // Correct coupled boundaries
    DxDb.correctBoundaryConditions();

    return tDxDb;
}


Foam::label Foam::NURBS3DVolume::nUSymmetry() const
{
    label nU(basisU_.nCPs());
    return label(nU % 2 == 0 ? 0.5*nU : 0.5*(nU - 1) + 1);
}


Foam::label Foam::NURBS3DVolume::nVSymmetry() const
{
    label nV(basisV_.nCPs());
    return label(nV % 2 == 0 ? 0.5*nV : 0.5*(nV - 1) + 1);
}


Foam::label Foam::NURBS3DVolume::nWSymmetry() const
{
    label nW(basisW_.nCPs());
    return label(nW % 2 == 0 ? 0.5*nW : 0.5*(nW - 1) + 1);
}


Foam::Vector<Foam::label> Foam::NURBS3DVolume::nSymmetry() const
{
    return Vector<label>(nUSymmetry(), nVSymmetry(), nWSymmetry());
}


void Foam::NURBS3DVolume::writeCps
(
    const fileName& baseName,
    const bool transform
) const
{
    const label nCPsU = basisU_.nCPs();
    const label nCPsV = basisV_.nCPs();

    vectorField cpsInCartesian(cps_);
    if (transform)
    {
        forAll(cpsInCartesian, cpI)
        {
            cpsInCartesian[cpI] = transformPointToCartesian(cps_[cpI]);
        }
    }

    Info<< "Writing control point positions to file" << endl;

    if (Pstream::master())
    {
        OFstream cpsFile("optimisation"/cpsFolder_/name_ + baseName + ".csv");
        // Write header
        cpsFile
            << "\"Points : 0\", \"Points : 1\", \"Points : 2\","
            << "\"i\", \"j\", \"k\","
            << "\"active : 0\", \"active : 1\", \"active : 2\"" << endl;

        forAll(cpsInCartesian, cpI)
        {
            const label iCPw = cpI/label(nCPsU*nCPsV);
            const label iCPv = (cpI - iCPw*nCPsU*nCPsV)/nCPsU;
            const label iCPu = (cpI - iCPw*nCPsU*nCPsV - iCPv*nCPsU);

            cpsFile
                << cpsInCartesian[cpI].x() << ", "
                << cpsInCartesian[cpI].y() << ", "
                << cpsInCartesian[cpI].z() << ", "
                << iCPu << ", "
                << iCPv << ", "
                << iCPw << ", "
                << activeDesignVariables_[3*cpI] << ", "
                << activeDesignVariables_[3*cpI + 1] << ", "
                << activeDesignVariables_[3*cpI + 2] << endl;
        }
    }
}


void Foam::NURBS3DVolume::writeParamCoordinates() const
{
    parametricCoordinatesPtr_().write();
}


bool Foam::NURBS3DVolume::writeData(Ostream& os) const
{
    cps_.writeEntry("controlPoints", os);
    return true;
}


// ************************************************************************* //
