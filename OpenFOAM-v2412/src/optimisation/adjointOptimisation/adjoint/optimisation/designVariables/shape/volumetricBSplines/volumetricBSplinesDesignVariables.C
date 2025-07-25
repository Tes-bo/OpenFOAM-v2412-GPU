/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "volumetricBSplinesDesignVariables.H"
#include "pointVolInterpolation.H"
#include "displacementMethodvolumetricBSplinesMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volumetricBSplinesDesignVariables, 0);
    addToRunTimeSelectionTable
    (
        shapeDesignVariables,
        volumetricBSplinesDesignVariables,
        dictionary
    );
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::label Foam::volumetricBSplinesDesignVariables::sensSize() const
{
    return 3*volBSplinesBase_.getTotalControlPointsNumber();
}


const Foam::labelList&
Foam::volumetricBSplinesDesignVariables::activeSensitivities() const
{
    return volBSplinesBase_.getActiveDesignVariables();
}


void Foam::volumetricBSplinesDesignVariables::setActiveDesignVariables()
{
    // Active design variables pertaining to the CPs numbering
    labelList activeVarsInCPs = volBSplinesBase_.getActiveDesignVariables();

    // Convert the aforementioned list to the numbering of the actual design
    // variables
    activeDesignVariables_ =
        constraint_().computeActiveDesignVariables(activeVarsInCPs);
}


void Foam::volumetricBSplinesDesignVariables::designVariablesToControlPoints()
{
    const scalarField& dvs = *this;

    // Convert design variables to CPs coordinates, stored in a scalarField
    scalarField cpsScalar(constraint_().designVariablesToControlPoints(dvs));

    // Convert the scalarField to vectorFields and transfer to morphing boxes
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    label varID(0);
    for (NURBS3DVolume& boxI : boxes)
    {
        vectorField cps(boxI.getControlPoints().size(), Zero);
        for (vector& cpI : cps)
        {
            cpI.x() = cpsScalar[varID++];
            cpI.y() = cpsScalar[varID++];
            cpI.z() = cpsScalar[varID++];
        }
        boxI.setControlPoints(cps);
    }
}


void Foam::volumetricBSplinesDesignVariables::controlPointsToDesignVariables()
{
    // Store CP coordinates to a scalarField
    scalarField cpsScalar(3*volBSplinesBase_.getTotalControlPointsNumber());
    const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxes();
    label varID(0);
    for (const NURBS3DVolume& boxI : boxes)
    {
        const vectorField& cps = boxI.getControlPoints();
        for (const vector& cpI : cps)
        {
            cpsScalar[varID++] = cpI.x();
            cpsScalar[varID++] = cpI.y();
            cpsScalar[varID++] = cpI.z();
        }
    }

    // Convert this scalarField to the design variables
    scalarField::operator=
        (constraint_().controlPointsToDesignVariables(cpsScalar));
}


void Foam::volumetricBSplinesDesignVariables::controlPointsToDesignVariables
(
    const vectorField& controlPoints
)
{
    // Store CP coordinates to a scalarField
    scalarField cpsScalar(3*volBSplinesBase_.getTotalControlPointsNumber());
    const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxes();
    label varID(0);
    for (const NURBS3DVolume& boxI : boxes)
    {
        const label nCPs(boxI.getControlPoints().size());
        for (label cpI = 0; cpI < nCPs; ++cpI)
        {
            cpsScalar[varID++] = controlPoints[cpI].x();
            cpsScalar[varID++] = controlPoints[cpI].y();
            cpsScalar[varID++] = controlPoints[cpI].z();
        }
    }

    // Convert this scalarField to the design variables
    scalarField::operator=
        (constraint_().controlPointsToDesignVariables(cpsScalar));
}


void Foam::volumetricBSplinesDesignVariables::readBounds
(
    autoPtr<scalar> lowerBoundPtr,
    autoPtr<scalar> upperBoundPtr
)
{
    designVariables::readBounds(lowerBoundPtr, upperBoundPtr);
    readBounds(lowerBounds_, "lower", -1);
    readBounds(upperBounds_, "upper",  1);

    // Update bounds based on the constraints - WIP
    constraint_().computeBounds(lowerBounds_, upperBounds_);
}


void Foam::volumetricBSplinesDesignVariables::readBounds
(
    autoPtr<scalarField>& bounds,
    const word& boundsName,
    const label sign
)
{
    // Read global bounds for the control points
    if (dict_.found(boundsName + "CPBounds"))
    {
        bounds.reset(new scalarField(getVars().size()));

        vector CPBounds(dict_.get<vector>(boundsName + "CPBounds"));
        const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        label varID(0);
        for (const NURBS3DVolume& boxI : boxes)
        {
            const label nCPs(boxI.getControlPoints().size());
            for (label iCP = 0; iCP < nCPs; ++iCP)
            {
                bounds()[varID++] = CPBounds.x();
                bounds()[varID++] = CPBounds.y();
                bounds()[varID++] = CPBounds.z();
            }
        }
    }
    // Read in bounds from the designVariables dictionary if present.
    // If nonOverlappingCPs is used, the current CPs are used to determine the
    // bounds of the CPs. If we continue from a previous solution, the current
    // CPs are different from the initial ones and, hence, different bounds
    // will be computed for the continuation run. Instead, read the bounds
    // from the designVariables dict, if present
    else if (localIOdictionary::found(boundsName + "Bounds"))
    {
        DebugInfo
            << "Reading " << boundsName << "Bounds from dict " << endl;
        bounds.reset
            (new scalarField(boundsName + "Bounds", *this, getVars().size()));

    }
    else if (nonOverlappingCPs_)
    {
        DebugInfo
            << "Setting " << boundsName << "Bounds from nonOverlappingCPs"
            << endl;
        bounds.reset(new scalarField(getVars().size()));
        const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        label varID(0);
        for (const NURBS3DVolume& boxI : boxes)
        {
            const vectorField& cps = boxI.getControlPoints();
            const Vector<label> nCPsDir = boxI.nCPsPerDirection();
            vector dists(Zero);
            for (label idir = 0; idir < 3; ++idir)
            {
                dists[idir] =
                    (max(cps.component(idir)) - min(cps.component(idir)))
                   /scalar(nCPsDir[idir] - 1);
            }
            const label nCPs(boxI.getControlPoints().size());
            for (label iCP = 0; iCP < nCPs; ++iCP)
            {
                const vector& cp = cps[iCP];
                bounds()[varID++] = cp.x() + sign*0.5*dists.x();
                bounds()[varID++] = cp.y() + sign*0.5*dists.y();
                bounds()[varID++] = cp.z() + sign*0.5*dists.z();
            }
        }
    }
}


void Foam::volumetricBSplinesDesignVariables::writeBounds
(
    const scalarField& bounds,
    const word& name
) const
{
    if (Pstream::master())
    {
        const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        label passed(0);
        for (const NURBS3DVolume& boxI : boxes)
        {
            OFstream file
            (
                word("optimisation")/word("controlPoints")/boxI.name()
              + name + mesh_.time().timeName() + ".csv"
            );
            // Write header
            file<< "\"Points : 0\", \"Points : 1\", \"Points : 2\","
                << "\"i\", \"j\", \"k\""<< endl;

            const vectorField& cps = boxI.getControlPoints();
            const label nCPsU = boxI.basisU().nCPs();
            const label nCPsV = boxI.basisV().nCPs();
            forAll(cps, cpI)
            {
                const label k = cpI/label(nCPsU*nCPsV);
                const label j = (cpI - k*nCPsU*nCPsV)/nCPsU;
                const label i = (cpI - k*nCPsU*nCPsV - j*nCPsU);

                file<< bounds[3*cpI + passed] << ", "
                    << bounds[3*cpI + 1 + passed] << ", "
                    << bounds[3*cpI + 2 + passed] << ", "
                    << i << ", "
                    << j << ", "
                    << k << endl;
            }
            passed += 3*cps.size();
        }
    }
}


void Foam::volumetricBSplinesDesignVariables::setDisplacement
(
    const vectorField& cpMovement
)
{
    displacementMethod& dm = displMethodPtr_.ref();
    // Are volumetric B-Splines also used to move the mesh ?
    if (isA<displacementMethodvolumetricBSplinesMotionSolver>(dm))
    {
        // Communicate the control points movement to the displacement method
        displMethodPtr_->setControlField(cpMovement);
    }
    else
    {
        // This will also update the control point positions
        tmp<vectorField> tnewPoints =
            volBSplinesBase_.computeBoundaryDisplacement
            (
                cpMovement,
                parametertisedPatches_.toc()
            );
        const vectorField& newPoints = tnewPoints();

        pointVectorField dx
        (
            IOobject
            (
                "dx",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            pointMesh::New(mesh_),
            dimensionedVector(dimless, Zero)
        );

        for (const label pI : parametertisedPatches_)
        {
            dx.boundaryField()[pI].setInInternalField
            (
                dx.primitiveFieldRef(),
                vectorField(newPoints, mesh_.boundaryMesh()[pI].meshPoints())
            );
        }

        // Set boundary movement of motion solver
        displMethodPtr_->setMotionField(dx);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volumetricBSplinesDesignVariables::volumetricBSplinesDesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    shapeDesignVariables(mesh, dict),
    localIOdictionary
    (
        IOobject
        (
            "volumetricBSplinesDesignVariables",
            mesh.time().timeName(),
            fileName("uniform"),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    volBSplinesBase_(const_cast<volBSplinesBase&>(volBSplinesBase::New(mesh))),
    nonOverlappingCPs_(dict_.getOrDefault<bool>("nonOverlappingCPs", false)),
    updateBounds_(dict_.getOrDefault<bool>("updateBounds", true)),
    constraint_(morphingBoxConstraint::New(mesh, dict, *this))
{
    // Read in design variables if present or initialise them
    if (localIOdictionary::found("designVariables"))
    {
        scalarField::operator=
            (scalarField("designVariables", *this, scalarField::size()));
    }
    else if (constraint_().initialiseVars())
    {
        controlPointsToDesignVariables();
    }

    // Set the active design variables
    setActiveDesignVariables();

    // Read bounds for design variables, if present
    readBounds();
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::volumetricBSplinesDesignVariables::controlPointMovement
(
    const scalarField& correction
)
{
    auto tcpMovement
    (
        tmp<vectorField>::New
        (
            volBSplinesBase_.getTotalControlPointsNumber(),
            Zero
        )
    );
    vectorField& cpMovement = tcpMovement.ref();

    // Convert the correction pertaining to the design variables to a
    // scalarField correction for the control points
    const scalarField correctionCPs(constraint_().correctionCPs(correction));

    // scalarField to vectorField conversion
    forAll(cpMovement, iCP)
    {
        cpMovement[iCP].x() = correctionCPs[3*iCP];
        cpMovement[iCP].y() = correctionCPs[3*iCP + 1];
        cpMovement[iCP].z() = correctionCPs[3*iCP + 2];
    }
    volBSplinesBase_.boundControlPointMovement(cpMovement);

    return tcpMovement;
}


void Foam::volumetricBSplinesDesignVariables::update(scalarField& correction)
{
    // Get controlPoint movement from correction
    tmp<vectorField> tcpMovement = controlPointMovement(correction);
    const vectorField& cpMovement = tcpMovement();

    // Set the field driving the displacement method
    setDisplacement(cpMovement);

    // Do the actual mesh movement
    // Updates also the control point positions
    moveMesh();

    // Update the design variables
    scalarField::operator+=(correction);
}


void Foam::volumetricBSplinesDesignVariables::resetDesignVariables()
{
    shapeDesignVariables::resetDesignVariables();
    designVariablesToControlPoints();
}


Foam::scalar
Foam::volumetricBSplinesDesignVariables::computeEta(scalarField& correction)
{
    return constraint_().computeEta(correction, maxInitChange_());
}


bool Foam::volumetricBSplinesDesignVariables::globalSum() const
{
    return false;
}


Foam::tmp<Foam::scalarField>
Foam::volumetricBSplinesDesignVariables::assembleSensitivities
(
    adjointSensitivity& adjointSens
)
{
    return
        constraint_().postProcessSens
        (
            shapeDesignVariables::assembleSensitivities(adjointSens)(),
            adjointSens.getAdjointSolver().solverName()
        );
}


void Foam::volumetricBSplinesDesignVariables::evolveNumber()
{
    constraint_().updateBounds(lowerBounds_, upperBounds_);
}


bool Foam::volumetricBSplinesDesignVariables::writeData(Ostream& os) const
{
    scalarField::writeEntry("designVariables", os);
    if (lowerBounds_)
    {
        lowerBounds_().writeEntry("lowerBounds", os);
        writeBounds(lowerBounds_(), "lowerBounds");
    }
    if (upperBounds_)
    {
        upperBounds_().writeEntry("upperBounds", os);
        writeBounds(upperBounds_(), "upperBounds");
    }
    return constraint_().writeData(os);
}


Foam::tmp<Foam::vectorField> Foam::volumetricBSplinesDesignVariables::dxdbVol
(
    const label varID
) const
{
    const displacementMethod& dm = displMethodPtr_();
    if (isA<displacementMethodvolumetricBSplinesMotionSolver>(dm))
    {
        Vector<label> decomposed = volBSplinesBase_.decomposeDV(varID);
        const label boxI = decomposed.x();
        const label cpILocal = decomposed.y();
        const label dir = decomposed.z();

        pointTensorField dxdb(volBSplinesBase_.boxRef(boxI).getDxDb(cpILocal));
        return unzipCol(dxdb, dir);
    }
    return tmp<vectorField>::New(0);
}


Foam::tmp<Foam::vectorField> Foam::volumetricBSplinesDesignVariables::dxdbFace
(
    const label patchI,
    const label varID
) const
{
    Vector<label> decomposed = volBSplinesBase_.decomposeDV(varID);
    const label boxI = decomposed.x();
    const label cpILocal = decomposed.y();
    const label dir = decomposed.z();

    tensorField dxdb
        (volBSplinesBase_.boxRef(boxI).patchDxDbFace(patchI, cpILocal));
    return unzipCol(dxdb, dir);
}


Foam::tmp<Foam::vectorField> Foam::volumetricBSplinesDesignVariables::dndb
(
    const label patchI,
    const label varID
) const
{
    Vector<label> decomposed = volBSplinesBase_.decomposeDV(varID);
    const label boxI = decomposed.x();
    const label cpILocal = decomposed.y();
    const label dir = decomposed.z();

    tensorField dndb
    (
        volBSplinesBase_.boxRef(boxI).
            dndbBasedSensitivities(patchI, cpILocal, false)
    );
    return unzipCol(dndb, dir);
}


Foam::tmp<Foam::vectorField> Foam::volumetricBSplinesDesignVariables::dSdb
(
    const label patchI,
    const label varID
) const
{
    Vector<label> decomposed = volBSplinesBase_.decomposeDV(varID);
    const label boxI = decomposed.x();
    const label cpILocal = decomposed.y();
    const label dir = decomposed.z();

    tensorField dndb
    (
        volBSplinesBase_.boxRef(boxI).dndbBasedSensitivities(patchI, cpILocal)
    );
    return unzipCol(dndb, dir);
}


Foam::tmp<Foam::volVectorField> Foam::volumetricBSplinesDesignVariables::dCdb
(
    const label varID
) const
{
    Vector<label> decomposed = volBSplinesBase_.decomposeDV(varID);
    const label boxI = decomposed.x();
    const label cpILocal = decomposed.y();
    const label dir = decomposed.z();
    NURBS3DVolume& box = volBSplinesBase_.boxRef(boxI);
    pointVolInterpolation volPointInter(pointMesh::New(mesh_), mesh_);
    // WIP: we compute the entire dxdb tensor corresponding to the contol point
    // and then extract the desired direction. This is quite expensive.
    // Specific functions returning what we want should be implemented in
    // NURBS3DVolume to reduce the cost
    tmp<volTensorField> dxdb = volPointInter.interpolate(box.getDxDb(cpILocal));
    auto tdxdbDir =
        tmp<volVectorField>::New
        (
            IOobject
            (
                "dxdbDir",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimless, Zero)
        );
    volVectorField& dxdbDir = tdxdbDir.ref();
    unzipCol(dxdb(), vector::components(dir), dxdbDir);
    return tdxdbDir;
}


// ************************************************************************* //
