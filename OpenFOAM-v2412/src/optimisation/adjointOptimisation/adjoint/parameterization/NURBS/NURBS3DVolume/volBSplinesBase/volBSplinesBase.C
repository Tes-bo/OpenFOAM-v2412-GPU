/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "volBSplinesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(volBSplinesBase, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volBSplinesBase::volBSplinesBase
(
    const fvMesh& mesh
)
:
    MeshObject_type(mesh),
    volume_(0),
    activeDesignVariables_(0)
{
    const dictionary NURBSdict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).subDict("volumetricBSplinesMotionSolverCoeffs")
    );

    // Populate NURBS volumes
    volume_.resize(NURBSdict.size());

    label iBox(0);

    for (const entry& dEntry : NURBSdict)
    {
        if (dEntry.isDict())
        {
            volume_.set
            (
                iBox,
                NURBS3DVolume::New(dEntry.dict(), mesh, true)
            );
            volume_[iBox].writeParamCoordinates();
            iBox++;
        }
    }
    volume_.resize(iBox);

    // Determine active design variables
    activeDesignVariables_.setSize(3*getTotalControlPointsNumber(), -1);
    label iActive(0);
    const labelList startCpID(getStartCpID());
    forAll(volume_, boxI)
    {
        const label start(3*startCpID[boxI]);
        const boolList& isActiveVar = volume_[boxI].getActiveDesignVariables();
        forAll(isActiveVar, varI)
        {
            if (isActiveVar[varI])
            {
                activeDesignVariables_[iActive++] = start + varI;
            }
        }
    }
    activeDesignVariables_.setSize(iActive);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const PtrList<NURBS3DVolume>& volBSplinesBase::boxes() const
{
    return volume_;
}


PtrList<NURBS3DVolume>& volBSplinesBase::boxesRef()
{
    return volume_;
}


const NURBS3DVolume& volBSplinesBase::box(const label boxI) const
{
    return volume_[boxI];
}


NURBS3DVolume& volBSplinesBase::boxRef(const label boxI)
{
    return volume_[boxI];
}


const vectorField& volBSplinesBase::getControlPoints(const label& iNURB) const
{
    return volume_[iNURB].getControlPoints();
}


vectorField volBSplinesBase::getAllControlPoints() const
{
    DynamicList<vector> totalCPs(0);
    forAll(volume_, iNURB)
    {
        totalCPs.push_back(volume_[iNURB].getControlPoints());
    }

    return vectorField(std::move(totalCPs));
}


Foam::label Foam::volBSplinesBase::getTotalControlPointsNumber() const
{
    label nCPs(0);
    forAll(volume_, iNURB)
    {
        nCPs += volume_[iNURB].getControlPoints().size();
    }

    return nCPs;
}


label Foam::volBSplinesBase::getNumberOfBoxes() const
{
    return volume_.size();
}


labelList volBSplinesBase::getStartCpID() const
{
    // Allocate an extra entry to track in which box a CP might be
    labelList startID(getNumberOfBoxes() + 1);
    startID[0] = 0;
    forAll(volume_, iNURB)
    {
        startID[iNURB+1] =
            startID[iNURB] + volume_[iNURB].getControlPoints().size();
    }

    return startID;
}


labelList volBSplinesBase::getStartVarID() const
{
    return 3*getStartCpID();
}


label volBSplinesBase::findBoxID(const label cpI) const
{
    const labelList startCPID(getStartCpID());
    for (label iBox = 0; iBox < startCPID.size() - 1 ; ++iBox)
    {
        if (cpI >= startCPID[iBox] || cpI < startCPID[iBox + 1])
        {
            return iBox;
        }
    }

    FatalErrorInFunction
        << "Invalid control point ID " << cpI << endl
        << exit(FatalError);
    return -1;
}


Vector<label> volBSplinesBase::decomposeDV(const label varID) const
{
    Vector<label> decomposed;
    labelList startVarID = getStartVarID();
    label boxID(-1);
    for (label iBox = 0; iBox < startVarID.size() - 1 ; ++iBox)
    {
        if (varID >= startVarID[iBox] && varID < startVarID[iBox + 1])
        {
            boxID = iBox;
            break;
        }
    }
    const label localVarID = varID - startVarID[boxID];
    decomposed.x() = boxID;
    decomposed.y() = localVarID/3;
    decomposed.z() = localVarID%3;
    DebugInfo
        << "varID " << varID
        << " belongs to box " << decomposed.x()
        << " cpLocal " << decomposed.y()
        << " dir " << decomposed.z()
        << endl;
    return decomposed;
}


const Foam::labelList& volBSplinesBase::getActiveDesignVariables() const
{
    return activeDesignVariables_;
}


Foam::scalar Foam::volBSplinesBase::computeMaxBoundaryDisplacement
(
    const vectorField& controlPointsMovement,
    const labelList& patchesToBeMoved
)
{
    scalar maxDisplacement(0);
    label pastControlPoints(0);
    forAll(volume_, iNURB)
    {
        const label nb(volume_[iNURB].getControlPoints().size());
        vectorField localControlPointsMovement(nb, Zero);

        // Set localControlPointsMovement
        forAll(localControlPointsMovement, iCPM)
        {
            localControlPointsMovement[iCPM] =
                controlPointsMovement[pastControlPoints + iCPM];
        }

        maxDisplacement = max
        (
            maxDisplacement,
            volume_[iNURB].computeMaxBoundaryDisplacement
            (
                localControlPointsMovement,
                patchesToBeMoved
            )
        );

        pastControlPoints += nb;
    }

    return maxDisplacement;
}


Foam::tmp<vectorField> Foam::volBSplinesBase::computeBoundaryDisplacement
(
    const vectorField& controlPointsMovement,
    const labelList& patchesToBeMoved
)
{
    auto tdisplacement(tmp<vectorField>::New(mesh_.nPoints(), Zero));
    vectorField& displacement = tdisplacement.ref();

    label pastControlPoints(0);
    forAll(volume_, iNURB)
    {
        const label nb(volume_[iNURB].getControlPoints().size());
        vectorField localControlPointsMovement(nb, Zero);

        // Set localControlPointsMovement
        forAll(localControlPointsMovement, iCPM)
        {
            localControlPointsMovement[iCPM] =
                controlPointsMovement[pastControlPoints + iCPM];
        }

        displacement +=
            volume_[iNURB].computeNewBoundaryPoints
            (
                localControlPointsMovement,
                patchesToBeMoved
            )
          - mesh_.points();

        pastControlPoints += nb;
    }

    return tdisplacement;
}


void Foam::volBSplinesBase::boundControlPointMovement
(
    vectorField& controlPointsMovement
) const
{
    label pastControlPoints(0);
    forAll(volume_, iNURB)
    {
        const label nb(volume_[iNURB].getControlPoints().size());
        vectorField localControlPointsMovement(nb, Zero);

        // Set localControlPointsMovement
        forAll(localControlPointsMovement, iCPM)
        {
            localControlPointsMovement[iCPM] =
                controlPointsMovement[pastControlPoints + iCPM];
        }

        volume_[iNURB].boundControlPointMovement(localControlPointsMovement);

        // Transfer bounding back to controlPointMovement
        forAll(localControlPointsMovement, iCPM)
        {
            controlPointsMovement[pastControlPoints + iCPM] =
                localControlPointsMovement[iCPM];
        }

        pastControlPoints += nb;
    }
}


void Foam::volBSplinesBase::moveControlPoints
(
    const vectorField& controlPointsMovement
)
{
    label pastControlPoints(0);
    forAll(volume_, iNURB)
    {
        const label nb(volume_[iNURB].getControlPoints().size());
        vectorField localControlPointsMovement(nb, Zero);

        // Set localControlPointsMovement
        forAll(localControlPointsMovement, iCPM)
        {
            localControlPointsMovement[iCPM] =
                controlPointsMovement[pastControlPoints + iCPM];
        }

        const vectorField newCps
        (
            volume_[iNURB].getControlPoints()
          + localControlPointsMovement
        );

        volume_[iNURB].setControlPoints(newCps);

        pastControlPoints += nb;
    }
}


void Foam::volBSplinesBase::writeControlPoints() const
{
    for (const NURBS3DVolume& box : volume_)
    {
        box.writeCps("cpsBsplines" + mesh_.time().timeName());
    }
}


bool volBSplinesBase::movePoints()
{
    // Does nothing
    return true;
}


void volBSplinesBase::updateMesh(const mapPolyMesh&)
{
    // Does nothing
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
