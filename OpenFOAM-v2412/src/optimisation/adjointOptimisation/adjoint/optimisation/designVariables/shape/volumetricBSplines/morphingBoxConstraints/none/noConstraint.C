/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

#include "noConstraint.H"
#include "volumetricBSplinesDesignVariables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noConstraint, 1);
    addToRunTimeSelectionTable
    (
        morphingBoxConstraint,
        noConstraint,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::noConstraint::computeDVsSensitivities
(
    scalarField& dvSens,
    const scalarField& cpSens
)
{
    dvSens = cpSens;
}


void Foam::noConstraint::updateInternalBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds,
    const NURBS3DVolume& boxI,
    const label passed
)
{
    const vectorField& cps = boxI.getControlPoints();
    const Vector<label> nCPsDir = boxI.nCPsPerDirection();
    // Internal points
    for (label k = 1; k < nCPsDir[2] - 1; ++k)
    {
        for (label j = 1; j < nCPsDir[1] - 1; ++j)
        {
            for (label i = 1; i < nCPsDir[0] - 1; ++i)
            {
                label cpID(boxI.getCPID(i, j, k));
                for (label idir = 0; idir < 3; ++idir)
                {
                    label iIncr(idir == 0);
                    label jIncr(idir == 1);
                    label kIncr(idir == 2);
                    label prevCP
                        (boxI.getCPID(i - iIncr, j - jIncr, k - kIncr));
                    label nextCP
                        (boxI.getCPID(i + iIncr, j + jIncr, k + kIncr));
                    lowerBounds()[3*cpID + idir + passed] =
                        0.5
                       *(
                            cps[prevCP].component(idir)
                          + cps[cpID].component(idir)
                        );
                    upperBounds()[3*cpID + idir + passed] =
                        0.5
                       *(
                            cps[nextCP].component(idir)
                          + cps[cpID].component(idir)
                        );
                }
            }
        }
    }
}


void Foam::noConstraint::updateBoundaryBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds,
    const NURBS3DVolume& boxI,
    const label passed
)
{
    const vectorField& cps = boxI.getControlPoints();
    const Vector<label> nCPsDir = boxI.nCPsPerDirection();
    // Loop over boundaries in all directions of the box
    for (label ibound = 0; ibound < 3; ++ibound)
    {
        // Start of iterators in the three directions
        Vector<label> minID(1, 1, 1);
        // End of iterators in the three directions
        Vector<label> maxID(nCPsDir[0] - 2, nCPsDir[1] - 2, nCPsDir[2] - 2);
        // Increment of iterators in the three directions
        Vector<label> incr(1, 1, 1);

        // Adjust looping in the direction we are examining
        minID[ibound] = 0;
        maxID[ibound] = nCPsDir[ibound];
        incr[ibound] = nCPsDir[ibound] - 1;
        Vector<label> indices(Zero);
        label& i = indices[0];
        label& j = indices[1];
        label& k = indices[2];

        for (k = minID[2]; k < maxID[2]; k += incr[2])
        {
            for (j = minID[1]; j < maxID[1]; j += incr[1])
            {
                for (i = minID[0]; i < maxID[0]; i += incr[0])
                {
                    label cpID(boxI.getCPID(i, j, k));
                    for (label dir = 0; dir < 3; ++dir)
                    {
                        Vector<label> incrMinus(dir == 0, dir == 1, dir == 2);
                        Vector<label> incrPlus(dir == 0, dir == 1, dir == 2);
                        // Adjust increment for the ibound direction
                        incrMinus[ibound] =
                            label
                            (
                                incrMinus[ibound]
                             && indices[ibound] == nCPsDir[ibound] - 1
                            );
                        incrPlus[ibound] =
                            label(incrMinus[ibound] && indices[ibound] == 0);
                        label prevCP =
                            boxI.getCPID
                            (
                                i - incrMinus[0],
                                j - incrMinus[1],
                                k - incrMinus[2]
                            );
                        label nextCP =
                            boxI.getCPID
                            (
                                i + incrPlus[0],
                                j + incrPlus[1],
                                k + incrPlus[2]
                            );
                        if (incrMinus[ibound])
                        {
                            lowerBounds()[3*cpID + dir + passed] =
                                0.5
                               *(
                                    cps[prevCP].component(dir)
                                  + cps[cpID].component(dir)
                                );
                        }
                        if (incrPlus[ibound])
                        {
                            upperBounds()[3*cpID + dir + passed] =
                                0.5
                               *(
                                    cps[nextCP].component(dir)
                                  + cps[cpID].component(dir)
                                );
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noConstraint::noConstraint
(
    const fvMesh& mesh,
    const dictionary& dict,
    volumetricBSplinesDesignVariables& designVariables
)
:
    morphingBoxConstraint(mesh, dict, designVariables)
{
    designVariables_.setSize(3*volBSplinesBase_.getTotalControlPointsNumber());
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::noConstraint::computeBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds
)
{
    // Does nothing
}


void Foam::noConstraint::updateBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds
)
{
    if (designVariables_.nonOverlappingCPs() && designVariables_.updateBounds())
    {
        DebugInfo
            << "Updating bounds for the design variables " << endl;
        const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        label passed(0);
        for (const NURBS3DVolume& boxI : boxes)
        {
            // Bounds for internal control points
            updateInternalBounds(lowerBounds, upperBounds, boxI, passed);
            // Bounds for boundary points.
            // Assumes that the boundary edges remain fixed
            updateBoundaryBounds(lowerBounds, upperBounds, boxI, passed);
            passed += 3*boxI.getControlPoints().size();
        }
        DebugInfo
            << "lower bounds " << lowerBounds() << endl;
        DebugInfo
            << "upper bounds " << upperBounds() << endl;
    }
}


Foam::labelList Foam::noConstraint::computeActiveDesignVariables
(
    const labelList& activeCPCoors
)
{
    return activeCPCoors;
}


Foam::tmp<Foam::scalarField> Foam::noConstraint::designVariablesToControlPoints
(
    const scalarField& designVariables
)
{
    return designVariables;
}


Foam::tmp<Foam::scalarField> Foam::noConstraint::controlPointsToDesignVariables
(
    const scalarField& cps
)
{
    return cps;
}


Foam::tmp<Foam::scalarField> Foam::noConstraint::correctionCPs
(
    const scalarField& correction
)
{
    return correction;
}


// ************************************************************************* //
