/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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
#include "transformBox.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(transformBox, 0);
    addToRunTimeSelectionTable
    (
        controlPointsDefinition,
        transformBox,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::transformBox::computeControlPoints()
{
    const label nCPsU(box_.basisU().nCPs());
    const label nCPsV(box_.basisV().nCPs());
    const label nCPsW(box_.basisW().nCPs());
    cps_.setSize(nCPsU*nCPsV*nCPsW, vector::zero);

    // Geometry bounds (the one on which clip was applied)
    const dictionary& dict = box_.dict();
    vector lowerCpBounds(dict.get<vector>("lowerCpBounds"));
    vector upperCpBounds(dict.get<vector>("upperCpBounds"));
    scalar spanU(upperCpBounds.x() - lowerCpBounds.x());
    scalar spanV(upperCpBounds.y() - lowerCpBounds.y());
    scalar spanW(upperCpBounds.z() - lowerCpBounds.z());

    // First, create morphing box based on lowerCpBounds, Max
    for (label iCPw = 0; iCPw < nCPsW; ++iCPw)
    {
        for (label iCPv = 0; iCPv < nCPsV; ++iCPv)
        {
            for (label iCPu = 0; iCPu < nCPsU; ++iCPu)
            {
                const label cpID(box_.getCPID(iCPu, iCPv, iCPw));
                cps_[cpID] = vector
                (
                    lowerCpBounds.x() + scalar(iCPu)/scalar(nCPsU - 1)*spanU,
                    lowerCpBounds.y() + scalar(iCPv)/scalar(nCPsV - 1)*spanV,
                    lowerCpBounds.z() + scalar(iCPw)/scalar(nCPsW - 1)*spanW
                );
            }
        }
    }

    // Transform control points
    transformControlPoints(lowerCpBounds, upperCpBounds);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transformBox::transformBox(NURBS3DVolume& box)
:
    controlPointsDefinition(box)
{
    computeControlPoints();
}


// ************************************************************************* //
