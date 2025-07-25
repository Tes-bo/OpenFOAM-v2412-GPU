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
#include "axisAligned.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(axisAligned, 0);
    addToRunTimeSelectionTable
    (
        controlPointsDefinition,
        axisAligned,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::axisAligned::computeControlPoints()
{
    const label nCPsU(box_.basisU().nCPs());
    const label nCPsV(box_.basisV().nCPs());
    const label nCPsW(box_.basisW().nCPs());

    cps_.setSize(nCPsU*nCPsV*nCPsW);

    vector lowerBounds(box_.dict().get<vector>("lowerCpBounds"));
    vector upperBounds(box_.dict().get<vector>("upperCpBounds"));
    scalar spanU(upperBounds.x() - lowerBounds.x());
    scalar spanV(upperBounds.y() - lowerBounds.y());
    scalar spanW(upperBounds.z() - lowerBounds.z());

    // Equidistribute cps in th u, v, w directions
    for (label iCPw = 0; iCPw < nCPsW; ++iCPw)
    {
        for (label iCPv = 0; iCPv < nCPsV; ++iCPv)
        {
            for (label iCPu = 0; iCPu < nCPsU; ++iCPu)
            {
                const label cpID(box_.getCPID(iCPu, iCPv, iCPw));
                cps_[cpID] = vector
                (
                    lowerBounds.x() + scalar(iCPu)/scalar(nCPsU - 1)*spanU,
                    lowerBounds.y() + scalar(iCPv)/scalar(nCPsV - 1)*spanV,
                    lowerBounds.z() + scalar(iCPw)/scalar(nCPsW - 1)*spanW
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::axisAligned::axisAligned(NURBS3DVolume& box)
:
    controlPointsDefinition(box)
{
    computeControlPoints();
}


// ************************************************************************* //
