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
#include "fromFile.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(fromFile, 0);
    addToRunTimeSelectionTable
    (
        controlPointsDefinition,
        fromFile,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fromFile::computeControlPoints()
{
    Info<< "Reading control points from file " << endl;
    const fvMesh& mesh = box_.mesh();
    const dictionary& dict = box_.dict();
    IOdictionary cpsDict
    (
        IOobject
        (
            dict.dictName() + "cpsBsplines" + mesh.time().timeName(),
            mesh.time().caseConstant(),
            "controlPoints",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    cpsDict.readEntry("controlPoints", cps_);
    const label nCPsU(box_.basisU().nCPs());
    const label nCPsV(box_.basisV().nCPs());
    const label nCPsW(box_.basisW().nCPs());
    if (cps_.size() != nCPsU*nCPsV*nCPsW)
    {
        FatalErrorInFunction
           << "Number of control points does not agree with "
           << "nCPsU*nCPv*nCPsW"
           << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fromFile::fromFile(NURBS3DVolume& box)
:
    controlPointsDefinition(box)
{
    computeControlPoints();
}


// ************************************************************************* //
