/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "porosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityModel, 0);
    defineRunTimeSelectionTable(porosityModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::porosityModel::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = cmptMax(resist.value());

    if (maxCmpt < 0)
    {
        FatalErrorInFunction
            << "Cannot have all resistances set to negative, resistance = "
            << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel::porosityModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const wordRe& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(true),
    zoneName_(cellZoneName),
    cellZoneIDs_(),
    csysPtr_
    (
        coordinateSystem::New(mesh, coeffs_, coordinateSystem::typeName)
    )
{
    if (zoneName_.empty())
    {
        dict.readIfPresent("active", active_);
        dict_.readEntry("cellZone", zoneName_);
    }

    cellZoneIDs_ = mesh_.cellZones().indices(zoneName_);

    Info<< "    creating porous zone: " << zoneName_ << endl;

    if (returnReduceAnd(cellZoneIDs_.empty()) && Pstream::master())
    {
        FatalErrorInFunction
            << "Cannot find porous cellZone " << zoneName_ << endl
            << "Valid zones : "
            << flatOutput(mesh_.cellZones().names()) << nl
            << "Valid groups: "
            << flatOutput(mesh_.cellZones().groupNames()) << nl
            << exit(FatalError);
    }

    Info<< incrIndent << indent << csys() << decrIndent << endl;

    const pointField& points = mesh_.points();
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];

        boundBox bb;

        for (const label celli : cZone)
        {
            const cell& c = cells[celli];
            const pointField cellPoints(c.points(faces, points));

            for (const point& pt : cellPoints)
            {
                bb.add(csys().localPosition(pt));
            }
        }

        bb.reduce();

        Info<< "    local bounds: " << bb.span() << nl << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModel::transformModelData()
{
    if (!mesh_.upToDatePoints(*this))
    {
        calcTransformModelData();

        // set model up-to-date wrt points
        mesh_.setUpToDatePoints(*this);
    }
}


Foam::tmp<Foam::vectorField> Foam::porosityModel::porosityModel::force
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    transformModelData();

    auto tforce = tmp<vectorField>::New(U.size(), Zero);

    if (!cellZoneIDs_.empty())
    {
        this->calcForce(U, rho, mu, tforce.ref());
    }

    return tforce;
}


void Foam::porosityModel::addResistance(fvVectorMatrix& UEqn)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn);
}


void Foam::porosityModel::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}


void Foam::porosityModel::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, AU);

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


bool Foam::porosityModel::writeData(Ostream& os) const
{
    return true;
}


bool Foam::porosityModel::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    coeffs_ = dict.optionalSubDict(type() + "Coeffs");

    dict.readEntry("cellZone", zoneName_);
    cellZoneIDs_ = mesh_.cellZones().indices(zoneName_);

    return true;
}


// ************************************************************************* //
