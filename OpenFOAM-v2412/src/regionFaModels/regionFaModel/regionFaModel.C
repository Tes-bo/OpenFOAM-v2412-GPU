/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "regionFaModel.H"
#include "faMesh.H"
#include "Time.H"
#include "mappedWallPolyPatch.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionFaModel, 0);
}
}

const Foam::word
Foam::regionModels::regionFaModel::regionFaModelName("regionFaModel");

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionFaModel::constructMeshObjects()
{
    regionMeshPtr_.reset
    (
        new faMesh(primaryMesh_)
    );
}


void Foam::regionModels::regionFaModel::initialise()
{
    if (debug)
    {
        Pout<< "regionFaModel::initialise()" << endl;
    }

    vsmPtr_.reset(new volSurfaceMapping(regionMeshPtr_()));

    if (!outputPropertiesPtr_)
    {
        const fileName uniformPath(word("uniform")/regionFaModelName);

        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    regionName_ + "OutputProperties",
                    time_.timeName(),
                    uniformPath/regionName_,
                    primaryMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


bool Foam::regionModels::regionFaModel::init(const dictionary& dict)
{
    if (active_)
    {
        if (const dictionary* dictptr = dict.findDict(modelName_ + "Coeffs"))
        {
            coeffs_ <<= *dictptr;
        }

        infoOutput_.readIfPresent("infoOutput", dict);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

const Foam::volSurfaceMapping& Foam::regionModels::regionFaModel::vsm() const
{
    return vsmPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionFaModel::regionFaModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(regionFaModelName, dict.get<word>("region")),
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(dict.get<Switch>("active")),
    infoOutput_(false),
    modelName_(modelName),
    regionMeshPtr_(nullptr),
    coeffs_(dict.subOrEmptyDict(modelName + "Coeffs")),
    outputPropertiesPtr_(nullptr),
    vsmPtr_(nullptr),
    regionName_(dict.get<word>("region"))
{
    constructMeshObjects();
    initialise();

    if (readFields)
    {
        init(dict);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::regionFaModel::evolve()
{
    if (active_)
    {
        Info<< "\nEvolving " << modelName_ << " for region "
            << regionMesh().name() << endl;

        preEvolveRegion();

        evolveRegion();

        postEvolveRegion();

        // Provide some feedback
        if (infoOutput_)
        {
            Info<< incrIndent;
            info();
            Info<< endl << decrIndent;
        }
    }
}


void Foam::regionModels::regionFaModel::preEvolveRegion()
{}


void Foam::regionModels::regionFaModel::evolveRegion()
{}


void Foam::regionModels::regionFaModel::postEvolveRegion()
{}


Foam::scalar Foam::regionModels::regionFaModel::CourantNumber() const
{
    return 0;
}


// ************************************************************************* //
