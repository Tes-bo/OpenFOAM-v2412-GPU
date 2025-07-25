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

#include "adjointEikonalSolver.H"
#include "runTimeSelectionTables.H"
#include "adjointSensitivity.H"
#include "adjointSolver.H"
#include "designVariables.H"
#include "fvOptions.H"
#include "reverseLinear.H"
#include "sensitivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointSensitivity, 0);
defineRunTimeSelectionTable(adjointSensitivity, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointSensitivity::adjointSensitivity
(
    const fvMesh& mesh,
    const dictionary& dict,
    class adjointSolver& adjointSolver
)
:
    sensitivity(mesh, dict),
    adjointSolver_(adjointSolver),
    derivatives_(0),
    suffix_(word::null),
    includeDistance_
    (
        this->dict().getOrDefault<bool>
        (
            "includeDistance",
            adjointSolver_.includeDistance()
        )
    ),
    eikonalSolver_(nullptr),
    gradDxDbMult_(nullptr),
    divDxDbMult_(nullptr),
    dxdbMult_(nullptr),
    dSfdbMult_(nullptr),
    dnfdbMult_(nullptr),
    dxdbDirectMult_(nullptr),
    pointDxDbDirectMult_(nullptr),
    bcDxDbMult_(nullptr),
    optionsDxDbMult_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<adjointSensitivity> adjointSensitivity::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    class adjointSolver& adjointSolver
)
{
    const word sensType =
        dict.optionalSubDict(mesh.name()).get<word>("sensitivityType");

    Info<< "adjointSensitivity type : " << sensType << endl;

    auto* ctorPtr = dictionaryConstructorTable(sensType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointSensitivity",
            sensType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointSensitivity>
    (
        ctorPtr(mesh, dict, adjointSolver)
    );
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool adjointSensitivity::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        // The adjoint eikonal solver requires the parameterized patches
        // as an argument, if they exist. Allocation will be managed by
        // derived classes that have access to them
        includeDistance_ = this->dict().getOrDefault<bool>
        (
            "includeDistance",
            adjointSolver_.includeDistance()
        );

        return true;
    }

    return false;
}


bool adjointSensitivity::computeDxDbInternalField() const
{
    return false;
}


void adjointSensitivity::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    derivatives_ = designVars->assembleSensitivities(*this);
}


const scalarField& adjointSensitivity::calculateSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    assembleSensitivities(designVars);
    write(type());
    return derivatives_;
}


const scalarField& adjointSensitivity::getSensitivities() const
{
    return derivatives_;
}


void adjointSensitivity::clearSensitivities()
{
    derivatives_ = Zero;
    if (fieldSensPtr_)
    {
        fieldSensPtr_().primitiveFieldRef() = Zero;
    }
    if (eikonalSolver_)
    {
        eikonalSolver_->reset();
    }
    if (adjointMeshMovementSolver_)
    {
        adjointMeshMovementSolver_->reset();
    }
}


void adjointSensitivity::write(const word& baseName)
{
    sensitivity::write(baseName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
