/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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

#include "incompressibleAdjointMeanFlowVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleAdjointMeanFlowVars, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void incompressibleAdjointMeanFlowVars::setFields()
{
    setField(paPtr_, mesh_, "pa", solverName_, useSolverNameForFields_);
    setField(UaPtr_, mesh_, "Ua", solverName_, useSolverNameForFields_);
    setFluxField
    (
        phiaPtr_,
        mesh_,
        UaInst(),
        "phia",
        solverName_,
        useSolverNameForFields_
    );

    mesh_.setFluxRequired(paPtr_->name());
}


void incompressibleAdjointMeanFlowVars::setMeanFields()
{
    // Allocate mean fields
    // Only mean flow here since turbulent quantities
    // are allocated automatically in RASModelVariables
    if (solverControl_.average())
    {
        Info<< "Allocating Mean Adjoint Fields" << endl;
        paMeanPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    paInst().name() + "Mean",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                paInst()
            )
        );
        UaMeanPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    UaInst().name() + "Mean",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                UaInst()
            )
        );
        phiaMeanPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiaInst().name() + "Mean",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                phiaInst()
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleAdjointMeanFlowVars::incompressibleAdjointMeanFlowVars
(
    fvMesh& mesh,
    solverControl& SolverControl,
    incompressibleVars& primalVars
)
:
    variablesSet(mesh, SolverControl.solverDict()),
    solverControl_(SolverControl),
    primalVars_(primalVars),
    paPtr_(nullptr),
    UaPtr_(nullptr),
    phiaPtr_(nullptr),
    paMeanPtr_(nullptr),
    UaMeanPtr_(nullptr),
    phiaMeanPtr_(nullptr)
{
    setFields();
    setMeanFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const incompressibleVars& incompressibleAdjointMeanFlowVars::primalVars() const
{
    return primalVars_;
}

const volScalarField& incompressibleAdjointMeanFlowVars::pa() const
{
    if (solverControl_.useAveragedFields())
    {
        return paMeanPtr_();
    }
    else
    {
        return paPtr_();
    }
}


volScalarField& incompressibleAdjointMeanFlowVars::pa()
{
    if (solverControl_.useAveragedFields())
    {
        return paMeanPtr_();
    }
    else
    {
        return paPtr_();
    }
}


const volVectorField& incompressibleAdjointMeanFlowVars::Ua() const
{
    if (solverControl_.useAveragedFields())
    {
        return UaMeanPtr_();
    }
    else
    {
        return UaPtr_();
    }
}


volVectorField& incompressibleAdjointMeanFlowVars::Ua()
{
    if (solverControl_.useAveragedFields())
    {
        return UaMeanPtr_();
    }
    else
    {
        return UaPtr_();
    }
}


const surfaceScalarField& incompressibleAdjointMeanFlowVars::phia() const
{
    if (solverControl_.useAveragedFields())
    {
        return phiaMeanPtr_();
    }
    else
    {
        return phiaPtr_();
    }
}


surfaceScalarField& incompressibleAdjointMeanFlowVars::phia()
{
    if (solverControl_.useAveragedFields())
    {
        return phiaMeanPtr_();
    }
    else
    {
        return phiaPtr_();
    }
}


bool incompressibleAdjointMeanFlowVars::computeMeanFields() const
{
    return solverControl_.average();
}


const solverControl& incompressibleAdjointMeanFlowVars::getSolverControl() const
{
    return solverControl_;
}


void incompressibleAdjointMeanFlowVars::nullify()
{
    variablesSet::nullifyField(paPtr_());
    variablesSet::nullifyField(UaPtr_());
    variablesSet::nullifyField(phiaPtr_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
