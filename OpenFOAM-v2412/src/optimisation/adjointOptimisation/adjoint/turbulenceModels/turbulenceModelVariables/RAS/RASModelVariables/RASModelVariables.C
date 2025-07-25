/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "RASModelVariables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModelVariables, 0);
defineRunTimeSelectionTable(RASModelVariables, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void RASModelVariables::allocateInitValues()
{
    if (solverControl_.storeInitValues())
    {
        Info<< "Storing initial values of turbulence variables" << endl;

        if (hasTMVar1())
        {
            TMVar1InitPtr_.reset
            (
                new volScalarField(TMVar1Inst().name()+"Init", TMVar1Inst())
            );
        }

        if (hasTMVar2())
        {
            TMVar2InitPtr_.reset
            (
                new volScalarField(TMVar2Inst().name()+"Init", TMVar2Inst())
            );
        }

        if (hasNut())
        {
            nutInitPtr_.reset
            (
                new volScalarField(nutRefInst().name()+"Init", nutRefInst())
            );
        }
    }
}


void RASModelVariables::allocateMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Allocating mean values of turbulence variables" << endl;

        if (hasTMVar1())
        {
            TMVar1MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        TMVar1Inst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    TMVar1Inst()
                )
            );
        }

        if (hasTMVar2())
        {
            TMVar2MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        TMVar2Inst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    TMVar2Inst()
                )
            );
        }

        if (hasNut())
        {
            nutMeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        nutRefInst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    nutRefInst()
                )
            );
        }
    }
}


Foam::refPtr<Foam::volScalarField>
RASModelVariables::cloneRefPtr(const refPtr<volScalarField>& obj) const
{
    if (obj)
    {
        const volScalarField& sf = obj();

        const word timeName = mesh_.time().timeName();

        return refPtr<volScalarField>::New(sf.name() + timeName, sf);
    }

    return nullptr;
}


void RASModelVariables::copyAndRename
(
    volScalarField& f1,
    volScalarField& f2
)
{
    f1 == f2;
    const word name1 = f1.name();
    const word name2 = f2.name();

    // Extra rename to avoid databese collision
    f2.rename("temp");
    f1.rename(name2);
    f2.rename(name1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASModelVariables::RASModelVariables
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    mesh_(mesh),
    solverControl_(SolverControl),

    TMVar1BaseName_(),
    TMVar2BaseName_(),
    nutBaseName_("nut"),

    TMVar1Ptr_(nullptr),
    TMVar2Ptr_(nullptr),
    nutPtr_(nullptr),
    distPtr_(nullptr),

    TMVar1InitPtr_(nullptr),
    TMVar2InitPtr_(nullptr),
    nutInitPtr_(nullptr),

    TMVar1MeanPtr_(nullptr),
    TMVar2MeanPtr_(nullptr),
    nutMeanPtr_(nullptr)
{}


RASModelVariables::RASModelVariables
(
    const RASModelVariables& rmv
)
:
    mesh_(rmv.mesh_),
    solverControl_(rmv.solverControl_),

    TMVar1BaseName_(rmv.TMVar1BaseName_),
    TMVar2BaseName_(rmv.TMVar2BaseName_),
    nutBaseName_(rmv.nutBaseName_),

    TMVar1Ptr_(cloneRefPtr(rmv.TMVar1Ptr_)),
    TMVar2Ptr_(cloneRefPtr(rmv.TMVar2Ptr_)),
    nutPtr_(cloneRefPtr(rmv.nutPtr_)),
    distPtr_(cloneRefPtr(rmv.distPtr_)),

    TMVar1InitPtr_(nullptr),
    TMVar2InitPtr_(nullptr),
    nutInitPtr_(nullptr),

    TMVar1MeanPtr_(nullptr),
    TMVar2MeanPtr_(nullptr),
    nutMeanPtr_(nullptr)
{}


autoPtr<RASModelVariables> RASModelVariables::clone() const
{
    return autoPtr<RASModelVariables>::New(*this);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<RASModelVariables> RASModelVariables::New
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            turbulenceModel::propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    word modelType("laminar"); // default to laminar

    const dictionary* dictptr = modelDict.findDict("RAS");

    if (dictptr)
    {
        // "RASModel" for v2006 and earlier
        dictptr->readCompat("model", {{"RASModel", -2006}}, modelType);
    }
    else
    {
        dictptr = &dictionary::null;
    }

    Info<< "Creating references for RASModel variables : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            *dictptr,
            "RASModelVariables",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<RASModelVariables>(ctorPtr(mesh, SolverControl));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> RASModelVariables::nutJacobianVar1
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "jutJacobianVar1 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    return volScalarField::New
    (
        "nutJacobianVar1",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


tmp<volScalarField> RASModelVariables::nutJacobianVar2
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "nutJacobianVar2 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    return volScalarField::New
    (
        "nutJacobianVar2",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


void RASModelVariables::restoreInitValues()
{
    if (solverControl_.storeInitValues())
    {
        if (hasTMVar1())
        {
            TMVar1Inst() == TMVar1InitPtr_();
        }
        if (hasTMVar2())
        {
            TMVar2Inst() == TMVar2InitPtr_();
        }
        if (hasNut())
        {
            nutRefInst() == nutInitPtr_();
        }
    }
}


void RASModelVariables::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Resetting mean turbulent fields to zero" << endl;

        // Reset fields to zero
        if (TMVar1MeanPtr_)
        {
            TMVar1MeanPtr_.ref() ==
                dimensionedScalar(TMVar1Inst().dimensions(), Zero);
        }
        if (TMVar2MeanPtr_)
        {
            TMVar2MeanPtr_.ref() ==
                dimensionedScalar(TMVar2Inst().dimensions(), Zero);
        }
        if (nutMeanPtr_)
        {
            nutMeanPtr_.ref() ==
                dimensionedScalar(nutRefInst().dimensions(), Zero);
        }
    }
}


void RASModelVariables::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        const label iAverageIter = solverControl_.averageIter();
        const scalar avIter(iAverageIter);
        const scalar oneOverItP1 = 1./(avIter + 1);
        const scalar mult = avIter*oneOverItP1;

        if (hasTMVar1())
        {
            TMVar1MeanPtr_.ref() ==
                (TMVar1MeanPtr_()*mult + TMVar1Inst()*oneOverItP1);
        }
        if (hasTMVar2())
        {
            TMVar2MeanPtr_.ref() ==
                (TMVar2MeanPtr_()*mult + TMVar2Inst()*oneOverItP1);
        }
        if (hasNut())
        {
            nutMeanPtr_.ref() ==
                (nutMeanPtr_()*mult + nutRefInst()*oneOverItP1);
        }
    }
}


tmp<volSymmTensorField> RASModelVariables::devReff
(
    const singlePhaseTransportModel& laminarTransport,
    const volVectorField& U
) const
{
    return volSymmTensorField::New
    (
        "devRhoReff",
        IOobject::NO_REGISTER,
      - (laminarTransport.nu() + nut())*devTwoSymm(fvc::grad(U))
    );
}


void RASModelVariables::correctBoundaryConditions
(
    const incompressible::turbulenceModel& turbulence
)
{
    if (hasTMVar1())
    {
        TMVar1Inst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar1MeanPtr_.ref().correctBoundaryConditions();
        }
    }

    if (hasTMVar2())
    {
        TMVar2Inst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar2MeanPtr_.ref().correctBoundaryConditions();
        }
    }

    if (hasNut())
    {
        nutRefInst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            nutMeanPtr_.ref().correctBoundaryConditions();
        }
    }
}


void RASModelVariables::transfer(RASModelVariables& rmv)
{
    if (rmv.hasTMVar1() && hasTMVar1())
    {
        copyAndRename(TMVar1Inst(), rmv.TMVar1Inst());
    }

    if (rmv.hasTMVar2() && hasTMVar2())
    {
        copyAndRename(TMVar2Inst(), rmv.TMVar2Inst());
    }

    if (rmv.hasNut() && hasNut())
    {
        copyAndRename(nutRefInst(), rmv.nutRefInst());
    }

    if (rmv.hasDist() && hasDist())
    {
        copyAndRename(d(), rmv.d());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
