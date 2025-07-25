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

#include "adjointRASModel.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointRASModel, 0);
defineRunTimeSelectionTable(adjointRASModel, dictionary);
addToRunTimeSelectionTable
(
    adjointTurbulenceModel,
    adjointRASModel,
    adjointTurbulenceModel
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void adjointRASModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


void adjointRASModel::restoreInitValues()
{
    const solverControl& solControl = adjointVars_.getSolverControl();
    if (solControl.storeInitValues())
    {
        if (adjointTMVariable1Ptr_)
        {
            volScalarField& var1 = adjointTMVariable1Ptr_.ref();
            var1 == dimensionedScalar(var1.dimensions(), Zero);
        }

        if (adjointTMVariable2Ptr_)
        {
            volScalarField& var2 = adjointTMVariable2Ptr_.ref();
            var2 == dimensionedScalar(var2.dimensions(), Zero);
        }
    }
}


void adjointRASModel::setMeanFields()
{
    const solverControl& solControl = adjointVars_.getSolverControl();
    if (solControl.average())
    {
        if (adjointTMVariable1Ptr_)
        {
            adjointTMVariable1MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        getAdjointTMVariable1Inst().name() + "Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    getAdjointTMVariable1Inst()
                )
            );
        }

        if (adjointTMVariable2Ptr_)
        {
            adjointTMVariable2MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        getAdjointTMVariable2Inst().name() + "Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    getAdjointTMVariable2Inst()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointRASModel::adjointRASModel
(
    const word& type,
    incompressibleVars& primalVars,
    incompressibleAdjointMeanFlowVars& adjointVars,
    objectiveManager& objManager,
    const word& adjointTurbulenceModelName
)
:
    adjointTurbulenceModel
    (
        primalVars,
        adjointVars,
        objManager,
        adjointTurbulenceModelName
    ),
    IOdictionary
    (
        IOobject
        (
            "adjointRASProperties",
            primalVars.U().time().constant(),
            primalVars.U().db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    objectiveManager_(objManager),

    adjointTurbulence_(get<word>("adjointTurbulence")),
    printCoeffs_(getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),

    y_(mesh_),

    adjointTMVariable1Ptr_(nullptr),
    adjointTMVariable2Ptr_(nullptr),
    adjointTMVariablesBaseNames_(0),
    adjointTMVariable1MeanPtr_(nullptr),
    adjointTMVariable2MeanPtr_(nullptr),
    adjMomentumBCSourcePtr_( createZeroBoundaryPtr<vector>(mesh_) ),
    wallShapeSensitivitiesPtr_( createZeroBoundaryPtr<vector>(mesh_) ),
    wallFloCoSensitivitiesPtr_( createZeroBoundaryPtr<vector>(mesh_) ),
    includeDistance_(false),
    changedPrimalSolution_(true)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<adjointRASModel> adjointRASModel::New
(
    incompressibleVars& primalVars,
    incompressibleAdjointMeanFlowVars& adjointVars,
    objectiveManager& objManager,
    const word& adjointTurbulenceModelName
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "adjointRASProperties",
            primalVars.U().time().constant(),
            primalVars.U().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word modelType(dict.get<word>("adjointRASModel"));

    Info<< "Selecting adjointRAS turbulence model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointRASModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointRASModel>
    (
        ctorPtr
        (
            primalVars,
            adjointVars,
            objManager,
            adjointTurbulenceModelName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjointRASModel::correct()
{
    adjointTurbulenceModel::correct();

    if (adjointTurbulence_ && mesh_.changing())
    {
        y_.correct();
    }
}


bool adjointRASModel::read()
{
    //if (regIOobject::read())

    // Bit of trickery : we are both IOdictionary ('adjointRASProperties') and
    // an regIOobject from the adjointTurbulenceModel level. Problem is to
    // distinguish between the two - we only want to reread the IOdictionary.

    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            IOdictionary::type()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        readEntry("adjointTurbulence", adjointTurbulence_);

        if (const dictionary* dictPtr = findDict(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        return true;
    }
    else
    {
        return false;
    }
}


volScalarField& adjointRASModel::getAdjointTMVariable1Inst()
{
    if (!adjointTMVariable1Ptr_)
    {
        // if pointer is not set, set it to a zero field
        adjointTMVariable1Ptr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "adjointTMVariable1" + type(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }

    return adjointTMVariable1Ptr_();
}


volScalarField& adjointRASModel::getAdjointTMVariable2Inst()
{
    if (!adjointTMVariable2Ptr_)
    {
        // if pointer is not set, set it to a zero field
        adjointTMVariable2Ptr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                   "adjointTMVariable2" + type(),
                   mesh_.time().timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }

    return adjointTMVariable2Ptr_();
}


volScalarField& adjointRASModel::getAdjointTMVariable1()
{
    if (adjointVars_.getSolverControl().useAveragedFields())
    {
        return adjointTMVariable1MeanPtr_();
    }
    else
    {
        return getAdjointTMVariable1Inst();
    }
}



volScalarField& adjointRASModel::getAdjointTMVariable2()
{
    if (adjointVars_.getSolverControl().useAveragedFields())
    {
        return adjointTMVariable2MeanPtr_();
    }
    else
    {
        return getAdjointTMVariable2Inst();
    }
}


autoPtr<volScalarField>& adjointRASModel::getAdjointTMVariable1InstPtr()
{
    return adjointTMVariable1Ptr_;
}


autoPtr<volScalarField>& adjointRASModel::getAdjointTMVariable2InstPtr()
{
    return adjointTMVariable2Ptr_;
}


const wordList& adjointRASModel::getAdjointTMVariablesBaseNames() const
{
    return adjointTMVariablesBaseNames_;
}


tmp<volScalarField> adjointRASModel::nutJacobianTMVar1() const
{
    dimensionSet dims
    (
        bool(adjointTMVariable1Ptr_)
      ? dimViscosity/adjointTMVariable1Ptr_().dimensions()
      : dimless
    );

    return
        tmp<volScalarField>::New
        (
            IOobject
            (
                "nutJacobianTMVar1"+type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dims, Zero)
        );
}


tmp<volScalarField> adjointRASModel::nutJacobianTMVar2() const
{
    dimensionSet dims
    (
        bool(adjointTMVariable2Ptr_)
      ? dimViscosity/adjointTMVariable2Ptr_().dimensions()
      : dimless
    );

    return
        tmp<volScalarField>::New
        (
            IOobject
            (
                "nutJacobianTMVar2"+type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dims, Zero)
        );
}


tmp<volVectorField> adjointRASModel::nutJacobianU
(
    tmp<volScalarField>& dNutdUMult
) const
{
    // Deliberately returning a null pointer in the base class
    return nullptr;
}


tmp<scalarField> adjointRASModel::diffusionCoeffVar1(label patchI) const
{
    return tmp<scalarField>::New(mesh_.boundary()[patchI].size(), Zero);
}


tmp<scalarField> adjointRASModel::diffusionCoeffVar2(label patchI) const
{
    return tmp<scalarField>::New(mesh_.boundary()[patchI].size(), Zero);
}


void adjointRASModel::setChangedPrimalSolution()
{
    changedPrimalSolution_ = true;
}


void adjointRASModel::resetMeanFields()
{
    const solverControl& solControl = adjointVars_.getSolverControl();
    if (solControl.average())
    {
        if (adjointTMVariable1MeanPtr_)
        {
            adjointTMVariable1MeanPtr_() ==
                dimensionedScalar(adjointTMVariable1Ptr_().dimensions(), Zero);
        }
        if (adjointTMVariable2MeanPtr_)
        {
            adjointTMVariable2MeanPtr_() ==
                dimensionedScalar(adjointTMVariable2Ptr_().dimensions(), Zero);
        }
    }
}


void adjointRASModel::computeMeanFields()
{
    const solverControl& solControl = adjointVars_.getSolverControl();
    if (solControl.doAverageIter())
    {
        const label iAverageIter = solControl.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter+1);
        scalar mult = avIter*oneOverItP1;
        if (adjointTMVariable1MeanPtr_)
        {
            adjointTMVariable1MeanPtr_() ==
                adjointTMVariable1Ptr_()*mult
              + getAdjointTMVariable1Inst()*oneOverItP1;
        }
        if (adjointTMVariable2MeanPtr_)
        {
            adjointTMVariable2MeanPtr_() ==
                adjointTMVariable2Ptr_()*mult
              + getAdjointTMVariable2Inst()*oneOverItP1;
        }
    }
}


bool adjointRASModel::includeDistance() const
{
    return includeDistance_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
