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

#include "objectiveManager.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveManager, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveManager::objectiveManager
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    regIOobject
    (
        IOobject
        (
            "objectiveManager" + adjointSolverName,
            mesh.time().system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    mesh_(mesh),
    dict_(dict),
    adjointSolverName_(adjointSolverName),
    primalSolverName_(primalSolverName),
    objectives_(0),
    weightedObjectiveFile_(nullptr)
{
    // Construct objectives
    //~~~~~~~~~~~~~~~~~~~~~
    Info << "Constructing objective functions " << nl << endl;
    const word objectiveType = dict.get<word>("type");
    const dictionary& objectiveNamesDict(dict.subDict("objectiveNames"));
    wordList objectiveNames(objectiveNamesDict.toc());
    objectives_.setSize(objectiveNames.size());

    forAll(objectiveNames, objectivei)
    {
        const word& objectiveName = objectiveNames[objectivei];

        objectives_.set
        (
            objectivei,
            objective::New
            (
                mesh_,
                objectiveNamesDict.subDict(objectiveName),
                objectiveType,
                adjointSolverName,
                primalSolverName
            )
        );
    }

    if (objectives_.empty())
    {
        FatalIOErrorInFunction(objectiveNamesDict)
            << "No objectives have been set - cannot perform an optimisation"
            << exit(FatalIOError);
    }

    if (Pstream::master())
    {
        if (objectives_.size() > 1)
        {
            const Time& time = mesh_.time();
            weightedObjectiveFile_.reset
            (
                new OFstream
                (
                    time.globalPath()/"optimisation"/"objective"
                        /time.timeName()/"weightedObjective"+adjointSolverName_
                )
            );

            unsigned int width = IOstream::defaultPrecision() + 5;
            weightedObjectiveFile_()
                << setw(4) << "#" << " "
                << setw(width) << "weightedObjective" << " ";
            for (objective& objI : objectives_)
            {
                weightedObjectiveFile_()
                    << setw(width) << objI.objectiveName() << " ";
            }
            weightedObjectiveFile_()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool objectiveManager::readDict(const dictionary& dict)
{
    dict_ = dict;
    for (objective& obj : objectives_)
    {
        obj.readDict
        (
            dict.subDict("objectiveNames").subDict(obj.objectiveName())
        );
    }

    return true;
}


void objectiveManager::updateNormalizationFactor()
{
    // Update normalization factors for all objectives
    for (objective& obj : objectives_)
    {
        if (obj.normalize())
        {
            obj.updateNormalizationFactor();
        }
    }
}


void objectiveManager::update()
{
    // Update all fields related to the objective function
    for (objective& obj : objectives_)
    {
        obj.update();
    }
}


void objectiveManager::updateOrNullify()
{
    // Update contributions to adjoint if true, otherwise return nulls
    for (objective& obj : objectives_)
    {
        if (obj.isWithinIntegrationTime())
        {
            obj.update();
        }
        else
        {
            obj.nullify();
        }
    }
}


void objectiveManager::incrementIntegrationTimes(const scalar timeSpan)
{
    // Update start and end integration times by adding the timeSpan
    // of the optimisation cycle
    for (objective& obj : objectives_)
    {
        obj.incrementIntegrationTimes(timeSpan);
    }
}


scalar objectiveManager::print(bool negate)
{
    scalar objValue(Zero);
    Info<< "Adjoint solver " << adjointSolverName_ << endl;
    for (objective& obj : objectives_)
    {
        // This function is used to obtain the value used to figure out if
        // line search is converged or not. If the objective is not updated
        // in each iteration of the primal solver, the old objective value
        // might be returned. Force the update of the objective the next
        // time objective::J() is called.
        obj.setComputed(false);
        const scalar cost = obj.JCycle(negate);
        objValue += cost;

        Info<< obj.objectiveName() << " : " << cost << endl;
    }

    Info<< "Weighted objective : " << objValue << nl << endl;

    return objValue;
}


void objectiveManager::setWrite(const bool shouldWrite)
{
    for (objective& obj : objectives_)
    {
        obj.setWrite(shouldWrite);
    }
}


bool objectiveManager::writeObjectives()
{
    for (const objective& obj : objectives_)
    {
        // Write objective function to file
        if (obj.shouldWrite())
        {
            obj.write();
            obj.writeMeanValue();
        }
    }

    return true;
}


bool objectiveManager::writeObjectives
(
    const scalar weightedObjective,
    const bool valid
)
{
    if (weightedObjectiveFile_)
    {
        auto& os = weightedObjectiveFile_();

        unsigned int width = IOstream::defaultPrecision() + 5;

        os
            << setw(4) << mesh_.time().timeName() << " "
            << setw(width) << weightedObjective << " ";

        for (objective& objI : objectives_)
        {
            os  << setw(width) << objI.JCycle() << " ";
        }
        os  << endl;
    }

    return writeObjectives();
}


void objectiveManager::updateAndWrite()
{
    updateNormalizationFactor();
    update();
    scalar weightedObjective = print();
    writeObjectives(weightedObjective);
}


PtrList<objective>& objectiveManager::getObjectiveFunctions()
{
    return objectives_;
}


const PtrList<objective>& objectiveManager::getObjectiveFunctions() const
{
    return objectives_;
}


const word& objectiveManager::adjointSolverName() const
{
    return adjointSolverName_;
}


const word& objectiveManager::primalSolverName() const
{
    return primalSolverName_;
}


void objectiveManager::checkIntegrationTimes() const
{
    for (const objective& obj : objectives_)
    {
        if (!obj.hasIntegrationStartTime() || !obj.hasIntegrationEndTime())
        {
            FatalErrorInFunction()
                << "Objective function " << obj.objectiveName()
                << " does not have a defined integration start or end time "
                << exit(FatalError);
        }
    }
}


void objectiveManager::addSource(fvVectorMatrix& matrix)
{
    for (objective& obj : objectives_)
    {
        obj.addSource(matrix);
    }
}


void objectiveManager::addSource(fvScalarMatrix& matrix)
{
    for (objective& obj : objectives_)
    {
        obj.addSource(matrix);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
