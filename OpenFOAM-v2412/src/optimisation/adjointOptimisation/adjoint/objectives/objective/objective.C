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

#include "objective.H"
#include "createZeroField.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objective, 0);
defineRunTimeSelectionTable(objective, objective);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void objective::makeFolder()
{
    if (Pstream::master())
    {
        const Time& time = mesh_.time();
        objFunctionFolder_ =
            time.globalPath()/"optimisation"/type()/time.timeName();
        if (mesh_.name() != polyMesh::defaultRegion)
        {
            objFunctionFolder_ /= mesh_.name();
        }

        mkDir(objFunctionFolder_);
    }
}


void objective::setObjectiveFilePtr() const
{
    objFunctionFilePtr_.reset
    (
        new OFstream(objFunctionFolder_/objectiveName_ + adjointSolverName_)
    );
}


void objective::setInstantValueFilePtr() const
{
    instantValueFilePtr_.reset
    (
        new OFstream
        (
            objFunctionFolder_/objectiveName_ + "Instant" + adjointSolverName_
        )
    );
}


void objective::setMeanValueFilePtr() const
{
    meanValueFilePtr_.reset
    (
        new OFstream
        (
            objFunctionFolder_/objectiveName_ + "Mean" + adjointSolverName_
        )
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const dictionary& objective::dict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objective::objective
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    localIOdictionary
    (
        IOobject
        (
            adjointSolverName + "_" + dict.dictName(),
            mesh.time().timeName(),
            fileName("uniform")/fileName("objectives")/adjointSolverName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        // avoid type checking since dictionary is read using the
        // derived type name and type() will result in "objective"
        // here
        word::null
    ),
    mesh_(mesh),
    dict_(dict),
    adjointSolverName_(adjointSolverName),
    primalSolverName_(primalSolverName),
    objectiveName_(dict.dictName()),
    computeMeanFields_(false), // is reset in derived classes
    nullified_(false),
    normalize_
    (
        dict.getOrDefaultCompat<bool>
        (
            "normalise", {{"normalize", 2406}}, false
        )
    ),
    shouldWrite_(true),

    J_(Zero),
    JMean_(this->getOrDefault<scalar>("JMean", Zero)),
    weight_(dict.get<scalar>("weight")),
    computed_(false),
    normFactor_(nullptr),
    target_
    (
        dict.found("target") ?
        autoPtr<scalar>::New(dict.get<scalar>("target")) :
        nullptr
    ),
    targetLeft_
    (
        dict.found("targetLeft") ?
        autoPtr<scalar>::New(dict.get<scalar>("targetLeft")) :
        nullptr
    ),
    integrationStartTimePtr_(nullptr),
    integrationEndTimePtr_(nullptr),
    fieldNames_(),

    // Initialize pointers to nullptr.
    // Not all of them are required for each objective function.
    // Each child should allocate whatever is needed.

    dJdbPtr_(nullptr),
    dJdbFieldPtr_(nullptr),
    bdJdbPtr_(nullptr),
    bdSdbMultPtr_(nullptr),
    bdndbMultPtr_(nullptr),
    bdxdbMultPtr_(nullptr),
    bdxdbDirectMultPtr_(nullptr),
    bEdgeContribution_(nullptr),
    divDxDbMultPtr_(nullptr),
    gradDxDbMultPtr_(nullptr),

    objFunctionFolder_("word"),
    objFunctionFilePtr_(nullptr),
    instantValueFilePtr_(nullptr),
    meanValueFilePtr_(nullptr),
    width_(IOstream::defaultPrecision() + 5)
{
    makeFolder();
    // Read integration start and end times, if present.
    // For unsteady runs only
    if (dict.found("integrationStartTime"))
    {
        integrationStartTimePtr_.reset
        (
            new scalar(dict.get<scalar>("integrationStartTime"))
        );
    }
    if (dict.found("integrationEndTime"))
    {
        integrationEndTimePtr_.reset
        (
            new scalar(dict.get<scalar>("integrationEndTime"))
        );
    }

    // Set normalization factor, if present
    if (normalize_)
    {
        scalar normFactor(Zero);
        if (dict.readIfPresent("normFactor", normFactor))
        {
            normFactor_.reset(new scalar(normFactor));
        }
        else if (this->readIfPresent("normFactor", normFactor))
        {
            normFactor_.reset(new scalar(normFactor));
        }
    }

    // Set the weight factor in case of continuation
    this->readIfPresent("weight", weight_);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<objective> objective::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& objectiveType,
    const word& adjointSolverName,
    const word& primalSolverName
)
{
    auto* ctorPtr = objectiveConstructorTable(objectiveType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "objective",
            objectiveType,
            *objectiveConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<objective>
    (
        ctorPtr
        (
            mesh,
            dict,
            adjointSolverName,
            primalSolverName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool objective::readDict(const dictionary& dict)
{
    dict_ = dict;
    return true;
}


scalar objective::JCycle(bool negate) const
{
    scalar J(J_);
    if
    (
        computeMeanFields_
     || (hasIntegrationStartTime() && hasIntegrationEndTime())
    )
    {
        J = JMean_;
    }

    // Subtract target, in case the objective is used as a constraint
    if (target_)
    {
        if (negate)
        {
            J = - J + targetLeft_();
        }
        else
        {
            J -= target_();
        }
    }

    // Normalize here, in order to get the correct value for line search
    if (normalize_ && normFactor_)
    {
        J /= normFactor_();
    }
    J *= weight_;

    return J;
}


void objective::updateNormalizationFactor()
{
    if (normalize_ && !normFactor_)
    {
        scalar J(JCycle()/weight_);
        normFactor_.reset(new scalar(J));
        DebugInfo
            << "objective " << name() << ":: updating norm factor "
            << "to " << normFactor_()
            << " for time = " << mesh_.time().timeName() << endl;
    }
}


void objective::accumulateJMean(solverControl& solverControl)
{
    if (solverControl.doAverageIter())
    {
        const label iAverageIter = solverControl.averageIter();
        if (iAverageIter == 0)
        {
            JMean_ = Zero;
        }
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter + 1);
        scalar mult = avIter*oneOverItP1;
        JMean_ = JMean_*mult + J_*oneOverItP1;
    }
}


void objective::accumulateJMean()
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        const scalar time = mesh_.time().value();
        if (isWithinIntegrationTime())
        {
            const scalar dt = mesh_.time().deltaTValue();
            const scalar elapsedTime = time - integrationStartTimePtr_();
            const scalar denom = elapsedTime + dt;
            JMean_ = (JMean_*elapsedTime + J_*dt)/denom;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time"
            << exit(FatalError);
    }
}


scalar objective::weight() const
{
    return weight_;
}


bool objective::normalize() const
{
    return normalize_;
}


void objective::doNormalization()
{
    if (normalize_ && normFactor_)
    {
        const scalar oneOverNorm(1./normFactor_());

        if (hasdJdb())
        {
            dJdbPtr_().primitiveFieldRef() *= oneOverNorm;
        }
        if (hasdJdbField())
        {
            dJdbFieldPtr_() *= oneOverNorm;
        }
        if (hasBoundarydJdb())
        {
            bdJdbPtr_() *= oneOverNorm;
        }
        if (hasdSdbMult())
        {
            bdSdbMultPtr_() *= oneOverNorm;
        }
        if (hasdndbMult())
        {
            bdndbMultPtr_() *= oneOverNorm;
        }
        if (hasdxdbMult())
        {
            bdxdbMultPtr_() *= oneOverNorm;
        }
        if (hasdxdbDirectMult())
        {
            bdxdbDirectMultPtr_() *= oneOverNorm;
        }
        if (hasBoundaryEdgeContribution())
        {
            bEdgeContribution_() *= oneOverNorm;
        }
        if (hasDivDxDbMult())
        {
            divDxDbMultPtr_() *= oneOverNorm;
        }
        if (hasGradDxDbMult())
        {
            gradDxDbMultPtr_() *= oneOverNorm;
        }
    }
}


bool objective::isWithinIntegrationTime() const
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        const scalar time = mesh_.time().value();
        return
            (
                time >= integrationStartTimePtr_()
             && time <= integrationEndTimePtr_()
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time for objective '"
            << objectiveName_ << "'"
            << exit(FatalError);
    }
    return false;
}


void objective::incrementIntegrationTimes(const scalar timeSpan)
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        integrationStartTimePtr_() += timeSpan;
        integrationEndTimePtr_() += timeSpan;
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time"
            << exit(FatalError);
    }
}


void objective::update()
{
    // Objective function value
    J();

    // volFields
    update_dJdb();
    update_dJdbField();
    update_divDxDbMultiplier();
    update_gradDxDbMultiplier();

    // boundaryFields
    update_boundarydJdb();
    update_dSdbMultiplier();
    update_dndbMultiplier();
    update_dxdbMultiplier();
    update_dxdbDirectMultiplier();
    update_boundaryEdgeContribution();
}


void objective::nullify()
{
    if (!nullified_)
    {
        if (hasdJdb())
        {
            dJdbPtr_() == dimensionedScalar(dJdbPtr_().dimensions(), Zero);
        }
        if (hasdJdbField())
        {
            dJdbFieldPtr_() = Zero;
        }
        if (hasBoundarydJdb())
        {
            bdJdbPtr_() == vector::zero;
        }
        if (hasdSdbMult())
        {
            bdSdbMultPtr_() == vector::zero;
        }
        if (hasdndbMult())
        {
            bdndbMultPtr_() == vector::zero;
        }
        if (hasdxdbMult())
        {
            bdxdbMultPtr_() == vector::zero;
        }
        if (hasdxdbDirectMult())
        {
            bdxdbDirectMultPtr_() == vector::zero;
        }
        if (hasBoundaryEdgeContribution())
        {
            for (Field<vectorField>& field : bEdgeContribution_())
            {
                field = vector::zero;
            }
        }
        if (hasDivDxDbMult())
        {
            divDxDbMultPtr_() ==
                dimensionedScalar(divDxDbMultPtr_().dimensions(), Zero);
        }
        if (hasGradDxDbMult())
        {
            gradDxDbMultPtr_() ==
                dimensionedTensor(gradDxDbMultPtr_().dimensions(), Zero);
        }

        nullified_ = true;
    }
}


bool objective::write(const bool valid) const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        if (!objFunctionFilePtr_)
        {
            setObjectiveFilePtr();
            OFstream& file = objFunctionFilePtr_();
            file.setf(std::ios_base::left);

            if (target_)
            {
                file<< setw(width_) << "#target" << " "
                    << setw(width_) << target_() << endl;
            }
            if (targetLeft_)
            {
                file<< setw(width_) << "#targetLeft" << " "
                    << setw(width_) << targetLeft_() << endl;
            }
            if (normalize_)
            {
                file<< setw(width_) << "#normFactor " << " "
                    << setw(width_) << normFactor_() << endl;
            }
            addHeaderInfo();
            file<< setw(4) << "#" << " ";
            file<< setw(width_) << "J" << " ";
            file<< setw(width_) << "JCycle" << " ";
            if (targetLeft_)
            {
                file<< setw(width_) << "JCycleLeft" << " ";
            }
            addHeaderColumns();
            file<< endl;
        }

        OFstream& file = objFunctionFilePtr_();
        file<< setw(4) << mesh_.time().value() << " ";
        file<< setw(width_) << J_ << " ";
        file<< setw(width_) << JCycle() << " ";
        if (targetLeft_)
        {
            file<< setw(width_) << JCycle(true) << " ";
        }
        addColumnValues();
        file<< endl;
    }

    return true;
}


void objective::writeInstantaneousValue() const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        unsigned int width = IOstream::defaultPrecision() + 6;
        if (!instantValueFilePtr_)
        {
            setInstantValueFilePtr();
        }

        instantValueFilePtr_()
            << setw(width) << mesh_.time().value() << tab << J_ << endl;
    }
}


void objective::writeInstantaneousSeparator() const
{
    if (Pstream::master())
    {
        if (instantValueFilePtr_)
        {
            instantValueFilePtr_() << endl;
        }
    }
}


void objective::writeMeanValue() const
{
    if (Pstream::master())
    {
        // Write mean value if necessary
        // Covers both steady and unsteady runs
        if
        (
            computeMeanFields_
         || (hasIntegrationStartTime() && hasIntegrationEndTime())
        )
        {
            // File is opened only upon invocation of the write function
            // in order to avoid various instantiations of the same objective
            // opening the same file
            if (!meanValueFilePtr_)
            {
                setMeanValueFilePtr();
            }

            meanValueFilePtr_()
                << mesh_.time().value() << tab << JMean_ << endl;
        }
    }
}


bool objective::writeData(Ostream& os) const
{
    os.writeEntry("JMean", JMean_);
    if (normFactor_)
    {
        os.writeEntry("normFactor", normFactor_());
    }
    os.writeEntry("weight", weight_);
    return os.good();
}


void objective::addHeaderInfo() const
{
    // Does nothing in base
}


void objective::addHeaderColumns() const
{
    // Does nothing in base
}


void objective::addColumnValues() const
{
    // Does nothing in base
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
