/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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

#include "SQPBase.H"
#include "IOmanip.H"
#include "updateMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SQPBase, 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SQPBase::SQPBase
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const updateMethod& UpdateMethod,
    const word& type
)
:
    constrainedOptimisationMethod
    (
        mesh,
        dict,
        designVars,
        UpdateMethod.nConstraints(),
        type
    ),
    LagrangianDerivatives_(designVars().getVars().size(), Zero),
    constraintDerivativesOld_
    (
        UpdateMethod.nConstraints(),
        scalarField(LagrangianDerivatives_.size(), Zero)
    ),
    lamdas_
    (
        UpdateMethod.found("lamdas") ?
        scalarField("lamdas", UpdateMethod, UpdateMethod.nConstraints()) :
        scalarField(UpdateMethod.nConstraints(), Zero)
    ),
    objFunctionFolder_
    (
        mesh.time().globalPath()/"optimisation"/"objective"/
        mesh.time().timeName()
    ),
    meritFunctionFile_(nullptr),
    mu_(Zero),
    delta_
    (
        UpdateMethod.coeffsDict(type).getOrDefault<scalar>("delta", 0.1)
    )
{
    // Read in old constraint derivatives if present
    forAll(lamdas_, cI)
    {
        if (UpdateMethod.found("constraintDerivativesOld" + Foam::name(cI)))
        {
            constraintDerivativesOld_[cI] =
                scalarField
                (
                    "constraintDerivativesOld" + Foam::name(cI),
                    UpdateMethod,
                    LagrangianDerivatives_.size()
                );
        }
    }
    // Create folder to merit function
    if (Pstream::master())
    {
        mkDir(objFunctionFolder_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SQPBase::addToFile(Ostream& os) const
{
    forAll(constraintDerivativesOld_, cI)
    {
        constraintDerivativesOld_[cI].
            writeEntry("constraintDerivativesOld" + Foam::name(cI), os);
    }
    lamdas_.writeEntry("lamdas", os);

    return true;
}


bool Foam::SQPBase::writeMeritFunction(const updateMethod& UpdateMethod)
{
    scalar objectivePart = UpdateMethod.getObjectiveValue();
    scalar constraintPart = mu_*meritFunctionConstraintPart();
    scalar merit = objectivePart + constraintPart;
    const scalarField& cValues = UpdateMethod.getConstraintValues();
    if (Pstream::master())
    {
        unsigned int width = IOstream::defaultPrecision() + 6;
        unsigned int constraintsSize = lamdas_.size();
        constraintsSize = constraintsSize*(width + 1) + 2;

        // Open file and write header
        if (!meritFunctionFile_)
        {
            meritFunctionFile_.reset
            (
                new OFstream(objFunctionFolder_/word("meritFunction"))
            );

            meritFunctionFile_()
                << setw(1) << "#" << " "
                << setw(width) << "merit" << " "
                << setw(width) << "J" << " "
                << setw(constraintsSize) << "lamdas" << " "
                << setw(constraintsSize) << "constraints" << " "
                << setw(width) << "mu" << " "
                << setw(width) << "constraintContr" << endl;

        }

        meritFunctionFile_()
            << setw(1) << UpdateMethod.getCycle() << " "
            << setw(width) << merit << " "
            << setw(width) << objectivePart << " "
            << setw(1) << "(";

        forAll(lamdas_, cI)
        {
            meritFunctionFile_()
                << setw(width) << lamdas_[cI] << setw(1) << " ";
        }
        meritFunctionFile_() << setw(3) << ")(";
        forAll(cValues, cI)
        {
            meritFunctionFile_()
                << setw(width) << cValues[cI] << setw(1) << " ";
        }
        meritFunctionFile_() << setw(2) << ") ";
        meritFunctionFile_() << setw(width) << mu_ << " ";
        meritFunctionFile_() << setw(width) << constraintPart << endl;
    }
    return true;
}


// ************************************************************************* //
