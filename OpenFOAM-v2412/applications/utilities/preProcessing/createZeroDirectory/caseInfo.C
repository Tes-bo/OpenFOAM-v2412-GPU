/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "caseInfo.H"
#include "Time.H"
#include "boundaryInfo.H"
#include "boundaryTemplates.H"

using namespace Foam;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::caseInfo::findPatchConditionID
(
    const label patchI,
    const word& patchName
) const
{
    const wordList& patchGroups = boundaryInfo_.groups()[patchI];

    // Assign condition according to last condition applied, wins
    forAllReverse(conditionNames_, conditionI)
    {
        const wordRes& patchNames = patchNames_[conditionI];

        for (const wordRe& select : patchNames)
        {
            if (select == patchName)
            {
                // Literal match
                return conditionI;
            }
            else if (select.match(patchName))
            {
                // Regex match
                return conditionI;
            }
            else
            {
                // Check for group match
                for (const word& groupName : patchGroups)
                {
                    if (select == groupName)
                    {
                        return conditionI;
                    }
                }
            }
        }
    }

    FatalErrorInFunction
        << "Boundary patch " << patchName << " not defined"
        << exit(FatalError);

    return -1;
}


void Foam::caseInfo::updateGeometricBoundaryField()
{
    forAll(boundaryInfo_.names(), i)
    {
        const word& patchName = boundaryInfo_.names()[i];

        if (!boundaryInfo_.constraint()[i])
        {
            // Condition ID to apply to mesh boundary patch name
            const label conditionI = findPatchConditionID(i, patchName);

            const word& category = patchCategories_[conditionI];

            boundaryInfo_.setType(i, category);
        }
    }

    boundaryInfo_.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::caseInfo::caseInfo(const Time& runTime, const word& regionName)
:
    properties_
    (
        IOobject
        (
            "caseProperties",
            runTime.system(),
            regionName,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundaryInfo_(runTime, regionName),
    bcDict_(properties_.subDict("boundaryConditions")),
    conditionNames_(bcDict_.toc()),
    patchNames_(conditionNames_.size()),
    patchCategories_(conditionNames_.size()),
    patchTypes_(conditionNames_.size())
{
    // Read the (user-supplied) boundary condition information
    Info<< "    Reading case properties" << endl;

    forAll(conditionNames_, i)
    {
        const dictionary& dict = bcDict_.subDict(conditionNames_[i]);
        dict.readEntry("category", patchCategories_[i]);
        dict.readEntry("type", patchTypes_[i]);
        dict.readEntry("patches", patchNames_[i]);
    }

    updateGeometricBoundaryField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::caseInfo::checkPatches
(
    const word& regionPrefix,
    const boundaryTemplates& bcTemplates
) const
{
    // Check that all conditions have been specified correctly wrt templates
    forAll(conditionNames_, i)
    {
        bcTemplates.checkPatch
        (
            regionPrefix,
            conditionNames_[i],
            patchCategories_[i],
            patchTypes_[i]
        );
    }
}


const Foam::List<Foam::wordRes>& Foam::caseInfo::patchNames() const
{
    return patchNames_;
}


const Foam::word& Foam::caseInfo::conditionName(const label patchI) const
{
    return conditionNames_[patchI];
}


const Foam::word& Foam::caseInfo::patchCategory(const label patchI) const
{
    return patchCategories_[patchI];
}


const Foam::word& Foam::caseInfo::patchType(const label patchI) const
{
    return patchTypes_[patchI];
}


Foam::dictionary Foam::caseInfo::generateBoundaryField
(
    const word& regionPrefix,
    const word& fieldName,
    const boundaryTemplates& bcTemplates
) const
{
    dictionary boundaryField;

    forAll(boundaryInfo_.names(), j)
    {
        const word& patchName = boundaryInfo_.names()[j];

        if (boundaryInfo_.constraint()[j])
        {
            dictionary patchDict;
            patchDict.add("type", boundaryInfo_.types()[j]);

            // Add value for processor patches
            patchDict.add("value", "${:internalField}");
            boundaryField.add(patchName.c_str(), patchDict);
        }
        else
        {
            // Condition ID to apply to mesh boundary patch name
            const label conditionI = findPatchConditionID(j, patchName);

            if (conditionI == -1)
            {
                FatalErrorInFunction
                    << "Unable to find patch " << patchName
                    << " in list of boundary conditions"
                    << exit(FatalError);
            }

            const word& condition = conditionNames_[conditionI];

            const word& category = patchCategories_[conditionI];

            const word& patchType = patchTypes_[conditionI];

            dictionary optionDict;
            if (bcTemplates.optionsRequired(regionPrefix, category, patchType))
            {
                optionDict = bcDict_.subDict(condition).subDict("options");
            }

            // Create the patch dictionary entry
            dictionary patchDict
            (
                bcTemplates.generatePatchDict
                (
                    regionPrefix,
                    fieldName,
                    condition,
                    category,
                    patchType,
                    optionDict
                )
            );

            boundaryField.add(patchName.c_str(), patchDict);
        }
    }

    return boundaryField;
}


// ************************************************************************* //
