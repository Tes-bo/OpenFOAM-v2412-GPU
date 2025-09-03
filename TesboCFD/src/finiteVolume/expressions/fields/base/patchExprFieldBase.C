/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider
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

#include "patchExprFieldBase.H"
#include "facePointPatch.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "pointMesh.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::expressions::patchExprFieldBase::readExpressions
(
    const dictionary& dict,
    enum expectedTypes expectedType,
    bool wantPointData
)
{
    if (debug_)
    {
        Info<< "Expression BC with " << dict << nl;
    }

    valueExpr_.clear();
    gradExpr_.clear();
    fracExpr_.clear();

    if (expectedTypes::VALUE_TYPE == expectedType)
    {
        // Mandatory
        valueExpr_.readEntry("valueExpr", dict);
    }
    else if (expectedTypes::GRADIENT_TYPE == expectedType)
    {
        // Mandatory
        gradExpr_.readEntry("gradientExpr", dict);
    }
    else
    {
        // MIXED_TYPE
        const bool evalValue = valueExpr_.readIfPresent("valueExpr", dict);
        const bool evalGrad = gradExpr_.readIfPresent("gradientExpr", dict);

        // Expect a fraction as well
        // - but allow it to be optional and defer treatment to inherited BC

        if (evalValue && evalGrad)
        {
            if
            (
                fracExpr_.readIfPresent("fractionExpr", dict)
             && !fracExpr_.empty()
            )
            {
                // Add function call wrapping for point data,
                // but not for 0/1 (handled as shortcuts later)
                if (wantPointData && fracExpr_ != "0" && fracExpr_ != "1")
                {
                    fracExpr_ = "point(" + fracExpr_ + ")";
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Entries 'valueExpr' and 'gradientExpr' "
                   "(mixed-conditon) not found in dictionary "
                << dict.name() << nl
                << exit(FatalIOError);
        }

        if (debug_)
        {
            if (!evalValue)
            {
                Info<< "Mixed with no valueExpr" << nl;
            }
            if (!evalGrad)
            {
                Info<< "Mixed with no gradientExpr" << nl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::patchExprFieldBase::patchExprFieldBase()
:
    debug_(false),
    evalOnConstruct_(false),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{}


Foam::expressions::patchExprFieldBase::patchExprFieldBase
(
    const dictionary& dict,
    enum expectedTypes expectedType,
    bool wantPointData
)
:
    debug_(dict.getOrDefault("debug", false)),
    evalOnConstruct_(dict.getOrDefault("evalOnConstruct", false)),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{
    readExpressions(dict, expectedType, wantPointData);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::patchExprFieldBase::write(Ostream& os) const
{
    os.writeEntryIfDifferent<bool>("evalOnConstruct", false, evalOnConstruct_);

    // Do not emit debug_ value

    // Write expression, but not empty ones
    valueExpr_.writeEntry("valueExpr", os, false);
    gradExpr_.writeEntry("gradientExpr", os, false);
    fracExpr_.writeEntry("fractionExpr", os, false);
}


// ************************************************************************* //
