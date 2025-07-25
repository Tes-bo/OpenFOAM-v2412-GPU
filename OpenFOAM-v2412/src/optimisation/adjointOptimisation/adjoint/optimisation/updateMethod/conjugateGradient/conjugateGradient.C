/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "conjugateGradient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conjugateGradient, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        conjugateGradient,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conjugateGradient::conjugateGradient
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    updateMethod(mesh, dict, designVars, nConstraints, type),

    dxOld_(readOrZeroField("dxOld", activeDesignVars_.size())),
    sOld_(readOrZeroField("sOld", activeDesignVars_.size())),
    betaType_(coeffsDict(type).getOrDefault<word>("betaType", "FletcherReeves"))
{
    // Check if beta type is valid
    if
    (
       !(betaType_ == "FletcherReeves")
    && !(betaType_ == "PolakRibiere")
    && !(betaType_ == "PolakRibiereRestarted")
    )
    {
        FatalErrorInFunction
           << "Invalid betaType " << betaType_ << ". Valid options are "
           << "FletcherReeves, PolakRibiere, PolakRibiereRestarted"
           << nl << nl
           << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conjugateGradient::computeCorrection()
{
    if (counter_ == 0)
    {
        Info<< "Using steepest descent for the first iteration" << endl;
        for (const label varI : activeDesignVars_)
        {
            correction_[varI] = -eta_*objectiveDerivatives_[varI];
        }

        dxOld_.map(-objectiveDerivatives_, activeDesignVars_);
        sOld_ = dxOld_;
    }
    else
    {
        scalarField dx(activeDesignVars_.size(), Zero);
        dx.map(-objectiveDerivatives_, activeDesignVars_);

        scalar beta(Zero);
        if (betaType_ == "FletcherReeves")
        {
            beta = globalSum(dx*dx)/globalSum(dxOld_ * dxOld_);
        }
        else if (betaType_ == "PolakRibiere")
        {
            beta = globalSum(dx*(dx - dxOld_))/globalSum(dxOld_ * dxOld_);
        }
        else if (betaType_ == "PolakRibiereRestarted")
        {
            beta =
            max
            (
                scalar(0),
                globalSum(dx*(dx - dxOld_))/globalSum(dxOld_ * dxOld_)
            );
            if (beta == scalar(0))
            {
                Info<< "Computed negative beta. Resetting to zero" << endl;
            }
        }

        scalarField s(dx + beta*sOld_);

        correction_ = Zero;
        forAll(activeDesignVars_, varI)
        {
            correction_[activeDesignVars_[varI]] = eta_*s[varI];
        }

        // Store fields for the next iteration
        dxOld_ = dx;
        sOld_ = s;
    }

    ++counter_;
}


void Foam::conjugateGradient::updateOldCorrection
(
    const scalarField& oldCorrection
)
{
    sOld_.map(oldCorrection, activeDesignVars_);
    sOld_ /= eta_;
    correction_ = oldCorrection;
}


bool Foam::conjugateGradient::writeData(Ostream& os) const
{
    dxOld_.writeEntry("dxOld", os);
    sOld_.writeEntry("sOld", os);

    return updateMethod::writeData(os);
}


// ************************************************************************* //
