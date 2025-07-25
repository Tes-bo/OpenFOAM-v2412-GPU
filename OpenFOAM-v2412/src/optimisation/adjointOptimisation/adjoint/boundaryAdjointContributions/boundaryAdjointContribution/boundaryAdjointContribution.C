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

#include "boundaryAdjointContribution.H"
#include "fvPatchFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryAdjointContribution, 0);
defineRunTimeSelectionTable(boundaryAdjointContribution, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryAdjointContribution::boundaryAdjointContribution
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
:
    patch_(patch)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<boundaryAdjointContribution> boundaryAdjointContribution::New
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
{
    auto* ctorPtr = dictionaryConstructorTable(simulationType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "boundaryAdjointContribution",
            simulationType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return
        autoPtr<boundaryAdjointContribution>
        (
            ctorPtr
            (
                managerName,
                adjointSolverName,
                simulationType,
                patch
            )
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> boundaryAdjointContribution::adjointTMVariable1Source()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::adjointTMVariable2Source()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::dJdnut()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<tensorField> boundaryAdjointContribution::dJdGradU()
{
    return tmp<tensorField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable1Diffusion()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable2Diffusion()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable1()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable2()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<fvPatchScalarField>
boundaryAdjointContribution::turbulentDiffusivity() const
{
    NotImplemented;
    return nullptr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
