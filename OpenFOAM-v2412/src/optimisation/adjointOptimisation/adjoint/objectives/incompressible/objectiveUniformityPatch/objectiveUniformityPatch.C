/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "objectiveUniformityPatch.H"
#include "createZeroField.H"
#include "coupledFvPatch.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveUniformityPatch, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveUniformityPatch,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void objectiveUniformityPatch::initialize()
{
    // If patches are prescribed, use them

    wordRes patchSelection;
    if (dict().readIfPresent("patches", patchSelection))
    {
        patches_ = mesh_.boundaryMesh().patchSet(patchSelection).sortedToc();
    }
    // Otherwise, pick them up based on the mass flow.
    // Note: a non-zero U initialisation should be used in order to pick up the
    // outlet patches correctly
    else
    {
        WarningInFunction
            << "No patches provided to " << type() << ". "
            << "Choosing them according to the patch mass flows" << nl;

        DynamicList<label> objectiveReportPatches(mesh_.boundary().size());
        const surfaceScalarField& phi = vars_.phiInst();
        forAll(mesh_.boundary(), patchI)
        {
           const fvsPatchScalarField& phiPatch = phi.boundaryField()[patchI];
           if (!isA<coupledFvPatch>(mesh_.boundary()[patchI]))
           {
                const scalar mass = gSum(phiPatch);
                if (mass > SMALL)
                {
                    objectiveReportPatches.push_back(patchI);
                }
            }
        }
        patches_.transfer(objectiveReportPatches);
    }
    UMean_.setSize(patches_.size(), Zero);
    UVar_.setSize(patches_.size(), Zero);

    if (patches_.empty())
    {
        FatalErrorInFunction
            << "No valid patch name on which to minimize " << type() << endl
            << exit(FatalError);
    }
    if (debug)
    {
        Info<< "Minimizing " << type() << " in patches:" << endl;
        forAll(patches_, pI)
        {
            Info<< "\t " << mesh_.boundary()[patches_[pI]].name() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveUniformityPatch::objectiveUniformityPatch
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    patches_(),
    UMean_(),
    UVar_()
{
    // Find inlet/outlet patches
    initialize();

    // Allocate boundary field pointers
    bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    bdJdvtPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveUniformityPatch::J()
{
    J_ = Zero;

    const volVectorField& U = vars_.UInst();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        const scalar sumMagSf = gSum(magSf);
        const fvPatchVectorField& Ub = U.boundaryField()[patchI];
        UMean_[oI] = gSum(Ub*magSf)/sumMagSf;
        UVar_[oI] = gSum(magSqr(Ub - UMean_[oI])*magSf)/sumMagSf;
        J_ += 0.5*UVar_[oI];
    }

    return J_;
}


void objectiveUniformityPatch::update_boundarydJdv()
{
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        const scalar sumMagSf = gSum(magSf);
        const fvPatchVectorField& Ub = U.boundaryField()[patchI];
        bdJdvPtr_()[patchI] = (Ub - UMean_[oI])/sumMagSf;
    }
}


void objectiveUniformityPatch::update_boundarydJdvn()
{
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        const scalar sumMagSf = gSum(magSf);
        const fvPatchVectorField& Ub = U.boundaryField()[patchI];
        tmp<vectorField> nf = mesh_.boundary()[patchI].nf();
        bdJdvnPtr_()[patchI] = ((Ub - UMean_[oI]) & nf)/sumMagSf;
    }
}


void objectiveUniformityPatch::update_boundarydJdvt()
{
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        const scalar sumMagSf = gSum(magSf);
        const fvPatchVectorField& Ub = U.boundaryField()[patchI];
        tmp<vectorField> nf = mesh_.boundary()[patchI].nf();
        vectorField UdiffTangent(Ub - UMean_[oI]);

        bdJdvtPtr_()[patchI] =
            (UdiffTangent - (UdiffTangent & nf())*nf())/sumMagSf;
    }
}


void objectiveUniformityPatch::addHeaderColumns() const
{
    OFstream& file = objFunctionFilePtr_();
    for (const label patchI : patches_)
    {
        const word patchName(mesh_.boundary()[patchI].name());
        file<< setw(width_) << word(patchName + "-" + "UMean") <<  " ";
        file<< setw(width_) << word(patchName + "-" + "UVar") <<  " ";
        file<< setw(width_) << word(patchName + "-" + "UStd") <<  " ";
    }
}


void objectiveUniformityPatch::addColumnValues() const
{
    OFstream& file = objFunctionFilePtr_();
    forAll(UMean_, pI)
    {
        file<< setw(width_) << mag(UMean_[pI]) << " ";
        file<< setw(width_) << UVar_[pI] << " ";
        file<< setw(width_) << sqrt(UVar_[pI]) << " ";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
