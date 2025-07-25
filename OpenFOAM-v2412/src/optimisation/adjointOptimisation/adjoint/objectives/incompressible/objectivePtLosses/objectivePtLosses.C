/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "objectivePtLosses.H"
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

defineTypeNameAndDebug(objectivePtLosses, 1);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectivePtLosses,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void objectivePtLosses::initialize()
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
            << "No patches provided to PtLosses. "
            << "Choosing them according to the patch mass flows" << nl;

        DynamicList<label> objectiveReportPatches(mesh_.boundary().size());
        const surfaceScalarField& phi = vars_.phiInst();
        forAll(mesh_.boundary(), patchI)
        {
           const fvsPatchScalarField& phiPatch = phi.boundaryField()[patchI];
           if (!isA<coupledFvPatch>(mesh_.boundary()[patchI]))
           {
                const scalar mass = gSum(phiPatch);
                if (mag(mass) > SMALL)
                {
                    objectiveReportPatches.push_back(patchI);
                }
            }
        }
        patches_.transfer(objectiveReportPatches);
    }
    patchPt_.setSize(patches_.size());

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

objectivePtLosses::objectivePtLosses
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    patches_(0),
    patchPt_(0)
{
    // Find inlet/outlet patches
    initialize();

    // Allocate boundary field pointers
    bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    bdJdvtPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectivePtLosses::J()
{
    J_ = Zero;

    // References
    const volScalarField& p = vars_.pInst();
    const volVectorField& U = vars_.UInst();

    // Inlet/outlet patches
    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];
        const vectorField& Sf = mesh_.boundary()[patchI].Sf();
        scalar pt = -gSum
        (
            (U.boundaryField()[patchI] & Sf)
           *(
                p.boundaryField()[patchI]
              + 0.5*magSqr(U.boundaryField()[patchI])
            )
        );
        patchPt_[oI] = mag(pt);
        J_ += pt;
    }

    return J_;
}


void objectivePtLosses::update_boundarydJdp()
{
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];

        tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
        const vectorField& nf = tnf();

        bdJdpPtr_()[patchI] = -(U.boundaryField()[patchI] & nf)*nf;
    }
}


void objectivePtLosses::update_boundarydJdv()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];

        tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
        const vectorField& nf = tnf();
        const fvPatchVectorField& Ub = U.boundaryField()[patchI];

        bdJdvPtr_()[patchI] =
          - (p.boundaryField()[patchI] + 0.5*magSqr(Ub))*nf
          - (Ub & nf)*Ub;
    }
}


void objectivePtLosses::update_boundarydJdvn()
{
    const volScalarField& p = vars_.p();
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];

        tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
        const vectorField& nf = tnf();

        bdJdvnPtr_()[patchI] =
          - p.boundaryField()[patchI]
          - 0.5*magSqr(U.boundaryField()[patchI])
          - sqr(U.boundaryField()[patchI] & nf);
    }
}


void objectivePtLosses::update_boundarydJdvt()
{
    const volVectorField& U = vars_.U();

    forAll(patches_, oI)
    {
        const label patchI = patches_[oI];

        tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
        const vectorField& nf = tnf();
        scalarField Un(U.boundaryField()[patchI] & nf);

        bdJdvtPtr_()[patchI] = -Un*(U.boundaryField()[patchI] - Un*nf);
    }
}


void objectivePtLosses::addHeaderColumns() const
{
    for (const label patchI : patches_)
    {
        objFunctionFilePtr_()
            << setw(width_) << mesh_.boundary()[patchI].name() <<  " ";
    }
}


void objectivePtLosses::addColumnValues() const
{
    for (const scalar pt : patchPt_)
    {
        objFunctionFilePtr_()
            << setw(width_) << pt << " ";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
