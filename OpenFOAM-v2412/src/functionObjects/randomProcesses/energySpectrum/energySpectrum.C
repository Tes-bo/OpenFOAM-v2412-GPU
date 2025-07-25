/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "energySpectrum.H"
#include "fft.H"
#include "fvMesh.H"
#include "boundBox.H"
#include "OFstream.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(energySpectrum, 0);
    addToRunTimeSelectionTable(functionObject, energySpectrum, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::energySpectrum::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Turbulence energy spectra");

    writeCommented(os, "kappa E(kappa)");

    os  << endl;
}


void Foam::functionObjects::energySpectrum::calcAndWriteSpectrum
(
    const vectorField& U,
    const vectorField& C,
    const vector& c0,
    const vector& deltaC,
    const Vector<int>& N,
    const scalar kappaNorm
)
{
    Log << "Computing FFT" << endl;
    const complexVectorField Uf
    (
        fft::forwardTransform
        (
            ComplexField(U, vector::zero),
            List<int>({N.x(), N.y(), N.z()})
        )
       /scalar(cmptProduct(N))
    );


    Log << "Computing wave numbers" << endl;
    const label nMax = cmptMax(N);
    scalarField kappa(nMax);
    forAll(kappa, kappai)
    {
        kappa[kappai] = kappai*kappaNorm;
    }

    Log << "Computing energy spectrum" << endl;
    scalarField E(nMax, Zero);
    const scalarField Ec(0.5*magSqr(Uf));
    forAll(C, celli)
    {
        point Cc(C[celli] - c0);
        label i = round((Cc.x() - 0.5*deltaC.x())/(deltaC.x())*(N.x() - 1));
        label j = round((Cc.y() - 0.5*deltaC.y())/(deltaC.y())*(N.y() - 1));
        label k = round((Cc.z() - 0.5*deltaC.z())/(deltaC.z())*(N.z() - 1));
        label kappai = round(Foam::sqrt(scalar(i*i + j*j + k*k)));

        E[kappai] += Ec[celli];
    }

    E /= kappaNorm;

    Log << "Writing spectrum" << endl;
    autoPtr<OFstream> osPtr = newFileAtTime(name(), time_.value());
    OFstream& os = osPtr.ref();
    writeFileHeader(os);

    forAll(kappa, kappai)
    {
        os  << kappa[kappai] << tab << E[kappai] << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::energySpectrum::energySpectrum
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    cellAddr_(mesh_.nCells()),
    UName_("U"),
    N_(Zero),
    c0_(Zero),
    deltaC_(Zero),
    kappaNorm_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::energySpectrum::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readIfPresent("U", UName_);

    const boundBox meshBb(mesh_.bounds());

    // Assume all cells are the same size...
    boundBox cellBb(mesh_.cellBb(0));

    const vector L(meshBb.span());
    const vector nCellXYZ(cmptDivide(L, cellBb.span()));

    N_ = Vector<int>
    (
        round(nCellXYZ.x()),
        round(nCellXYZ.z()),
        round(nCellXYZ.z())
    );

    // Check that the mesh is a structured box
    vector cellDx(cellBb.span());
    vector expectedMax(N_.x()*cellDx.x(), N_.y()*cellDx.y(), N_.z()*cellDx.z());
    vector relativeSize(cmptDivide(L, expectedMax));
    for (direction i = 0; i < 3; ++i)
    {
        if (mag(relativeSize[i] - 1) > 1e-3)
        {
            FatalErrorInFunction
                << name() << " function object is only appropriate for "
                << "isotropic structured IJK meshes. Mesh extents: " << L
                << ", computed IJK mesh extents: " << expectedMax
                << exit(FatalError);
        }
    }
    Log << "Mesh extents (deltax,deltay,deltaz): " << L << endl;
    Log << "Number of cells (Nx,Ny,Nz): " << N_ << endl;

    // Map into i-j-k co-ordinates
    const vectorField& C = mesh_.C();
    c0_ = returnReduce(min(C), minOp<vector>());
    deltaC_ = returnReduce(max(C), maxOp<vector>());
    deltaC_ -= c0_;

    forAll(C, celli)
    {
        label i = round((C[celli].x() - c0_.x())/(deltaC_.x())*(N_.x() - 1));
        label j = round((C[celli].y() - c0_.y())/(deltaC_.y())*(N_.y() - 1));
        label k = round((C[celli].z() - c0_.z())/(deltaC_.z())*(N_.z() - 1));

        cellAddr_[celli] = k + j*N_.y() + i*N_.y()*N_.z();
    }

    kappaNorm_ = constant::mathematical::twoPi/cmptMax(L);

    return true;
}


bool Foam::functionObjects::energySpectrum::execute()
{
    return true;
}


bool Foam::functionObjects::energySpectrum::write()
{
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const vectorField& Uc = U.primitiveField();
    const vectorField& C = mesh_.C();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs;

        {
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster << Uc << C << cellAddr_;
        }

        pBufs.finishedGathers();

        if (Pstream::master())
        {
            const label nGlobalCells(cmptProduct(N_));
            vectorField Uijk(nGlobalCells);
            vectorField Cijk(nGlobalCells);

            for (const int proci : Pstream::allProcs())
            {
                UIPstream fromProc(proci, pBufs);
                vectorField Up;
                vectorField Cp;
                labelList cellAddrp;
                fromProc >> Up >> Cp >> cellAddrp;
                UIndirectList<vector>(Uijk, cellAddrp) = Up;
                UIndirectList<vector>(Cijk, cellAddrp) = Cp;
            }

            calcAndWriteSpectrum(Uijk, Cijk, c0_, deltaC_, N_, kappaNorm_);
        }

    }
    else
    {
        vectorField Uijk(Uc, cellAddr_);
        vectorField Cijk(C, cellAddr_);

        calcAndWriteSpectrum(Uijk, Cijk, c0_, deltaC_, N_, kappaNorm_);
    }

    return true;
}


// ************************************************************************* //
