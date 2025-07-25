/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "viewFactor.H"
#include "surfaceFields.H"
#include "constants.H"
#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "boundaryRadiationProperties.H"
#include "lduCalculatedProcessorField.H"



using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(viewFactor, 0);
        addToRadiationRunTimeSelectionTables(viewFactor);
    }
}

const Foam::word Foam::radiation::viewFactor::viewFactorWalls
    = "viewFactorWall";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::viewFactor::initialise()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (!finalAgglom_.typeHeaderOk<labelListIOList>())
    {
        finalAgglom_.setSize(patches.size());
        for (label patchi=0;  patchi < patches.size(); patchi++)
        {
            finalAgglom_[patchi] = identity(patches[patchi].size());
        }
    }

    coarseMesh_.reset
    (
        new singleCellFvMesh
        (
            IOobject
            (
                "coarse:" + mesh_.name(),
                mesh_.polyMesh::instance(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            finalAgglom_
        )
    );

    const polyBoundaryMesh& coarsePatches = coarseMesh_->boundaryMesh();

    selectedPatches_ = mesh_.boundaryMesh().indices(viewFactorWalls);

    for (const label patchi : selectedPatches_)
    {
        nLocalCoarseFaces_ += coarsePatches[patchi].size();
    }

    if (debug)
    {
        Pout<< "radiation::viewFactor::initialise() Selected patches:"
            << selectedPatches_ << endl;
        Pout<< "radiation::viewFactor::initialise() Number of coarse faces:"
            << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = returnReduce(nLocalCoarseFaces_, sumOp<label>());

    DebugInFunction
        << "Total number of clusters : " << totalNCoarseFaces_ << endl;

    useDirect_ = coeffs_.getOrDefault<bool>("useDirectSolver", true);



    map_.reset
    (
        new IOmapDistribute
        (
            IOobject
            (
                "mapDist",
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );

    FmyProc_.reset
    (
        new scalarListIOList
        (
            IOobject
            (
                "F",
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );

    globalFaceFaces_.reset
    (
        new labelListIOList
        (
            IOobject
            (
                "globalFaceFaces",
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );


    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Size coarse qr
    qrBandI_.setSize(nBands_);
    forAll (qrBandI_, bandI)
    {
        qrBandI_[bandI].setSize(nLocalCoarseFaces_, 0.0);
    }

    if (!useDirect_)
    {
        DynamicList<label> dfaceJ;

        // Per processor to owner (local)/neighbour (remote)
        List<DynamicList<label>> procOwner(Pstream::nProcs());
        List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());

        forAll (globalFaceFaces_(), iFace)
        {
            const labelList& visFaces = globalFaceFaces_()[iFace];
            forAll (visFaces, j)
            {
                label gFacej = visFaces[j];
                label proci = globalNumbering.whichProcID(gFacej);
                label faceJ = globalNumbering.toLocal(proci, gFacej);

                if (Pstream::myProcNo() == proci)
                {
                    edge e(iFace, faceJ);
                    if (rays_.insert(e))
                    {
                        dfaceJ.append(j);
                    }
                }
                else
                {
                    label gFaceI =
                        globalNumbering.toGlobal(Pstream::myProcNo(), iFace);

                    label proci = globalNumbering.whichProcID(gFacej);

                    label facei =
                        globalNumbering.toLocal(Pstream::myProcNo(), gFaceI);

                    label remoteFacei = globalNumbering.toLocal(proci, gFacej);

                    procOwner[proci].append(facei);
                    dynProcNeighbour[proci].append(remoteFacei);
                }
            }
        }

        mapRayToFmy_.transfer(dfaceJ);

        labelList upper(rays_.size(), -1);
        labelList lower(rays_.size(), -1);

        const edgeList raysLst(rays_.sortedToc());
        label rayI = 0;
        for (const auto& e : raysLst)
        {
            label faceI = e.start();
            label faceJ = e.end();
            upper[rayI] = max(faceI, faceJ);
            lower[rayI] = min(faceI, faceJ);
            rayI++;
        }

        labelListList procNeighbour(dynProcNeighbour.size());
        forAll(procNeighbour, i)
        {
            procNeighbour[i] = std::move(dynProcNeighbour[i]);
        }
        labelListList mySendCells;
        Pstream::exchange<labelList, label>(procNeighbour, mySendCells);

        label nbri = 0;
        forAll(procOwner, proci)
        {
            if (procOwner[proci].size())
            {
                nbri++;
            }
            if (mySendCells[proci].size())
            {
                nbri++;
            }
        }

        DebugInFunction<< "Number of procBound : " <<  nbri << endl;

        PtrList<const lduPrimitiveProcessorInterface> primitiveInterfaces(nbri);
        internalCoeffs_.clear();
        boundaryCoeffs_.resize_null(nbri);

        nbri = 0;

        procToInterface_.setSize(Pstream::nProcs(), -1);

        forAll(procOwner, proci)
        {
            if (proci < Pstream::myProcNo() && procOwner[proci].size())
            {
                if (debug)
                {
                    Pout<< "Adding interface " << nbri
                        << " to receive my " << procOwner[proci].size()
                        << " from " << proci << endl;
                }
                procToInterface_[proci] = nbri;
                primitiveInterfaces.set
                (
                    nbri++,
                    new lduPrimitiveProcessorInterface
                    (
                        procOwner[proci],
                        Pstream::myProcNo(),
                        proci,
                        tensorField(0),
                        Pstream::msgType()+2
                    )
                );
            }
            else if (proci > Pstream::myProcNo() && mySendCells[proci].size())
            {
                if (debug)
                {
                    Pout<< "Adding interface " << nbri
                        << " to send my " << mySendCells[proci].size()
                        << " to " << proci << endl;
                }
                primitiveInterfaces.set
                (
                    nbri++,
                    new lduPrimitiveProcessorInterface
                    (
                        mySendCells[proci],
                        Pstream::myProcNo(),
                        proci,
                        tensorField(0),
                        Pstream::msgType()+2
                    )
                );
            }
        }
        forAll(procOwner, proci)
        {
            if (proci > Pstream::myProcNo() && procOwner[proci].size())
            {
                if (debug)
                {
                    Pout<< "Adding interface " << nbri
                        << " to receive my " << procOwner[proci].size()
                        << " from " << proci << endl;
                }
                procToInterface_[proci] = nbri;
                primitiveInterfaces.set
                (
                    nbri++,
                    new lduPrimitiveProcessorInterface
                    (
                        procOwner[proci],
                        Pstream::myProcNo(),
                        proci,
                        tensorField(0),
                        Pstream::msgType()+3
                    )
                );
            }
            else if (proci < Pstream::myProcNo() && mySendCells[proci].size())
            {
                if (debug)
                {
                    Pout<< "Adding interface " << nbri
                        << " to send my " << mySendCells[proci].size()
                        << " to " << proci << endl;
                }
                primitiveInterfaces.set
                (
                    nbri++,
                    new lduPrimitiveProcessorInterface
                    (
                        mySendCells[proci],
                        Pstream::myProcNo(),
                        proci,
                        tensorField(0),
                        Pstream::msgType()+3
                    )
                );
            }
        }

        forAll (boundaryCoeffs_, proci)
        {
            boundaryCoeffs_.set
            (
                proci,
                new Field<scalar>
                (
                    primitiveInterfaces[proci].faceCells().size(),
                    Zero
                )
            );
        }

        lduInterfacePtrsList allInterfaces(primitiveInterfaces.size());

        forAll(primitiveInterfaces, i)
        {
            const lduPrimitiveProcessorInterface& pp = primitiveInterfaces[i];

            allInterfaces.set(i, &pp);
        }

        const lduSchedule ps
        (
            lduPrimitiveMesh::nonBlockingSchedule
                <lduPrimitiveProcessorInterface>(allInterfaces)
        );

        PtrList<const lduInterface> allInterfacesPtr(allInterfaces.size());
        forAll (allInterfacesPtr, i)
        {
            const lduPrimitiveProcessorInterface& pp = primitiveInterfaces[i];

            allInterfacesPtr.set
            (
                i,
                new lduPrimitiveProcessorInterface(pp)
            );
        }

        lduPtr_.reset
        (
            new lduPrimitiveMesh
            (
                globalFaceFaces_().size(),    //rays_.size(),
                lower,
                upper,
                allInterfacesPtr,
                ps,
                UPstream::worldComm
            )
        );

        // Set size for local lduMatrix
        matrixPtr_.reset(new lduMatrix(lduPtr_()));

        scalarListList& myF = FmyProc_();

        if (coeffs_.get<bool>("smoothing"))
        {
            scalar maxDelta = 0;
            scalar totalDelta = 0;

            if (myF.size())
            {
                forAll (myF, i)
                {
                    scalar sumF = 0.0;
                    scalarList& myFij = myF[i];
                    forAll (myFij, j)
                    {
                        sumF += myFij[j];
                    }
                    const scalar delta = sumF - 1.0;
                    forAll (myFij, j)
                    {
                        myFij[j] *= (1.0 - delta/(sumF + 0.001));
                    }
                    totalDelta += delta;
                    if (delta > maxDelta)
                    {
                        maxDelta = delta;
                    }
                }
                totalDelta /= myF.size();
            }
            reduce(totalDelta, sumOp<scalar>());
            reduce(maxDelta, maxOp<scalar>());
            Info << "Smoothing average delta : " << totalDelta << endl;
            Info << "Smoothing maximum delta : " << maxDelta << nl << endl;
        }
    }

    if (useDirect_)
    {
        List<labelListList> globalFaceFacesProc(Pstream::nProcs());
        globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces_();
        Pstream::gatherList(globalFaceFacesProc);

        List<scalarListList> F(Pstream::nProcs());
        F[Pstream::myProcNo()] = FmyProc_();
        Pstream::gatherList(F);

        if (Pstream::master())
        {
            Fmatrix_.reset
            (
                new scalarSquareMatrix(totalNCoarseFaces_, Zero)
            );

            DebugInFunction
                << "Insert elements in the matrix..." << endl;

            for (const int procI : Pstream::allProcs())
            {
                insertMatrixElements
                (
                    globalNumbering,
                    procI,
                    globalFaceFacesProc[procI],
                    F[procI],
                    Fmatrix_()
                );
            }

            if (coeffs_.get<bool>("smoothing"))
            {
                DebugInFunction << "Smoothing the matrix..." << endl;

                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    scalar sumF = 0.0;
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        sumF += Fmatrix_()(i, j);
                    }

                    const scalar delta = sumF - 1.0;
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        Fmatrix_()(i, j) *= (1.0 - delta/(sumF + 0.001));
                    }
                }
            }

            coeffs_.readEntry("constantEmissivity", constEmissivity_);
            if (constEmissivity_)
            {
                CLU_.reset
                (
                    new scalarSquareMatrix(totalNCoarseFaces_, Zero)
                );

                pivotIndices_.setSize(CLU_().m());
            }
        }
    }

    coeffs_.readIfPresent("useSolarLoad", useSolarLoad_);

    if (useSolarLoad_)
    {
        const dictionary& solarDict = this->subDict("solarLoadCoeffs");
        solarLoad_.reset(new solarLoad(solarDict, T_));

        if (solarLoad_->nBands() != nBands_)
        {
            FatalErrorInFunction
                << "Solar radiation and view factor band numbers "
                << "are different"
                << abort(FatalError);
        }

        Info<< "Creating Solar Load Model " << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::viewFactor::viewFactor(const volScalarField& T)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    map_(),
    coarseMesh_(),
//     (
//         IOobject
//         (
//             "coarse:" + mesh_.name(),
//             mesh_.polyMesh::instance(),
//             mesh_.time(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE,
//             IOobject::NO_REGISTER
//         ),
//         mesh_,
//         finalAgglom_
//     ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0),
    useSolarLoad_(false),
    solarLoad_(),
    nBands_(coeffs_.getOrDefault<label>("nBands", 1)),
    FmyProc_()
{
    initialise();
}


Foam::radiation::viewFactor::viewFactor
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    map_(),
    coarseMesh_(),
//     (
//         IOobject
//         (
//             "coarse:" + mesh_.name(),
//             mesh_.polyMesh::instance(),
//             mesh_.time(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE,
//             IOobject::NO_REGISTER
//         ),
//         mesh_,
//         finalAgglom_
//     ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0),
    useSolarLoad_(false),
    solarLoad_(),
    nBands_(coeffs_.getOrDefault<label>("nBands", 1)),
    FmyProc_()
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::viewFactor::read()
{
    if (radiationModel::read())
    {
        return true;
    }

    return false;
}


void Foam::radiation::viewFactor::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(procI, faceI);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


void Foam::radiation::viewFactor::calculate()
{
    // Store previous iteration
    qr_.storePrevIter();

    if (useSolarLoad_)
    {
        solarLoad_->calculate();
    }

    // Global net radiation
    scalarField qNet(totalNCoarseFaces_, 0.0);

    // Local net radiation
    scalarField qTotalCoarse(nLocalCoarseFaces_, 0.0);

    // Referen to fvMesh qr field
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    globalIndex globalNumbering(nLocalCoarseFaces_);

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(mesh_);

    for (label bandI = 0; bandI < nBands_; bandI++)
    {
        // Global bandI radiation
        scalarField qBandI(totalNCoarseFaces_, 0.0);

        scalarField compactCoarseT4(map_->constructSize(), 0.0);
        scalarField compactCoarseE(map_->constructSize(), 0.0);
        scalarField compactCoarseHo(map_->constructSize(), 0.0);

        // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
        DynamicList<scalar> localCoarseT4ave(nLocalCoarseFaces_);
        DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
        DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

        forAll(selectedPatches_, i)
        {
            label patchID = selectedPatches_[i];

            const scalarField& Tp = T_.boundaryField()[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

            fvPatchScalarField& qrPatch = qrBf[patchID];

            auto& qrp =
                refCast
                <
                    greyDiffusiveViewFactorFixedValueFvPatchScalarField
                >(qrPatch);

            const tmp<scalarField> teb =
                boundaryRadiation.emissivity(patchID, bandI, nullptr, &Tp);

            const scalarField& eb = teb();

            const tmp<scalarField> tHoi = qrp.qro(bandI);
            const scalarField& Hoi = tHoi();

            const polyPatch& pp = coarseMesh_->boundaryMesh()[patchID];

            const labelList& coarsePatchFace =
                coarseMesh_->patchFaceMap()[patchID];

            scalarList T4ave(pp.size(), 0.0);
            scalarList Eave(pp.size(), 0.0);
            scalarList Hoiave(pp.size(), 0.0);

            if (pp.size() > 0)
            {
                const labelList& agglom = finalAgglom_[patchID];
                label nAgglom = max(agglom) + 1;

                labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

                forAll(coarseToFine, coarseI)
                {
                    const label coarseFaceID = coarsePatchFace[coarseI];
                    const labelList& fineFaces = coarseToFine[coarseFaceID];
                    UIndirectList<scalar> fineSf
                    (
                        sf,
                        fineFaces
                    );

                    const scalar area = sum(fineSf());

                    // Temperature, emissivity and external flux area weighting
                    forAll(fineFaces, j)
                    {
                        label facei = fineFaces[j];
                        T4ave[coarseI] += (pow4(Tp[facei])*sf[facei])/area;
                        Eave[coarseI] += (eb[facei]*sf[facei])/area;
                        Hoiave[coarseI] += (Hoi[facei]*sf[facei])/area;
                    }
                }
            }

            localCoarseT4ave.append(T4ave);
            localCoarseEave.append(Eave);
            localCoarseHoave.append(Hoiave);

        }

        // Fill the local values to distribute
        SubList<scalar>(compactCoarseT4, nLocalCoarseFaces_) =
            localCoarseT4ave;
        SubList<scalar>(compactCoarseE, nLocalCoarseFaces_) = localCoarseEave;
        SubList<scalar>(compactCoarseHo, nLocalCoarseFaces_) =
            localCoarseHoave;


        // Distribute data
        map_->distribute(compactCoarseT4);
        map_->distribute(compactCoarseE);
        map_->distribute(compactCoarseHo);

        // Distribute local global ID
        labelList compactGlobalIds(map_->constructSize(), Zero);

        SubList<label>
        (
            compactGlobalIds,
            nLocalCoarseFaces_
        ) = identity
            (
                globalNumbering.localSize(),
                globalNumbering.localStart()
            );

        map_->distribute(compactGlobalIds);

        if (!useDirect_)
        {
            const labelList globalToCompact
            (
                invert(totalNCoarseFaces_, compactGlobalIds)
            );

            scalarField& diag = matrixPtr_->diag(localCoarseEave.size());
            scalarField& upper = matrixPtr_->upper(rays_.size());
            scalarField& lower = matrixPtr_->lower(rays_.size());

            scalarField source(nLocalCoarseFaces_, 0);


            // Local diag and source
            forAll(source, i)
            {
                const scalar sigmaT4 =
                    physicoChemical::sigma.value()*localCoarseT4ave[i];

                diag[i] = 1/localCoarseEave[i];

                source[i] += -sigmaT4 + localCoarseHoave[i];
            }

            // Local matrix coefficients
            if (!constEmissivity_ || iterCounter_ == 0)
            {
                const edgeList raysLst(rays_.sortedToc());

                label rayI = 0;
                for (const auto& e : raysLst)
                {
                    label facelJ = e.end();
                    label faceI = e.start();

                    label faceJ = mapRayToFmy_[rayI];

                    lower[rayI] =
                        (1-1/localCoarseEave[faceI])*FmyProc_()[faceI][faceJ];

                    upper[rayI] =
                        (1-1/localCoarseEave[facelJ])*FmyProc_()[faceI][faceJ];

                    rayI++;
                }
                iterCounter_++;
            }

            // Extra local contribution to the source and extra matrix coefs
            // from other procs (boundaryCoeffs)
            label nInterfaces = lduPtr_().interfaces().size();
            labelList boundCoeffI(nInterfaces, Zero);
            forAll (globalFaceFaces_(), iFace)
            {
                const labelList& visFaces = globalFaceFaces_()[iFace];
                forAll (visFaces, jFace)
                {
                    label gFacej = visFaces[jFace];
                    label proci = globalNumbering.whichProcID(gFacej);
                    if (Pstream::myProcNo() == proci)
                    {
                        label lFacej =
                            globalNumbering.toLocal
                            (
                                Pstream::myProcNo(),
                                gFacej
                            );

                        const scalar sigmaT4 =
                            physicoChemical::sigma.value()
                            *localCoarseT4ave[lFacej];

                        source[iFace] += FmyProc_()[iFace][jFace]*sigmaT4;
                    }
                    else
                    {
                        label compactFaceJ = globalToCompact[gFacej];
                        const scalar sigmaT4 =
                            physicoChemical::sigma.value()
                            *compactCoarseT4[compactFaceJ];

                        source[iFace] += FmyProc_()[iFace][jFace]*sigmaT4;

                        label interfaceI = procToInterface_[proci];

                        boundaryCoeffs_
                            [interfaceI][boundCoeffI[interfaceI]++] =
                                -(1-1/compactCoarseE[compactFaceJ])
                                *FmyProc_()[iFace][jFace];
                    }
                }
            }

            PtrList<const lduInterfaceField> interfaces(nInterfaces);
            for(label i = 0; i < interfaces.size(); i++)
            {
                interfaces.set
                (
                    i,
                    new lduCalculatedProcessorField<scalar>
                    (
                        lduPtr_().interfaces()[i]
                    )
                );
            }

            const dictionary& solverControls =
                qr_.mesh().solverDict
                (
                    qr_.select(qr_.mesh().data().isFinalIteration())
                );

            // Solver call
            solverPerformance solverPerf = lduMatrix::solver::New
            (
                "qr",
                matrixPtr_(),
                boundaryCoeffs_,
                internalCoeffs_,
                interfaces,
                solverControls
            )->solve(qrBandI_[bandI], source);

            solverPerf.print(Info.masterStream(qr_.mesh().comm()));

            qTotalCoarse += qrBandI_[bandI];
        }

        if (useDirect_)
        {
            // Create global size vectors
            scalarField T4(totalNCoarseFaces_, 0.0);
            scalarField E(totalNCoarseFaces_, 0.0);
            scalarField qrExt(totalNCoarseFaces_, 0.0);

            // Fill lists from compact to global indexes.
            forAll(compactCoarseT4, i)
            {
                T4[compactGlobalIds[i]] = compactCoarseT4[i];
                E[compactGlobalIds[i]] = compactCoarseE[i];
                qrExt[compactGlobalIds[i]] = compactCoarseHo[i];
            }

            Pstream::listCombineReduce(T4, maxEqOp<scalar>());
            Pstream::listCombineReduce(E, maxEqOp<scalar>());
            Pstream::listCombineReduce(qrExt, maxEqOp<scalar>());

            if (Pstream::master())
            {
                // Variable emissivity
                if (!constEmissivity_)
                {
                    scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

                    for (label i=0; i<totalNCoarseFaces_; i++)
                    {
                        for (label j=0; j<totalNCoarseFaces_; j++)
                        {
                            const scalar invEj = 1.0/E[j];
                            const scalar sigmaT4 =
                                physicoChemical::sigma.value()*T4[j];

                            if (i==j)
                            {
                                C(i, j) = invEj;
                                qBandI[i] += -sigmaT4 +  qrExt[j];
                            }
                            else
                            {
                                C(i, j) = (1.0 - invEj)*Fmatrix_()(i, j);
                                qBandI[i] += Fmatrix_()(i, j)*sigmaT4;
                            }
                        }
                    }

                    Info<< "Solving view factor equations for band :"
                        << bandI << endl;

                    // Negative coming into the fluid
                    LUsolve(C, qBandI);
                }
                else //Constant emissivity
                {
                    // Initial iter calculates CLU and caches it
                    if (iterCounter_ == 0)
                    {
                        for (label i=0; i<totalNCoarseFaces_; i++)
                        {
                            for (label j=0; j<totalNCoarseFaces_; j++)
                            {
                                const scalar invEj = 1.0/E[j];
                                if (i==j)
                                {
                                    CLU_()(i, j) = invEj;
                                }
                                else
                                {
                                    CLU_()(i, j) =
                                        (1.0-invEj)*Fmatrix_()(i, j);
                                }
                            }
                        }

                        DebugInFunction << "\nDecomposing C matrix..." << endl;

                        LUDecompose(CLU_(), pivotIndices_);
                    }

                    for (label i=0; i<totalNCoarseFaces_; i++)
                    {
                        for (label j=0; j<totalNCoarseFaces_; j++)
                        {
                            const scalar sigmaT4 =
                                constant::physicoChemical::sigma.value()*T4[j];

                            if (i==j)
                            {
                                qBandI[i] += -sigmaT4 + qrExt[j];
                            }
                            else
                            {
                                qBandI[i] += Fmatrix_()(i, j)*sigmaT4;
                            }
                        }
                    }

                    Info<< "Solving view factor equations for band : "
                        << bandI  << endl;

                    LUBacksubstitute(CLU_(), pivotIndices_, qBandI);
                    iterCounter_ ++;
                }
            }
            // Broadcast qBandI and fill qr
            Pstream::broadcast(qBandI);

            qNet += qBandI;
        }
    }

    label globCoarseId = 0;
    for (const label patchID : selectedPatches_)
    {
        const polyPatch& pp = coarseMesh_->boundaryMesh()[patchID];

        if (pp.size() > 0)
        {
            scalarField& qrp = qrBf[patchID];
            //const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_->patchFaceMap()[patchID];

            //scalar heatFlux = 0.0;
            forAll(coarseToFine, coarseI)
            {
                label globalCoarse = globalNumbering.toGlobal
                (
                    Pstream::myProcNo(),
                    globCoarseId
                );

                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];

                for (const label facei : fineFaces)
                {
                    if (useDirect_)
                    {
                        qrp[facei] = qNet[globalCoarse];
                    }
                    else
                    {
                        qrp[facei] = qTotalCoarse[globCoarseId];
                    }
                    //heatFlux += qrp[facei]*sf[facei];
                }
                globCoarseId++;
            }
        }
    }

    if (debug)
    {
        for (const label patchID : selectedPatches_)
        {
            const scalarField& qrp = qrBf[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qrp*magSf);

            Info<< "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qr if necessary
    qr_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::viewFactor::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar
        (
            dimMass/dimLength/pow3(dimTime)/pow4(dimTemperature), Zero
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::viewFactor::Ru() const
{
    return DimensionedField<scalar, volMesh>::New
    (
        "Ru",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    );
}


// ************************************************************************* //
