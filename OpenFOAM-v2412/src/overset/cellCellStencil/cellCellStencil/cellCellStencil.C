/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "cellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "oversetFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCellStencil, 0);
    defineRunTimeSelectionTable(cellCellStencil, mesh);
}

const Foam::Enum
<
    Foam::cellCellStencil::cellType
>
Foam::cellCellStencil::cellTypeNames_
({
    { cellType::CALCULATED, "calculated" },
    { cellType::INTERPOLATED, "interpolated" },
    { cellType::HOLE, "hole" },
    { cellType::SPECIAL, "special" },
    { cellType::POROUS, "porous" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencil::cellCellStencil(const fvMesh& mesh)
:
    mesh_(mesh),
    nonInterpolatedFields_({"zoneID"})
{}


Foam::autoPtr<Foam::cellCellStencil> Foam::cellCellStencil::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool update
)
{
    DebugInFunction << "Constructing cellCellStencil" << endl;

    const word stencilType(dict.get<word>("method"));

    auto* ctorPtr = meshConstructorTable(stencilType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "cellCellStencil",
            stencilType,
            *meshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<cellCellStencil>(ctorPtr(mesh, dict, update));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencil::~cellCellStencil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::cellCellStencil::baseName(const word& name)
{
    if (name.ends_with("_0"))
    {
        return baseName(name.substr(0, name.size() - 2));
    }

    return name;
}


void Foam::cellCellStencil::suppressMotionFields()
{
    // Protect local fields from interpolation
    nonInterpolatedFields_.insert("cellInterpolationWeight");
    nonInterpolatedFields_.insert("cellTypes");
    nonInterpolatedFields_.insert("maxMagWeight");

    // For convenience also suppress frequently used displacement field
    {
        nonInterpolatedFields_.insert("cellDisplacement");
        nonInterpolatedFields_.insert("grad(cellDisplacement)");
        const word w("snGradCorr(cellDisplacement)");
        const word d("((viscosity*faceDiffusivity)*magSf)");
        nonInterpolatedFields_.insert("surfaceIntegrate(("+d+"*"+w+"))");
    }

    // For convenience also suppress frequently used displacement field
    {
        nonInterpolatedFields_.insert("cellMotionU");
        nonInterpolatedFields_.insert("grad(cellMotionU)");
        const word w("snGradCorr(cellMotionU)");
        const word d("((viscosity*faceDiffusivity)*magSf)");
        nonInterpolatedFields_.insert("surfaceIntegrate(("+d+"*"+w+"))");
    }
}


const Foam::labelIOList& Foam::cellCellStencil::zoneID(const fvMesh& mesh)
{
    labelIOList* zoneIDPtr = mesh.getObjectPtr<labelIOList>("zoneID");

    if (!zoneIDPtr)
    {
        zoneIDPtr = new labelIOList
        (
            IOobject
            (
                "zoneID",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh.nCells()
        );

        zoneIDPtr->store();

        auto& zoneID = *zoneIDPtr;

        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                mesh.time().findInstance(mesh.dbDir(), "zoneID"),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh
        );
        forAll(volZoneID, celli)
        {
            zoneID[celli] = label(volZoneID[celli]);
        }
    }

    return *zoneIDPtr;
}


Foam::labelList Foam::cellCellStencil::count
(
    const label size,
    const labelUList& lst
)
{
    labelList count(size, Zero);
    forAll(lst, i)
    {
        count[lst[i]]++;
    }
    Pstream::listCombineGather(count, plusEqOp<label>());
    return count;
}


const Foam::wordHashSet& Foam::cellCellStencil::nonInterpolatedFields() const
{
    return nonInterpolatedFields_;
}


Foam::wordHashSet& Foam::cellCellStencil::nonInterpolatedFields()
{
    return nonInterpolatedFields_;
}


bool Foam::cellCellStencil::localStencil(const labelUList& slots) const
{
    forAll(slots, i)
    {
        if (slots[i] >= mesh_.nCells())
        {
            return false;
        }
    }
    return true;
}


void Foam::cellCellStencil::globalCellCells
(
    const globalIndex& gi,
    const polyMesh& mesh,
    const boolList& isValidCell,
    const labelList& selectedCells,
    labelListList& cellCells,
    pointListList& cellCellCentres
)
{
    // For selected cells determine the face neighbours (in global numbering)

    const pointField& cellCentres = mesh.cellCentres();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const cellList& cells = mesh.cells();


    // 1. Determine global cell number on other side of coupled patches

    labelList globalCellIDs(Foam::identity(gi.range()));

    labelList nbrGlobalCellIDs;
    syncTools::swapBoundaryCellList
    (
        mesh,
        globalCellIDs,
        nbrGlobalCellIDs
    );
    pointField nbrCellCentres;
    syncTools::swapBoundaryCellList
    (
        mesh,
        cellCentres,
        nbrCellCentres
    );

    boolList nbrIsValidCell;
    syncTools::swapBoundaryCellList
    (
        mesh,
        isValidCell,
        nbrIsValidCell
    );


    // 2. Collect cell and all its neighbours

    cellCells.setSize(mesh.nCells());
    cellCellCentres.setSize(cellCells.size());

    forAll(selectedCells, i)
    {
        label celli = selectedCells[i];

        const cell& cFaces = cells[celli];
        labelList& stencil = cellCells[celli];
        pointList& stencilPoints = cellCellCentres[celli];
        stencil.setSize(cFaces.size()+1);
        stencilPoints.setSize(stencil.size());
        label compacti = 0;

        // First entry is cell itself
        if (isValidCell[celli])
        {
            stencil[compacti] = globalCellIDs[celli];
            stencilPoints[compacti++] = cellCentres[celli];
        }

        // Other entries are cell neighbours
        forAll(cFaces, i)
        {
            label facei = cFaces[i];
            label bFacei = facei-mesh.nInternalFaces();
            label own = faceOwner[facei];
            label nbrCelli;
            point nbrCc;
            bool isValid = false;
            if (bFacei >= 0)
            {
                nbrCelli = nbrGlobalCellIDs[bFacei];
                nbrCc = nbrCellCentres[bFacei];
                isValid = nbrIsValidCell[bFacei];
            }
            else
            {
                if (own != celli)
                {
                    nbrCelli = gi.toGlobal(own);
                    nbrCc = cellCentres[own];
                    isValid = isValidCell[own];
                }
                else
                {
                    label nei = faceNeighbour[facei];
                    nbrCelli = gi.toGlobal(nei);
                    nbrCc = cellCentres[nei];
                    isValid = isValidCell[nei];
                }
            }

            if (isValid)
            {
                SubList<label> current(stencil, compacti);
                if (!current.found(nbrCelli))
                {
                    stencil[compacti] = nbrCelli;
                    stencilPoints[compacti++] = nbrCc;
                }
            }
        }
        stencil.setSize(compacti);
        stencilPoints.setSize(compacti);
    }
}

void Foam::cellCellStencil::seedCell
(
    const label celli,
    const scalar wantedFraction,
    bitSet& isFront,
    scalarField& fraction
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    forAll(cFaces, i)
    {
        const label nbrFacei = cFaces[i];
        if (fraction[nbrFacei] < wantedFraction)
        {
            fraction[nbrFacei] = wantedFraction;
            isFront.set(nbrFacei);
        }
    }
}


void Foam::cellCellStencil::setUpFront
(
    const labelList& allCellTypes,
    bitSet& isFront
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        const label ownType = allCellTypes[own[facei]];
        const label neiType = allCellTypes[nei[facei]];
        if
        (
            (ownType == HOLE && neiType != HOLE)
         || (ownType != HOLE && neiType == HOLE)
        )
        {
            isFront.set(facei);
        }
    }

    labelList nbrCellTypes;
    syncTools::swapBoundaryCellList(mesh_, allCellTypes, nbrCellTypes);

    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        const label ownType = allCellTypes[own[facei]];
        const label neiType = nbrCellTypes[facei-mesh_.nInternalFaces()];

        if
        (
            (ownType == HOLE && neiType != HOLE)
         || (ownType != HOLE && neiType == HOLE)
        )
        {
            isFront.set(facei);
        }
    }
}


void Foam::cellCellStencil::setUpFrontOnOversetPatch
(
    const labelList& allCellTypes,
    bitSet& isFront
) const
{
    const fvBoundaryMesh& fvm = mesh_.boundary();
    // 'overset' patches
    forAll(fvm, patchi)
    {
        if (isA<oversetFvPatch>(fvm[patchi]))
        {
            const labelList& fc = fvm[patchi].faceCells();
            forAll(fc, i)
            {
                const label celli = fc[i];
                if (allCellTypes[celli] == INTERPOLATED)
                {
                    // Note that acceptors might have been marked hole if
                    // there are no donors in which case we do not want to
                    // walk this out. This is an extreme situation.
                    isFront.set(fvm[patchi].start()+i);
                }
            }
        }
    }
}


void Foam::cellCellStencil::walkFront
(
    const globalIndex& globalCells,
    const scalar layerRelax,
    const labelListList& allStencil,
    labelList& allCellTypes,
    scalarField& allWeight,
    const scalarList& compactCellVol,
    const labelListList& compactStencil,
    const labelList& zoneID,
    const label holeLayers,
    const label useLayer
) const
{
    if (useLayer > holeLayers)
    {
        FatalErrorInFunction<< "useLayer: " << useLayer
            << " is larger than : " <<  holeLayers
            << abort(FatalError);
    }

    if (useLayer == 0)
    {
        FatalErrorInFunction<< "useLayer: " << useLayer
            << " can not be zero."
            << abort(FatalError);
    }

    // Current front
    bitSet isFront(mesh_.nFaces());

    // List of cellTYpes
    DynamicList<labelList> dataSet;
    // List of cellTypes for increasing layers
    DynamicList<scalarField> dAlllWeight;
    // List of average volumen ration interpolated/donor
    DynamicList<scalar> daverageVolRatio;

    // Counting holes of different layers
    DynamicList<label> dHoles;

    const scalarField& V = mesh_.V();

    // Current layer
    labelField nLayer(allCellTypes.size(), 0);

    for (label currLayer = 1; currLayer < holeLayers+1; ++currLayer)
    {
        // Set up original front
        isFront = false;

        setUpFront(allCellTypes, isFront);

        labelList allCellTypesWork(allCellTypes);

        bitSet isFrontWork(isFront);
        label nCurrLayer = currLayer;

        while (nCurrLayer > 1 && returnReduce(isFrontWork.any(), orOp<bool>()))
        {
            bitSet newIsFront(mesh_.nFaces());
            forAll(isFrontWork, facei)
            {
                if (isFrontWork.test(facei))
                {
                    const label own = mesh_.faceOwner()[facei];

                    if (allCellTypesWork[own] != HOLE)
                    {
                        allCellTypesWork[own] = HOLE;
                        newIsFront.set(mesh_.cells()[own]);
                    }

                    if (mesh_.isInternalFace(facei))
                    {
                        const label nei = mesh_.faceNeighbour()[facei];
                        if (allCellTypesWork[nei] != HOLE)
                        {
                            allCellTypesWork[nei] = HOLE;
                            newIsFront.set(mesh_.cells()[nei]);
                        }
                    }
                }
            }
            syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());

            isFrontWork.transfer(newIsFront);

            nCurrLayer--;
        }


        if ((debug & 2) && mesh_.time().writeTime())
        {
            tmp<volScalarField> tfld
            (
                createField
                (
                    mesh_,
                    "allCellTypesWork_Holes" + name(currLayer),
                    allCellTypesWork
                )
            );
            tfld().write();
        }

        if (currLayer == 1)
        {
            setUpFrontOnOversetPatch(allCellTypes, isFront);
        }

        if (currLayer > 1)
        {
            isFront = false;
            setUpFrontOnOversetPatch(allCellTypesWork, isFront);
            setUpFront(allCellTypesWork, isFront);
        }

        // Current interpolation fraction
        scalarField fraction(mesh_.nFaces(), Zero);

        // Ratio between Inter/donor
        scalarField volRatio(allCellTypes.size(), Zero);

        forAll(isFront, facei)
        {
            if (isFront.test(facei))
            {
                fraction[facei] = 1.0;
            }
        }

        scalarField allWeightWork(allCellTypes.size(), Zero);
        bitSet nHoles(allCellTypes.size());

        while (returnReduce(isFront.any(), orOp<bool>()))
        {
            // Interpolate cells on front
            bitSet newIsFront(mesh_.nFaces());
            scalarField newFraction(fraction);

            forAll(isFront, facei)
            {
                if (isFront.test(facei))
                {
                    const label own = mesh_.faceOwner()[facei];

                    if (allCellTypesWork[own] != HOLE)
                    {
                        if (allWeightWork[own] < fraction[facei])
                        {
                            // Cell wants to become interpolated (if sufficient
                            // stencil, otherwise becomes hole)
                            if (allStencil[own].size())
                            {
                                nLayer[own] = currLayer;

                                allWeightWork[own] = fraction[facei];
                                allCellTypesWork[own] = INTERPOLATED;

                                const label donorId = compactStencil[own][0];

                                volRatio[own] = V[own]/compactCellVol[donorId];

                                seedCell
                                (
                                    own,
                                    fraction[facei]-layerRelax,
                                    newIsFront,
                                    newFraction
                                );

                                nHoles[own] = true;
                            }
                            else
                            {
                                allWeightWork[own] = 0.0;
                                allCellTypesWork[own] = HOLE;
                                // Add faces of cell as new front
                                seedCell
                                (
                                    own,
                                    1.0,
                                    newIsFront,
                                    newFraction
                                );
                            }
                        }
                    }
                    if (mesh_.isInternalFace(facei))
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        if (allCellTypesWork[nei] != HOLE)
                        {
                            if (allWeightWork[nei] < fraction[facei])
                            {
                                if (allStencil[nei].size())
                                {
                                    nLayer[nei] = currLayer;

                                    allWeightWork[nei] = fraction[facei];
                                    allCellTypesWork[nei] = INTERPOLATED;

                                    const label donorId =
                                        compactStencil[nei][0];

                                    volRatio[nei] =
                                        V[nei]/compactCellVol[donorId];

                                    seedCell
                                    (
                                        nei,
                                        fraction[facei]-layerRelax,
                                        newIsFront,
                                        newFraction
                                    );
                                }
                                else
                                {
                                    allWeightWork[nei] = 0.0;
                                    allCellTypesWork[nei] = HOLE;
                                    nHoles[nei] = true;
                                    seedCell
                                    (
                                        nei,
                                        1.0,
                                        newIsFront,
                                        newFraction
                                    );
                                }
                            }
                        }
                    }
                }
            }

            syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());
            syncTools::syncFaceList(mesh_, newFraction, maxEqOp<scalar>());

            isFront.transfer(newIsFront);
            fraction.transfer(newFraction);
        }

        if ((debug & 2) && mesh_.time().writeTime())
        {
            tmp<volScalarField> tfld
            (
                createField
                (
                    mesh_,
                    "allCellTypesWork_Layers" + name(currLayer),
                    allCellTypesWork
                )
            );
            tfld().write();
        }

        dAlllWeight.append(allWeightWork);
        dataSet.append(allCellTypesWork);

        // Counting interpolated cells
        label count = 1;
        forAll (volRatio, i)
        {
            if (volRatio[i] > 0)
            {
                count++;
            }
        }
        label nCount(count);
        reduce(nCount, sumOp<label>());

        const scalar aveVol = mag(gSum(volRatio));

        daverageVolRatio.append(aveVol/nCount);

        // Check holes number. A sudden increase occurs when the walk leaks
        // out from the obstacle
        label nTotalHoles(nHoles.count());
        reduce(nTotalHoles, sumOp<label>());
        dHoles.append(nTotalHoles);

        // Check the increase between this layer and the last one
        // if over 50% breaks the layer loop.
        if
        (
            (currLayer > 1)
          & (nTotalHoles > 2.0*dHoles[currLayer - 1])
        )
        {
            break;
        }
    }


    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "walkFront_layers", nLayer)
        );
        tfld().write();
    }

    // Try to find the best averageVolRatio the further from the initial set
    // As this one is next to HOLES
    scalarList averageVolRatio;
    averageVolRatio.transfer(daverageVolRatio);
    label minVolId = findMin(averageVolRatio);

    if (useLayer > -1)
    {
        Info<< nl << " Number of layers : " << averageVolRatio.size() << nl
            << " Average volumetric ratio : " << averageVolRatio << nl
            << " Number of holes cells : " << dHoles << nl
            << " Using layer : " << useLayer <<  nl << endl;
    }

    if (useLayer != -1)
    {
        minVolId = useLayer - 1;
    }

    allCellTypes.transfer(dataSet[minVolId]);
    allWeight.transfer(dAlllWeight[minVolId]);
}


// * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<cellCellStencil>& iproxy
)
{
    const auto& e = *iproxy;

    const labelUList& cellTypes = e.cellTypes();
    const labelUList& interpolationCells = e.interpolationCells();
    const labelListList& cellStencil = e.cellStencil();

    labelList nCells(cellCellStencil::count(5, cellTypes));

    label nInvalidInterpolated = 0;
    label nLocal = 0;
    label nMixed = 0;
    label nRemote = 0;
    forAll(interpolationCells, i)
    {
        label celli = interpolationCells[i];
        const labelList& slots = cellStencil[celli];

        if (slots.empty())
        {
            nInvalidInterpolated++;
        }

        bool hasLocal = false;
        bool hasRemote = false;

        forAll(slots, sloti)
        {
            if (slots[sloti] >= cellTypes.size())
            {
                hasRemote = true;
            }
            else
            {
                hasLocal = true;
            }
        }

        if (hasRemote)
        {
            if (!hasLocal)
            {
                nRemote++;
            }
            else
            {
                nMixed++;
            }
        }
        else if (hasLocal)
        {
            nLocal++;
        }
    }
    reduce(nLocal, sumOp<label>());
    reduce(nMixed, sumOp<label>());
    reduce(nRemote, sumOp<label>());
    reduce(nInvalidInterpolated, sumOp<label>());

    os  << "Overset analysis : nCells : "
        << returnReduce(cellTypes.size(), sumOp<label>()) << nl
        << incrIndent
        << indent << "calculated   : " << nCells[cellCellStencil::CALCULATED]
        << nl
        << indent << "interpolated : " << nCells[cellCellStencil::INTERPOLATED]
        << " (interpolated from local:" << nLocal
        << "  mixed local/remote:" << nMixed
        << "  remote:" << nRemote
        << "  special:" << nInvalidInterpolated << ")" << nl
        << indent << "hole         : " << nCells[cellCellStencil::HOLE] << nl;

    if (nCells[cellCellStencil::SPECIAL] || nCells[cellCellStencil::POROUS])
    {
        os  << indent << "special      : " << nCells[cellCellStencil::SPECIAL]
            << nl
            << indent << "porous       : " << nCells[cellCellStencil::POROUS]
            << nl;
    }
    os  << decrIndent << endl;

    return os;
}


// ************************************************************************* //
