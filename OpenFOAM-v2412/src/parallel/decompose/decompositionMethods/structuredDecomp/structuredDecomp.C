/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "structuredDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structuredDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        structuredDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structuredDecomp::structuredDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict),
    methodDict_(findCoeffsDict(typeName + "Coeffs", selectionType::MANDATORY)),
    patches_(methodDict_.get<wordRes>("patches"))
{
    methodDict_.set("numberOfSubdomains", nDomains());
    method_ = decompositionMethod::New(methodDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::structuredDecomp::parallelAware() const
{
    return method_().parallelAware();
}


Foam::labelList Foam::structuredDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& cc,
    const scalarField& cWeights
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelHashSet patchIDs(pbm.patchSet(patches_));

    label nFaces = 0;
    for (const label patchi : patchIDs)
    {
        nFaces += pbm[patchi].size();
    }

    // Extract a submesh.
    labelHashSet patchCells(2*nFaces);
    for (const label patchi : patchIDs)
    {
        patchCells.insert(pbm[patchi].faceCells());
    }

    // Subset the layer of cells next to the patch
    fvMeshSubset subsetter
    (
        dynamic_cast<const fvMesh&>(mesh),
        patchCells
    );
    const fvMesh& subMesh = subsetter.subMesh();
    pointField subCc(cc, subsetter.cellMap());
    scalarField subWeights;
    if (cWeights.size() == cc.size())
    {
        subWeights = scalarField(cWeights, subsetter.cellMap());
    }

    // Decompose the layer of cells
    labelList subDecomp(method_().decompose(subMesh, subCc, subWeights));


    // Transfer to final decomposition
    labelList finalDecomp(cc.size(), -1);
    forAll(subDecomp, i)
    {
        finalDecomp[subsetter.cellMap()[i]] = subDecomp[i];
    }

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh.nCells());
    List<topoDistanceData<label>> faceData(mesh.nFaces());

    // Start of changes
    labelList patchFaces(nFaces);
    List<topoDistanceData<label>> patchData(nFaces);
    nFaces = 0;
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];
        const labelUList& fc = pp.faceCells();
        forAll(fc, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData<label>(0, finalDecomp[fc[i]]);
            nFaces++;
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );

    // And extract
    bool haveWarned = false;
    forAll(finalDecomp, celli)
    {
        if (!cellData[celli].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some cells, e.g. cell " << celli
                    << " at " << mesh.cellCentres()[celli] << nl
                    << "Assigning  these cells to domain 0." << endl;
                haveWarned = true;
            }
            finalDecomp[celli] = 0;
        }
        else
        {
            finalDecomp[celli] = cellData[celli].data();
        }
    }

    return finalDecomp;
}


// ************************************************************************* //
