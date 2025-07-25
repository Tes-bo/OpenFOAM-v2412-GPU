/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "AMIWeights.H"
#include "fvMesh.H"
#include "foamVtkSurfaceWriter.H"
#include "PatchTools.H"
#include "cyclicACMIPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(AMIWeights, 0);
    addToRunTimeSelectionTable(functionObject, AMIWeights, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::AMIWeights::writeFileHeader(Ostream& os)
{
    writeHeader(os, "AMI");

    writeCommented(os, "Time");
    forAll(patchIDs_, patchi)
    {
        writeTabbed(os, "Patch");
        writeTabbed(os, "nbr_patch");

        if (UPstream::parRun())
        {
            writeTabbed(os, "distributed");
        }

        writeTabbed(os, "src_min_weight");
        writeTabbed(os, "src_max_weight");
        writeTabbed(os, "src_average_weight");
        writeTabbed(os, "src_min_neighbours");
        writeTabbed(os, "src_max_neighbours");
        writeTabbed(os, "src_average_neighbours");
        writeTabbed(os, "tgt_min_weight");
        writeTabbed(os, "tgt_max_weight");
        writeTabbed(os, "tgt_average_weight");
        writeTabbed(os, "tgt_min_neighbours");
        writeTabbed(os, "tgt_max_neighbours");
        writeTabbed(os, "tgt_average_neighbours");
    }

    os  << endl;
}


void Foam::functionObjects::AMIWeights::reportPatch
(
    const cyclicAMIPolyPatch& pp
)
{
    const word& nbrPatchName = pp.neighbPatchName();

    const Switch distributed = pp.AMI().distributed();

    const scalarField& srcWeightsSum = pp.AMI().srcWeightsSum();
    const scalar srcMinWeight = gMin(srcWeightsSum);
    const scalar srcMaxWeight = gMax(srcWeightsSum);
    const scalar srcAveWeight = gAverage(srcWeightsSum);

    const labelListList& srcAddress = pp.AMI().srcAddress();
    label srcMinNbr = labelMax;
    label srcMaxNbr = labelMin;
    scalar srcAveNbr = 0;
    for (const labelList& srcFace : srcAddress)
    {
        const label n = srcFace.size();
        srcAveNbr += n;
        srcMinNbr = min(srcMinNbr, n);
        srcMaxNbr = max(srcMaxNbr, n);
    }

    reduce(srcMinNbr, minOp<label>());
    reduce(srcMaxNbr, maxOp<label>());

    srcAveNbr =
        returnReduce(srcAveNbr, sumOp<scalar>())
       /(returnReduce(srcAddress.size(), sumOp<scalar>()) + ROOTVSMALL);

    const scalarField& tgtWeightsSum = pp.AMI().tgtWeightsSum();
    const scalar tgtMinWeight = gMin(tgtWeightsSum);
    const scalar tgtMaxWeight = gMax(tgtWeightsSum);
    const scalar tgtAveWeight = gAverage(tgtWeightsSum);

    const labelListList& tgtAddress = pp.AMI().tgtAddress();
    label tgtMinNbr = labelMax;
    label tgtMaxNbr = labelMin;
    scalar tgtAveNbr = 0;
    for (const labelList& tgtFace : tgtAddress)
    {
        const label n = tgtFace.size();
        tgtAveNbr += n;
        tgtMinNbr = min(tgtMinNbr, n);
        tgtMaxNbr = max(tgtMaxNbr, n);
    }

    reduce(tgtMinNbr, minOp<label>());
    reduce(tgtMaxNbr, maxOp<label>());

    tgtAveNbr =
        returnReduce(tgtAveNbr, sumOp<scalar>())
       /(returnReduce(tgtAddress.size(), sumOp<scalar>()) + ROOTVSMALL);

    file()
        << mesh_.time().timeName() << tab
        << pp.name() << tab
        << nbrPatchName << tab;


    if (UPstream::parRun())
    {
        file() << distributed << tab;
    }

    file()
        << srcMinWeight << tab
        << srcMaxWeight << tab
        << srcAveWeight << tab
        << srcMinNbr << tab
        << srcMaxNbr << tab
        << srcAveNbr << tab
        << tgtMinWeight << tab
        << tgtMaxWeight << tab
        << tgtAveWeight << tab
        << tgtMinNbr << tab
        << tgtMaxNbr << tab
        << tgtAveNbr << tab
        << endl;

    Log << "    Patches: " << nl
        << "        Source: " << pp.name() << nl
        << "        Target: " << nbrPatchName << nl;

    if (UPstream::parRun())
    {
        Log << "        Parallel distributed: " << distributed << nl;
    }

    Log << nl;

    const label w = IOstream::defaultPrecision() + 8;

    Log << "                     | " << setw(w) << pp.name()
        << " | " << setw(w) << nbrPatchName << " | " << nl
        << "        min(weight)  | " << setw(w) << srcMinWeight
        << " | " << setw(w) << tgtMinWeight << " | " << nl
        << "        max(weight)  | " << setw(w) << srcMaxWeight
        << " | " << setw(w) << tgtMaxWeight << " | " << nl
        << "        ave(weight)  | " << setw(w) << srcAveWeight
        << " | " << setw(w) << tgtAveWeight << " | " << nl
        << "        min(address) | " << setw(w) << srcMinNbr
        << " | " << setw(w) << tgtMinNbr << " | " << nl
        << "        max(address) | " << setw(w) << srcMaxNbr
        << " | " << setw(w) << tgtMaxNbr << " | " << nl
        << "        ave(address) | " << setw(w) << srcAveNbr
        << " | " << setw(w) << tgtAveNbr << " | " << nl
        << endl;

    setResult(pp.name() + ":src", pp.name());
    setResult(pp.name() + ":tgt", nbrPatchName);
    setResult(pp.name() + ":src:min(weight)", srcMinWeight);
    setResult(pp.name() + ":src:max(weight)", srcMaxWeight);
    setResult(pp.name() + ":src:ave(weight)", srcAveWeight);
    setResult(pp.name() + ":src:min(address)", srcMinNbr);
    setResult(pp.name() + ":src:max(address)", srcMaxNbr);
    setResult(pp.name() + ":src:ave(address)", srcAveNbr);
    setResult(pp.name() + ":tgt:min(weight)", tgtMinWeight);
    setResult(pp.name() + ":tgt:max(weight)", tgtMaxWeight);
    setResult(pp.name() + ":tgt:ave(weight)", tgtAveWeight);
    setResult(pp.name() + ":tgt:min(address)", tgtMinNbr);
    setResult(pp.name() + ":tgt:max(address)", tgtMaxNbr);
    setResult(pp.name() + ":tgt:ave(address)", tgtAveNbr);
}


void Foam::functionObjects::AMIWeights::writeWeightField
(
    const cyclicAMIPolyPatch& cpp,
    const scalarField& weightSum,
    const word& side
) const
{
    // Collect geometry
    labelList pointToGlobal;
    labelList uniqueMeshPointLabels;
    autoPtr<globalIndex> globalPoints;
    autoPtr<globalIndex> globalFaces;
    faceList mergedFaces;
    pointField mergedPoints;
    Foam::PatchTools::gatherAndMerge
    (
        mesh_,
        cpp.localFaces(),
        cpp.meshPoints(),
        cpp.meshPointMap(),

        pointToGlobal,
        uniqueMeshPointLabels,
        globalPoints,
        globalFaces,

        mergedFaces,
        mergedPoints
    );

    // Collect field
    scalarField mergedWeights;
    globalFaces().gather(weightSum, mergedWeights);

    const bool isACMI = isA<cyclicACMIPolyPatch>(cpp);

    scalarField mergedMask;
    if (isACMI)
    {
        const cyclicACMIPolyPatch& pp = refCast<const cyclicACMIPolyPatch>(cpp);

        globalFaces().gather(pp.mask(), mergedMask);
    }

    if (Pstream::master())
    {
        instant inst(mesh_.time().value(), mesh_.time().timeName());

        vtk::surfaceWriter writer
        (
            mergedPoints,
            mergedFaces,
            (baseTimeDir()/cpp.name() + "_" + side),
            false // serial: master-only
        );

        writer.setTime(inst);
        writer.writeTimeValue();
        writer.writeGeometry();

        writer.beginCellData(1 + (isACMI ? 1 : 0));
        writer.write("weightsSum", mergedWeights);

        if (isACMI)
        {
            writer.write("mask", mergedMask);
        }
    }
}


void Foam::functionObjects::AMIWeights::writeWeightFields
(
    const cyclicAMIPolyPatch& cpp
) const
{
    if (cpp.owner())
    {
        writeWeightField(cpp, cpp.AMI().srcWeightsSum(), "src");
        writeWeightField(cpp.neighbPatch(), cpp.AMI().tgtWeightsSum(), "tgt");
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::AMIWeights::AMIWeights
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    writeFields_(false),
    patchIDs_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::AMIWeights::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        patchIDs_.clear();
        labelHashSet ids;

        for (const polyPatch& pp : mesh_.boundaryMesh())
        {
            const auto* amicpp = isA<cyclicAMIPolyPatch>(pp);

            if (amicpp && amicpp->owner())
            {
                ids.insert(pp.index());
            }
        }

        patchIDs_ = ids.sortedToc();

        writeFileHeader(file());

        writeFields_ = dict.get<bool>("writeFields");

        return true;
    }

    return false;
}


bool Foam::functionObjects::AMIWeights::execute()
{
    return true;
}


bool Foam::functionObjects::AMIWeights::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    for (const label patchi : patchIDs_)
    {
        const polyPatch& pp = pbm[patchi];
        const auto& cpp = static_cast<const cyclicAMIPolyPatch&>(pp);

        reportPatch(cpp);

        if (writeFields_)
        {
            writeWeightFields(cpp);
        }
    }

    return true;
}


// ************************************************************************* //
