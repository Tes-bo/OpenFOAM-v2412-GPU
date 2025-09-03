/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "uniformBin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(uniformBin, 0);
    addToRunTimeSelectionTable(binModel, uniformBin, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::binModels::uniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Use geometry limits if not specified by the user
    {
        // Determine extents of patches/cells
        boundBox geomLimits;

        for (const label patchi : patchIDs_)
        {
            vectorField pts
            (
                coordSysPtr_->localPosition(pbm[patchi].faceCentres())
            );

            MinMax<vector> limits(pts);

            geomLimits.add(limits.min());
            geomLimits.add(limits.max());
        }

        for (const label zonei : cellZoneIDs_)
        {
            const cellZone& cZone = mesh_.cellZones()[zonei];
            const vectorField pts
            (
                coordSysPtr_->localPosition(vectorField(mesh_.C(), cZone))
            );

            MinMax<vector> limits(pts);

            geomLimits.add(limits.min());
            geomLimits.add(limits.max());
        }

        // Globally consistent
        geomLimits.reduce();

        // Slightly boost max so that region of interest is fully within bounds
        // TBD: could also adjust min?
        const vector adjust(1e-4*geomLimits.span());
        geomLimits.max() += adjust;

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            // Use geometry limits if not specified by the user
            if (binLimits_.min()[cmpt] == GREAT)
            {
                binLimits_.min()[cmpt] = geomLimits.min()[cmpt];
            }
            if (binLimits_.max()[cmpt] == GREAT)
            {
                binLimits_.max()[cmpt] = geomLimits.max()[cmpt];
            }
        }
    }

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (binLimits_.min()[cmpt] > binLimits_.max()[cmpt])
        {
            FatalErrorInFunction
                << "Max bounds must be greater than min bounds" << nl
                << "    direction   = " << cmpt << nl
                << "    min         = " << binLimits_.min()[cmpt] << nl
                << "    max         = " << binLimits_.max()[cmpt] << nl
                << exit(FatalError);
        }

        //- Compute bin widths in binning directions
        binWidth_[cmpt] =
        (
            (binLimits_.max()[cmpt] - binLimits_.min()[cmpt])
          / scalar(nBins_[cmpt])
        );

        if (binWidth_[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Bin widths must be greater than zero" << nl
                << "    direction = " << cmpt << nl
                << "    min bound = " << binLimits_.min()[cmpt] << nl
                << "    max bound = " << binLimits_.max()[cmpt] << nl
                << "    bin width = " << binWidth_[cmpt] << nl
                << exit(FatalError);
        }
    }

    setBinsAddressing();
}


Foam::labelList Foam::binModels::uniformBin::binAddr(const vectorField& d) const
{
    labelList binIndices(d.size(), -1);

    forAll(d, i)
    {
        // Avoid elements outside of the bin
        bool faceInside = true;
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if
            (
                d[i][cmpt] < binLimits_.min()[cmpt]
             || d[i][cmpt] > binLimits_.max()[cmpt]
            )
            {
                faceInside = false;
                break;
            }
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the element
            Vector<label> n(Zero);
            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                label bini = floor
                (
                    (d[i][cmpt] - binLimits_.min()[cmpt])/binWidth_[cmpt]
                );

                n[cmpt] = min(max(bini, 0), nBins_[cmpt] - 1);
            }

            // Order: (e1, e2, e3), the first varies the fastest
            binIndices[i] = n.x() + nBins_[0]*n.y() + nBins_[0]*nBins_[1]*n.z();
        }
        else
        {
            binIndices[i] = -1;
        }
    }

    return binIndices;
}


void Foam::binModels::uniformBin::setBinsAddressing()
{
    faceToBin_.resize_nocopy(mesh_.nBoundaryFaces());
    faceToBin_ = -1;

    for (const label patchi : patchIDs_)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const label i0 = pp.start() - mesh_.nInternalFaces();

        SubList<label>(faceToBin_, pp.size(), i0) =
            binAddr(coordSysPtr_->localPosition(pp.faceCentres()));
    }

    cellToBin_.resize_nocopy(mesh_.nCells());
    cellToBin_ = -1;

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];
        labelList bins
        (
            binAddr(coordSysPtr_->localPosition(vectorField(mesh_.C(), cZone)))
        );

        forAll(cZone, i)
        {
            const label celli = cZone[i];
            cellToBin_[celli] = bins[i];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::uniformBin::uniformBin
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& outputPrefix
)
:
    binModel(dict, mesh, outputPrefix),
    nBins_(Zero),
    binWidth_(Zero),
    binLimits_(vector::uniform(GREAT))
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::uniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    Info<< "    Activating a set of uniform bins" << endl;

    const dictionary& binDict = dict.subDict("binData");

    nBins_ = binDict.get<Vector<label>>("nBin");

    for (const label n : nBins_)
    {
        nBin_ *= n;
    }

    if (nBin_ <= 0)
    {
        FatalIOErrorInFunction(binDict)
            << "Number of bins must be greater than zero" << nl
            << "    e1 bins = " << nBins_.x() << nl
            << "    e2 bins = " << nBins_.y() << nl
            << "    e3 bins = " << nBins_.z()
            << exit(FatalIOError);
    }

    Info<< "    - Employing:" << nl
        << "        " << nBins_.x() << " e1 bins," << nl
        << "        " << nBins_.y() << " e2 bins," << nl
        << "        " << nBins_.z() << " e3 bins"
        << endl;

    cumulative_ = binDict.getOrDefault<bool>("cumulative", false);
    Info<< "    - cumulative    : " << cumulative_ << endl;
    Info<< "    - decomposePatchValues    : " << decomposePatchValues_ << endl;

    const dictionary* minMaxDictPtr = binDict.findDict("minMax");

    if (minMaxDictPtr)
    {
        const auto& minMaxDict = *minMaxDictPtr;

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            const word ei("e" + Foam::name(cmpt));

            scalarMinMax range;

            if (minMaxDict.readIfPresent(ei, range))
            {
                binLimits_.min()[cmpt] = range.min();
                binLimits_.max()[cmpt] = range.max();

                Info<< "    - " << ei << " min/max    : " << range << nl;
            }
        }
    }
    Info<< endl;

    initialise();

    return true;
}


void Foam::binModels::uniformBin::apply()
{
    forAll(fieldNames_, i)
    {
        const bool ok =
        (
            processField<scalar>(i)
         || processField<vector>(i)
         || processField<sphericalTensor>(i)
         || processField<symmTensor>(i)
         || processField<tensor>(i)
        );

        if (!ok)
        {
            WarningInFunction
                << "Unable to find field " << fieldNames_[i] << endl;
        }
    }

    writtenHeader_ = true;
}


void Foam::binModels::uniformBin::updateMesh(const mapPolyMesh& mpm)
{}


void Foam::binModels::uniformBin::movePoints(const polyMesh& mesh)
{
    setBinsAddressing();
}


// ************************************************************************* //
