/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "targetVolumeToCell.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "plane.H"
#include "bitSet.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(targetVolumeToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, targetVolumeToCell, word);
    addToRunTimeSelectionTable(topoSetSource, targetVolumeToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, targetVolumeToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, targetVolumeToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        targetVolumeToCell,
        word,
        targetVolume
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        targetVolumeToCell,
        istream,
        targetVolume
    );
}


Foam::topoSetSource::addToUsageTable Foam::targetVolumeToCell::usage_
(
    targetVolumeToCell::typeName,
    "\n    Usage: targetVolumeToCell (nx ny nz)\n\n"
    "    Adjust plane until obtained selected volume\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::targetVolumeToCell::volumeOfSet
(
    const bitSet& selected
) const
{
    scalar sumVol = 0.0;

    // Loop over selected cells only
    for (const label celli : selected)
    {
        sumVol += mesh_.cellVolumes()[celli];
    }

    return returnReduce(sumVol, sumOp<scalar>());
}


Foam::label Foam::targetVolumeToCell::selectCells
(
    const scalar normalComp,
    const bitSet& maskSet,
    bitSet& selected
) const
{
    selected.resize(mesh_.nCells());
    selected = false;

    label nSelected = 0;

    forAll(mesh_.cellCentres(), celli)
    {
        const point& cc = mesh_.cellCentres()[celli];

        if (maskSet.test(celli) && ((cc & normal_) < normalComp))
        {
            selected.set(celli);
            ++nSelected;
        }
    }

    return returnReduce(nSelected, sumOp<label>());
}


void Foam::targetVolumeToCell::combine(topoSet& set, const bool add) const
{
    if (vol_ <= 0)
    {
        // Select no cells
        return;
    }

    bitSet maskSet(mesh_.nCells(), true);
    label nTotCells = mesh_.globalData().nTotalCells();
    if (maskSetName_.size())
    {
        // Read cellSet
        if (verbose_)
        {
            Info<< "    Operating on subset defined by cellSet "
                << maskSetName_ << endl;
        }

        cellSet subset(mesh_, maskSetName_, IOobject::NO_REGISTER);
        const labelHashSet& cellLabels = subset;

        maskSet = false;
        maskSet.setMany(cellLabels.begin(), cellLabels.end());

        nTotCells = returnReduce(cellLabels.size(), sumOp<label>());
    }


    // Get plane for min,max volume.
    // Planes all have base (0 0 0) and fixed normal so work only on normal
    // component.

    scalar maxComp = -GREAT;
    label maxCells = 0;
    //scalar maxVol = 0;
    scalar minComp = GREAT;
    {
        const boundBox& bb = mesh_.bounds();
        pointField points(bb.points());

        //label minPointi = -1;
        label maxPointi = -1;
        forAll(points, pointi)
        {
            const scalar c = (points[pointi] & normal_);
            if (c > maxComp)
            {
                maxComp = c;
                maxPointi = pointi;
            }
            else if (c < minComp)
            {
                minComp = c;
                //minPointi = pointi;
            }
        }

        bitSet maxSelected(mesh_.nCells());
        maxCells = selectCells(maxComp, maskSet, maxSelected);
        //maxVol = volumeOfSet(maxSelected);

        // Check that maxPoint indeed selects all cells
        if (maxCells != nTotCells)
        {
            WarningInFunction
                << "Plane " << plane(points[maxPointi], normal_)
                << " selects " << maxCells
                << " cells instead of all " << nTotCells
                << " cells. Results might be wrong." << endl;
        }
    }


    // Bisection
    // ~~~~~~~~~

    bitSet selected(mesh_.nCells());
    label nSelected = -1;
    scalar selectedVol = 0.0;
    //scalar selectedComp = 0.0;


    scalar low = minComp;
    scalar high = maxComp;

    const scalar tolerance = SMALL*100*(maxComp-minComp);

    while ((high-low) > tolerance)
    {
        const scalar mid = 0.5*(low + high);

        nSelected = selectCells(mid, maskSet, selected);
        selectedVol = volumeOfSet(selected);

        //Pout<< "High:" << high << " low:" << low << " mid:" << mid << nl
        //    << "    nSelected:" << nSelected << nl
        //    << "    vol      :" << selectedVol << nl
        //    << endl;

        if (selectedVol < vol_)
        {
            low = mid;

            bitSet highSelected(mesh_.nCells());
            label nHigh = selectCells(high, maskSet, selected);
            if (nSelected == nHigh)
            {
                break;
            }
        }
        else
        {
            high = mid;

            bitSet lowSelected(mesh_.nCells());
            label nLow = selectCells(low, maskSet, selected);
            if (nSelected == nLow)
            {
                break;
            }
        }
    }

    nSelected = selectCells(high, maskSet, selected);
    selectedVol = volumeOfSet(selected);

    if (selectedVol < vol_)
    {
        //selectedComp = high;
    }
    else
    {
        nSelected = selectCells(low, maskSet, selected);
        selectedVol = volumeOfSet(selected);

        if (selectedVol < vol_)
        {
            //selectedComp = low;
        }
        else
        {
            WarningInFunction
                << "Did not converge onto plane. " << nl
                << "high plane:"
                << plane(high*normal_, normal_)
                << nl
                << "low plane :"
                << plane(low*normal_, normal_)
                << endl;
        }
    }


    if (verbose_)
    {
        Info<< "    Selected " << nSelected << " with actual volume "
            << selectedVol << endl;
    }

    addOrDelete(set, selected, add);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::targetVolumeToCell::targetVolumeToCell
(
    const polyMesh& mesh,
    const scalar vol,
    const vector& normal,
    const word& maskSetName
)
:
    topoSetCellSource(mesh),
    vol_(vol),
    normal_(normal),
    maskSetName_(maskSetName)
{}


Foam::targetVolumeToCell::targetVolumeToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    vol_(dict.getCheck<scalar>("volume", scalarMinMax::ge(0))),
    normal_(dict.get<vector>("normal")),
    maskSetName_(dict.getOrDefault<word>("set", ""))
{}


Foam::targetVolumeToCell::targetVolumeToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    vol_(readScalar(checkIs(is))),
    normal_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::targetVolumeToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells up to target volume " << vol_
                << " out of total volume "
                << gSum(mesh_.cellVolumes()) << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells up to target volume " << vol_
                << " out of total volume "
                << gSum(mesh_.cellVolumes()) << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
