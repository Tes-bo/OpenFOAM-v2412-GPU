/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "decompositionInformation.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionInformation::populate
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const labelUList& decomp,
    const label nDomain
)
{
    nDomains_ = nDomain;

    distrib_.resize_nocopy(nDomain);

    for (labelList& subdist : distrib_)
    {
        subdist.resize_nocopy(nDomain);
        subdist = Foam::zero{};
    }

    // Protect against zero-sized offset list
    const label numCells = Foam::max(0, (xadj.size()-1));

    for (label celli = 0; celli < numCells; ++celli)
    {
        const label ownProc = decomp[celli];

        labelList& subdist = distrib_[ownProc];

        // Number of cells
        ++subdist[ownProc];

        for (label i = xadj[celli]; i < xadj[celli+1]; ++i)
        {
            const label neiProc = decomp[adjncy[i]];

            if (neiProc != ownProc)
            {
                // Number of processor faces
                ++subdist[neiProc];
            }
        }
    }

    // Build summary

    labelList cellsCount(nDomains_, Foam::zero{});
    labelList neighCount(nDomains_, Foam::zero{});
    labelList facesCount(nDomains_, Foam::zero{});

    forAll(distrib_, ownProc)
    {
        const labelList& subdist = distrib_[ownProc];

        cellsCount[ownProc] = subdist[ownProc];

        forAll(subdist, neiProc)
        {
            const label n = subdist[neiProc];

            if (n && ownProc != neiProc)
            {
                ++neighCount[ownProc];
                facesCount[ownProc] += n;
            }
        }
    }

    const label n2 = (nDomains_ / 2);

    sort(cellsCount);
    cellsInfo_.min    = cellsCount.first();
    cellsInfo_.max    = cellsCount.last();
    cellsInfo_.median = cellsCount[n2];

    sort(neighCount);
    neighInfo_.min    = neighCount.first();
    neighInfo_.max    = neighCount.last();
    neighInfo_.median = neighCount[n2];

    sort(facesCount);
    facesInfo_.min    = facesCount.first();
    facesInfo_.max    = facesCount.last();
    facesInfo_.median = facesCount[n2];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionInformation::decompositionInformation
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const labelUList& decomp,
    const label nDomains
)
:
    distrib_(),
    nDomains_(0)
{
    populate(adjncy, xadj, decomp, nDomains);
}


Foam::decompositionInformation::decompositionInformation
(
    const CompactListList<label>& cellCells,
    const labelUList& decomp,
    const label nDomains
)
:
    distrib_(),
    nDomains_(0)
{
    populate(cellCells.values(), cellCells.offsets(), decomp, nDomains);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::decompositionInformation::clear()
{
    distrib_.clear();
    cellsInfo_.clear();
    neighInfo_.clear();
    facesInfo_.clear();
}


void Foam::decompositionInformation::printSummary(Ostream& os) const
{
    os  << "Cells "; cellsInfo_.print(os) << nl;
    os  << "Neigh "; neighInfo_.print(os)<< nl;
    os  << "Faces "; facesInfo_.print(os)<< nl;
}


void Foam::decompositionInformation::printDetails(Ostream& os) const
{
    os  << "Decomposition details with (proc faces) "
        "for each processor connection" << nl << nl;

    forAll(distrib_, ownProc)
    {
        const labelList& subdist = distrib_[ownProc];

        // First pass:
        label neighCount = 0;
        label facesCount = 0;

        forAll(subdist, neiProc)
        {
            const label n = subdist[neiProc];

            if (n && ownProc != neiProc)
            {
                ++neighCount;
                facesCount += n;
            }
        }

        os  << "Part[" << ownProc << "] cells:" << subdist[ownProc]
            << " neigh:" << neighCount
            << " faces:" << facesCount;

        // Second pass with details:
        if (facesCount)
        {
            os  << ' ';

            forAll(subdist, neiProc)
            {
                const label n = subdist[neiProc];

                if (n && ownProc != neiProc)
                {
                    os  << " (" << neiProc << ' ' << n << ')';
                }
            }
        }

        os  << nl;
    }
}


void Foam::decompositionInformation::printAll(Ostream& os) const
{
    printDetails(os);
    printSummary(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::decompositionInformation::stats::print(Ostream& os) const
{
    os  << "min:" << this->min
        << " max:" << this->max
        << " median:" << this->median;

    if (this->median)
    {
        const scalar ratio = scalar(100*this->max)/this->median;

        os  << " (" << ratio << "%)";
    }

    return os;
}


// ************************************************************************* //
