/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "labelRanges.H"
#include <numeric>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Print range for debugging purposes
static Ostream& printRange(Ostream& os, const labelRange& range)
{
    if (range.empty())
    {
        os  << "empty";
    }
    else
    {
        os  << range << " = " << range.min() << ':' << range.max();
    }
    return os;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::labelRanges::insertBefore
(
    const label insert,
    const labelRange& range
)
{
    auto& list = ranges_;

    // Insert via copying up
    label nElem = list.size();

    if (labelRange::debug)
    {
        Info<< "before insert "
            << nElem << " elements, insert at " << insert << nl
            << list << endl;
    }

    list.resize(nElem+1);

    if (labelRange::debug)
    {
        Info<< "copy between " << nElem << " and " << insert << nl;
    }

    for (label i = nElem-1; i >= insert; --i)
    {
        if (labelRange::debug)
        {
            Info<< "copy from " << (i) << " to " << (i+1) << nl;
        }

        list[i+1] = list[i];
    }

    // Finally insert the range
    if (labelRange::debug)
    {
        Info<< "finally insert the range at " << insert << nl;
    }

    list[insert] = range;
}


void Foam::labelRanges::purgeEmpty()
{
    auto& list = ranges_;

    // Purge empty ranges by copying down
    label nElem = 0;

    forAll(list, elemi)
    {
        if (!list[elemi].empty())
        {
            if (nElem != elemi)
            {
                list[nElem] = list[elemi];
            }
            ++nElem;
        }
    }

    // Truncate
    list.resize(nElem);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelRanges::labelRanges(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::labelRanges::add(const labelRange& range)
{
    auto& list = ranges_;

    if (range.empty())
    {
        return false;
    }
    else if (list.empty())
    {
        list.push_back(range);
        return true;
    }

    // Find the correct place for insertion
    forAll(list, elemi)
    {
        labelRange& currRange = list[elemi];

        if (currRange.overlaps(range, true))
        {
            // absorb into the existing (adjacent/overlapping) range
            currRange.join(range);

            // Might connect with the next following range(s)
            for (; elemi < list.size()-1; ++elemi)
            {
                labelRange& nextRange = list[elemi+1];
                if (currRange.overlaps(nextRange, true))
                {
                    currRange.join(nextRange);
                    nextRange.reset();
                }
                else
                {
                    break;
                }
            }

            // Done - remove any empty ranges that might have been created
            purgeEmpty();
            return true;
            break;
        }
        else if (range < currRange)
        {
            insertBefore(elemi, range);
            return true;
            break;
        }
    }


    // not found: simply append
    list.push_back(range);

    return true;
}


bool Foam::labelRanges::remove(const labelRange& range)
{
    auto& list = ranges_;
    bool status = false;

    if (range.empty() || list.empty())
    {
        return status;
    }

    forAll(list, elemi)
    {
        labelRange& currRange = list[elemi];

        if (range.min() > currRange.min())
        {
            if (range.max() < currRange.max())
            {
                // Removal of range fragments of currRange

                if (labelRange::debug)
                {
                    Info<< "Fragment removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                // left-hand-side fragment: insert before current range
                label lower = currRange.min();
                label upper = range.min() - 1;

                labelRange fragment(lower, upper - lower + 1);

                // right-hand-side fragment
                lower = range.max() + 1;
                upper = currRange.max();

                currRange.reset(lower, upper - lower + 1);
                currRange.clampSize();
                status = true;
                insertBefore(elemi, fragment);

                if (labelRange::debug)
                {
                    Info<< "fragment ";
                    printRange(Info, fragment) << endl;
                    Info<< "yields ";
                    printRange(Info, currRange) << endl;
                }

                // fragmentation can only affect a single range
                // thus we are done
                break;
            }
            else if (range.min() <= currRange.max())
            {
                // keep left-hand-side, remove right-hand-side

                if (labelRange::debug)
                {
                    Info<< "RHS removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                const label lower = currRange.min();
                const label upper = range.min() - 1;

                currRange.reset(lower, upper - lower + 1);
                currRange.clampSize();
                status = true;

                if (labelRange::debug)
                {
                    Info<< "yields ";
                    printRange(Info, currRange) << endl;
                }
            }
        }
        else if (range.min() <= currRange.min())
        {
            if (range.max() >= currRange.min())
            {
                // remove left-hand-side, keep right-hand-side

                if (labelRange::debug)
                {
                    Info<< "LHS removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                const label lower = range.max() + 1;
                const label upper = currRange.max();

                currRange.reset(lower, upper - lower + 1);
                currRange.clampSize();
                status = true;

                if (labelRange::debug)
                {
                    Info<< "yields ";
                    printRange(Info, currRange) << endl;
                }
            }
        }
    }

    purgeEmpty();

    return status;
}


Foam::List<Foam::label> Foam::labelRanges::labels() const
{
    label total = 0;
    for (const labelRange& range : ranges_)
    {
        if (range.size() > 0)  // Ignore negative size (paranoid)
        {
            total += range.size();
        }
    }

    if (!total)
    {
        // Skip this check?
        return List<label>();
    }

    List<label> result(total);

    auto* iter = result.begin();

    for (const labelRange& range : ranges_)
    {
        const label len = range.size();

        if (len > 0)  // Ignore negative size (paranoid)
        {
            std::iota(iter, (iter + len), range.start());
            iter += len;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::label Foam::labelRanges::operator[](const label i) const
{
    if (i < 0) return -1;

    label subIdx = i;

    for (const labelRange& range : ranges_)
    {
        if (subIdx < range.size())
        {
            return (range.start() + subIdx);
        }
        else
        {
            subIdx -= range.size();
        }
    }

    return -1;  // Not found
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::labelRanges::readList(Istream& is)
{
    return ranges_.readList(is);
}


Foam::Ostream& Foam::labelRanges::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    return ranges_.writeList(os, shortLen);
}


Foam::Istream& Foam::operator>>(Istream& is, labelRanges& list)
{
    return list.readList(is);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const labelRanges& list)
{
    return list.writeList
    (
        os,
        Detail::ListPolicy::short_length<labelRange>::value
    );
}


// ************************************************************************* //
