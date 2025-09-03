/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "ensightMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// DIY flatOutput without a leading size indicator
static Ostream& printPatterns(Ostream& os, const wordRes& list)
{
    os << token::BEGIN_LIST;

    bool sep = false;

    for (const wordRe& item : list)
    {
        if (sep) os << token::SPACE;
        sep = true;

        os << item;
    }
    os << token::END_LIST;

    return os;
}

} // End namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::options::options()
:
    lazy_(false),
    internal_(true),
    boundary_(true),
    cellZones_(true),
    patchInclude_(),
    patchExclude_(),
    cellZoneInclude_(),
    faceZoneInclude_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightMesh::options::lazy() const noexcept
{
    return lazy_;
}


bool Foam::ensightMesh::options::useInternalMesh() const noexcept
{
    return internal_;
}


bool Foam::ensightMesh::options::useBoundaryMesh() const noexcept
{
    return boundary_;
}


bool Foam::ensightMesh::options::useCellZones() const noexcept
{
    return cellZones_;
}


bool Foam::ensightMesh::options::useFaceZones() const noexcept
{
    return faceZoneInclude_.size();
}


void Foam::ensightMesh::options::reset()
{
    internal_ = true;
    boundary_ = true;
    cellZones_ = true;
    patchInclude_.clear();
    patchExclude_.clear();
    faceZoneInclude_.clear();
    cellZoneInclude_.clear();
}


bool Foam::ensightMesh::options::lazy(bool on) noexcept
{
    bool old(lazy_);
    lazy_ = on;
    return old;
}


bool Foam::ensightMesh::options::useInternalMesh(bool on) noexcept
{
    bool old(internal_);
    internal_ = on;
    return old;
}


bool Foam::ensightMesh::options::useBoundaryMesh(bool on)
{
    bool old(boundary_);
    boundary_ = on;

    if (!boundary_)
    {
        if (patchInclude_.size())
        {
            patchInclude_.clear();

            WarningInFunction
                << "Deactivating boundary, removed old patch selection"
                << endl;
        }
    }

    return old;
}


bool Foam::ensightMesh::options::useCellZones(bool on)
{
    bool old(cellZones_);
    cellZones_ = on;

    if (!cellZones_ && cellZoneInclude_.size())
    {
        cellZoneInclude_.clear();

        WarningInFunction
            << "Deactivating cellZones, removed old zone selection"
            << endl;
    }

    return old;
}


void Foam::ensightMesh::options::patchSelection
(
    const UList<wordRe>& patterns
)
{
    patchInclude_ = wordRes(patterns);

    if (!boundary_ && patchInclude_.size())
    {
        patchInclude_.clear();

        WarningInFunction
            << "Ignoring patch selection, boundary is disabled"
            << endl;
    }
}


void Foam::ensightMesh::options::patchSelection
(
    List<wordRe>&& patterns
)
{
    patchInclude_ = wordRes(std::move(patterns));

    if (!boundary_ && patchInclude_.size())
    {
        patchInclude_.clear();

        WarningInFunction
            << "Ignoring patch selection, boundary is disabled"
            << endl;
    }
}


void Foam::ensightMesh::options::patchExclude
(
    const UList<wordRe>& patterns
)
{
    patchExclude_ = wordRes(patterns);
}


void Foam::ensightMesh::options::patchExclude
(
    List<wordRe>&& patterns
)
{
    patchExclude_ = wordRes(std::move(patterns));
}


void Foam::ensightMesh::options::faceZoneSelection
(
    const UList<wordRe>& patterns
)
{
    faceZoneInclude_ = wordRes(patterns);
}


void Foam::ensightMesh::options::faceZoneSelection
(
    List<wordRe>&& patterns
)
{
    faceZoneInclude_ = wordRes(std::move(patterns));
}


void Foam::ensightMesh::options::cellZoneSelection
(
    const UList<wordRe>& patterns
)
{
    cellZoneInclude_ = wordRes(patterns);

    if (!cellZones_ && cellZoneInclude_.size())
    {
        cellZoneInclude_.clear();

        WarningInFunction
            << "Ignoring cellZone selection, cellZones are disabled"
            << endl;
    }
}


void Foam::ensightMesh::options::cellZoneSelection
(
    List<wordRe>&& patterns
)
{
    cellZoneInclude_ = wordRes(std::move(patterns));

    if (!cellZones_ && cellZoneInclude_.size())
    {
        cellZoneInclude_.clear();

        WarningInFunction
            << "Ignoring cellZone selection, cellZones are disabled"
            << endl;
    }
}


void Foam::ensightMesh::options::print(Ostream& os) const
{
    os << "internal: " << Switch::name(internal_) << nl;
    os << "cellZones: " << Switch::name(cellZones_) << nl;
    if (cellZones_)
    {
        os.incrIndent();
        if (!cellZoneInclude_.empty())
        {
            os.writeKeyword("include");
            printPatterns(os, cellZoneInclude_) << nl;
        }
        os.decrIndent();
    }

    os << "boundary: " << Switch::name(boundary_) << nl;
    if (boundary_)
    {
        os.incrIndent();
        if (!patchInclude_.empty())
        {
            os.writeKeyword("include");
            printPatterns(os, patchInclude_) << nl;
        }
        if (!patchExclude_.empty())
        {
            os.writeKeyword("exclude");
            printPatterns(os, patchExclude_) << nl;
        }
        os.decrIndent();
    }

    os << "faceZones: " << Switch::name(useFaceZones()) << nl;
    if (useFaceZones())
    {
        os.incrIndent();
        if (!faceZoneInclude_.empty())
        {
            os.writeKeyword("include");
            printPatterns(os, faceZoneInclude_) << nl;
        }
        os.decrIndent();
    }
}


// ************************************************************************* //
