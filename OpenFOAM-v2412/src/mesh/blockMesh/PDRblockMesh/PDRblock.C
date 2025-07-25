/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "ListOps.H"
#include "emptyPolyPatch.H"
#include "gradingDescriptors.H"

// Output when verbosity is enabled
#undef  Log
#define Log if (verbose_) Info


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::PDRblock::expansionType
>
Foam::PDRblock::expansionNames_
({
    { expansionType::EXPAND_UNIFORM, "uniform" },
    { expansionType::EXPAND_RATIO, "ratio" },
    { expansionType::EXPAND_RELATIVE, "relative"}
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Calculate geometric expansion factor from expansion ratio
inline scalar calcGexp(const scalar expRatio, const label nDiv)
{
    return nDiv > 1 ? pow(expRatio, 1.0/(nDiv - 1)) : 0.0;
}

//- Calculate geometric ratio from relative ratio
inline scalar relativeToGeometricRatio
(
    const scalar expRatio,
    const label nDiv
)
{
    return nDiv > 1 ? pow(expRatio, (nDiv - 1)) : 1.0;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::PDRblock::checkMonotonic
(
    const direction cmpt,
    const UList<scalar>& pts
)
{
    const label len = pts.size();

    if (!len)
    {
        return false;
    }

    const scalar& minVal = pts[0];

    for (label i=1; i < len; ++i)
    {
        if (pts[i] <= minVal)
        {
            FatalErrorInFunction
                << "Points in " << vector::componentNames[cmpt]
                << " direction do not increase monotonically" << nl
                << flatOutput(pts) << nl << nl
                << exit(FatalError);
        }
    }

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRblock::addDefaultPatches()
{
    // Default boundaries with patchi == shapeFacei
    patches_.resize(6);
    for (label patchi=0; patchi < 6; ++patchi)
    {
        boundaryEntry& bentry = patches_.emplace_set(patchi);

        bentry.name_ = "patch" + Foam::name(patchi);
        bentry.type_ = "patch";
        bentry.size_ = 0;
        bentry.faces_ = labelList(one{}, patchi);
    }
}


void Foam::PDRblock::adjustSizes()
{
    // Adjust i-j-k addressing
    sizes().x() = grid_.x().nCells();
    sizes().y() = grid_.y().nCells();
    sizes().z() = grid_.z().nCells();

    if (sizes().x() <= 0 || sizes().y() <= 0 || sizes().z() <= 0)
    {
        // Sanity check. Silently disallow bad sizing
        ijkMesh::clear();

        grid_.x().clear();
        grid_.y().clear();
        grid_.z().clear();

        bounds_.reset();
        edgeLimits_.min() = 0;
        edgeLimits_.max() = 0;
        return;
    }

    // Adjust boundBox
    bounds_ = bounds(grid_.x(), grid_.y(), grid_.z());

    // Min/max edge lengths
    edgeLimits_.clear();

    edgeLimits_.add(grid_.x().edgeLimits());
    edgeLimits_.add(grid_.y().edgeLimits());
    edgeLimits_.add(grid_.z().edgeLimits());
}


void Foam::PDRblock::readGridControl
(
    const direction cmpt,
    const dictionary& dict,
    const scalar scaleFactor,
    expansionType expandType
)
{
    gridControl& ctrl = control_[cmpt];

    // Begin/end nodes for each segment
    scalarList& knots = static_cast<scalarList&>(ctrl);

    // The number of division per segment
    labelList& divisions = ctrl.divisions_;

    // The expansion ratio per segment
    scalarList& expansion = ctrl.expansion_;

    expansion.clear();  // expansion is optional

    Log << "Reading grid control for "
        << vector::componentNames[cmpt] << " direction" << nl;

    // Points
    // ~~~~~~

    dict.readEntry("points", knots);

    const label nSegments = (knots.size()-1);

    if (nSegments < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Must define at least two control points:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }

    // Apply scaling
    if (scaleFactor > 0)
    {
        for (scalar& pt : knots)
        {
            pt *= scaleFactor;
        }
    }

    // Do points increase monotonically?
    checkMonotonic(cmpt, knots);

    Log << "    points : " << flatOutput(knots) << nl;


    // Divisions
    // ~~~~~~~~~

    dict.readEntry("nCells", divisions);

    label nTotalDivs = 0;
    for (const label ndiv : divisions)
    {
        nTotalDivs += ndiv;

        if (ndiv < 1)
        {
            FatalIOErrorInFunction(dict)
                << "Negative or zero divisions:"
                << " in " << vector::componentNames[cmpt]
                << " direction" << nl
                << exit(FatalIOError);
        }
    }

    if (divisions.size() != nSegments)
    {
        FatalIOErrorInFunction(dict)
            << "Mismatch in number of divisions and number of points:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }
    else if (!nTotalDivs)
    {
        FatalIOErrorInFunction(dict)
            << "No grid divisions at all:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }

    Log << "    nCells : " << flatOutput(divisions) << nl;


    // Expansion ratios
    // ~~~~~~~~~~~~~~~~

    dict.readIfPresent("ratios", expansion);

    // Correct expansion ratios - negative is the same as inverse.
    for (scalar& expRatio : expansion)
    {
        if (expRatio < 0)
        {
            expRatio = 1.0/(-expRatio);
        }
    }

    if (expansion.size() && divisions.size() != expansion.size())
    {
        FatalIOErrorInFunction(dict)
            << "Mismatch in number of expansion ratios and divisions:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }

    if (expansion.empty())
    {
        expansion.resize(nSegments, scalar(1));

        if (expandType != expansionType::EXPAND_UNIFORM)
        {
            expandType = expansionType::EXPAND_UNIFORM;

            Log << "Warning: no 'ratios', use uniform spacing" << nl;
        }
    }
    else
    {
        switch (expandType)
        {
            case expansionType::EXPAND_UNIFORM:
            {
                expansion = scalar(1);
                break;
            }

            case expansionType::EXPAND_RATIO:
            {
                Log << "    ratios : " << flatOutput(expansion) << nl;
                break;
            }

            case expansionType::EXPAND_RELATIVE:
            {
                Log << "  relative : " << flatOutput(expansion) << nl;

                auto divIter = divisions.cbegin();

                for (scalar& expRatio : expansion)
                {
                    expRatio = relativeToGeometricRatio(expRatio, *divIter);
                    ++divIter;
                }

                Log << "    ratios : " << flatOutput(expansion) << nl;
                break;
            }
        }
    }


    // Define the grid points

    scalarList& gridPoint = grid_[cmpt];
    gridPoint.resize(nTotalDivs+1);

    label start = 0;

    for (label segmenti=0; segmenti < nSegments; ++segmenti)
    {
        const label nDiv = divisions[segmenti];
        const scalar expRatio = expansion[segmenti];

        SubList<scalar> subPoint(gridPoint, nDiv+1, start);
        start += nDiv;

        subPoint[0] = knots[segmenti];
        subPoint[nDiv] = knots[segmenti+1];

        const scalar dist = (subPoint.last() - subPoint.first());

        if
        (
            expandType == expansionType::EXPAND_UNIFORM
         || equal(expRatio, 1)
        )
        {
            for (label i=1; i < nDiv; ++i)
            {
                subPoint[i] = (subPoint[0] + (dist * i)/nDiv);
            }
        }
        else
        {
            // Geometric expansion factor from the expansion ratio
            const scalar expFact = calcGexp(expRatio, nDiv);

            for (label i=1; i < nDiv; ++i)
            {
                subPoint[i] =
                (
                    subPoint[0]
                  + dist * (1.0 - pow(expFact, i))/(1.0 - pow(expFact, nDiv))
                );
            }
        }
    }
}


void Foam::PDRblock::readBoundary(const dictionary& dict)
{
    Log << "Reading boundary entries" << nl;

    PtrList<entry> patchEntries;

    if (dict.found("boundary"))
    {
        PtrList<entry> newEntries(dict.lookup("boundary"));
        patchEntries.transfer(newEntries);
    }

    patches_.clear();
    patches_.resize(patchEntries.size());

    // Hex cell has 6 sides
    const labelRange validFaceId(0, 6);

    // Track which sides have been handled
    labelHashSet handled;

    forAll(patchEntries, patchi)
    {
        if (!patchEntries.set(patchi))
        {
            continue;
        }

        const entry& e = patchEntries[patchi];

        if (!e.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Entry " << e
                << " in boundary section is not a valid dictionary."
                << exit(FatalIOError);
        }

        const dictionary& patchDict = e.dict();

        const word& patchName = e.keyword();

        const word patchType(patchDict.get<word>("type"));

        labelList faceLabels(patchDict.get<labelList>("faces"));

        labelHashSet patchFaces(faceLabels);

        if (faceLabels.size() != patchFaces.size())
        {
            WarningInFunction
                << "Patch: " << patchName
                << ": Ignoring duplicate face ids" << nl;
        }

        label nErased = patchFaces.filterKeys(validFaceId);
        if (nErased)
        {
            WarningInFunction
                << "Patch: " << patchName << ": Ignoring " << nErased
                << " faces with invalid identifiers" << nl;
        }

        nErased = patchFaces.erase(handled);
        if (nErased)
        {
            WarningInFunction
                << "Patch: " << patchName << ": Ignoring " << nErased
                << " faces, which were already handled" << nl;
        }

        // Mark these faces as having been handled
        handled += patchFaces;


        // Save information for later access during mesh creation.

        boundaryEntry& bentry = patches_.emplace_set(patchi);

        bentry.name_ = patchName;
        bentry.type_ = patchType;
        bentry.size_ = 0;
        bentry.faces_ = patchFaces.sortedToc();

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }
    }


    // Deal with unhandled faces

    labelHashSet missed(identity(6));
    missed.erase(handled);

    if (missed.size())
    {
        boundaryEntry& bentry = patches_.emplace_back();

        bentry.name_ = "defaultFaces";
        bentry.type_ = emptyPolyPatch::typeName;
        bentry.size_ = 0;
        bentry.faces_ = missed.sortedToc();

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }


        // Use name/type from "defaultPatch" entry if available
        // - be generous with error handling

        const dictionary* defaultEntry = dict.findDict("defaultPatch");
        if (defaultEntry)
        {
            defaultEntry->readIfPresent("name", bentry.name_);
            defaultEntry->readIfPresent("type", bentry.type_);
        }
    }

    if (verbose_)
    {
        label patchi = 0;
        for (const boundaryEntry& bentry : patches_)
        {
            Info<< "    patch: " << patchi
                << " (size: " << bentry.size_
                << " type: " << bentry.type_
                << ") name: " << bentry.name_
                << " faces: " << flatOutput(bentry.faces_) << nl;

            ++patchi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRblock::PDRblock()
:
    PDRblock(dictionary::null, false)
{}


Foam::PDRblock::PDRblock
(
    const UList<scalar>& xgrid,
    const UList<scalar>& ygrid,
    const UList<scalar>& zgrid
)
:
    PDRblock(dictionary::null, false)
{
    addDefaultPatches();
    reset(xgrid, ygrid, zgrid);
}


Foam::PDRblock::PDRblock(const boundBox& box, const labelVector& nCells)
:
    PDRblock(dictionary::null, false)
{
    addDefaultPatches();
    reset(box, nCells);
}


Foam::PDRblock::PDRblock(const dictionary& dict, bool verboseOutput)
:
    ijkMesh(),
    meshDict_(dict),
    grid_(),
    outer_(),
    bounds_(),
    patches_(),
    edgeLimits_(0,0),
    verbose_(verboseOutput)
{
    if (!dict.isNullDict())
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PDRblock::read(const dictionary& dict)
{
    const scalar scaleFactor(dict.getOrDefault<scalar>("scale", -1));

    expansionType expandType
    (
        expansionNames_.getOrDefault
        (
            "expansion",
            dict,
            expansionType::EXPAND_RATIO
        )
    );

    readGridControl(0, dict.subDict("x"), scaleFactor, expandType);
    readGridControl(1, dict.subDict("y"), scaleFactor, expandType);
    readGridControl(2, dict.subDict("z"), scaleFactor, expandType);

    adjustSizes();

    readBoundary(dict);

    // Outer treatment: (none | extend | box | sphere)
    outer_.clear();

    const dictionary* outerDictPtr = dict.findDict("outer");
    if (outerDictPtr)
    {
        outer_.read(*outerDictPtr);
    }
    outer_.report(Info);

    return true;
}


void Foam::PDRblock::reset
(
    const UList<scalar>& xgrid,
    const UList<scalar>& ygrid,
    const UList<scalar>& zgrid
)
{
    static_cast<scalarList&>(grid_.x()) = xgrid;
    static_cast<scalarList&>(grid_.y()) = ygrid;
    static_cast<scalarList&>(grid_.z()) = zgrid;

    #ifdef FULLDEBUG
    for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
    {
        checkMonotonic(cmpt, grid_[cmpt]);
    }
    #endif

    adjustSizes();

    // Adjust boundaries
    for (boundaryEntry& bentry : patches_)
    {
        bentry.size_ = 0;

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }
    }
}


void Foam::PDRblock::reset(const boundBox& box, const labelVector& nCells)
{
    grid_.x().reset(box.min().x(), box.max().x(), nCells.x());
    grid_.y().reset(box.min().y(), box.max().y(), nCells.y());
    grid_.z().reset(box.min().z(), box.max().z(), nCells.z());

    adjustSizes();

    // Adjust boundaries
    for (boundaryEntry& bentry : patches_)
    {
        bentry.size_ = 0;

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }
    }
}


bool Foam::PDRblock::findCell(const point& pt, labelVector& pos) const
{
    // Out-of-bounds is handled explicitly, for efficiency and consistency,
    // but principally to ensure that findLower() returns a valid
    // result when the point is to the right of the bounds.

    // Since findLower returns the lower index, it corresponds to the
    // cell in which the point is found

    if (!bounds_.contains(pt))
    {
        return false;
    }

    for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
    {
        // Binary search
        pos[cmpt] = findLower(grid_[cmpt], pt[cmpt]);
    }

    return true;
}


bool Foam::PDRblock::gridIndex
(
    const point& pt,
    labelVector& pos,
    const scalar relTol
) const
{
    const scalar tol = relTol * edgeLimits_.min();

    for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
    {
        // Linear search
        pos[cmpt] = grid_[cmpt].findIndex(pt[cmpt], tol);

        if (pos[cmpt] < 0) return false;
    }

    return true;
}


Foam::labelVector Foam::PDRblock::findCell(const point& pt) const
{
    labelVector pos;

    if (findCell(pt, pos))
    {
        return pos;
    }

    return labelVector(-1,-1,-1);
}


Foam::labelVector Foam::PDRblock::gridIndex
(
    const point& pt,
    const scalar relTol
) const
{
    labelVector pos;

    if (gridIndex(pt, pos, relTol))
    {
        return pos;
    }

    return labelVector(-1,-1,-1);
}


Foam::Vector<Foam::gradingDescriptors> Foam::PDRblock::grading() const
{
    return grading(control_);
}


Foam::gradingDescriptors Foam::PDRblock::grading(const direction cmpt) const
{
    switch (cmpt)
    {
        case vector::X :
        case vector::Y :
        case vector::Z :
        {
            return control_[cmpt].grading();
            break;
        }

        default :
            FatalErrorInFunction
                << "Not gridControl for direction " << label(cmpt) << endl
                << exit(FatalError);
            break;
    }

    return gradingDescriptors();
}


// ************************************************************************* //
