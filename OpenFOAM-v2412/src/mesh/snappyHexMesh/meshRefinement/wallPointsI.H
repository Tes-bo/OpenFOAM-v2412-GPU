/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallPoints::update
(
    const point& pt,
    const label index1,
    const wallPoints& w2,
    const label index2,

    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 = magSqr(pt - w2.origin_[index2]);

    if (!valid(td))
    {
        // currently not yet set so use any value
        distSqr_[index1] = dist2;
        origin_[index1] = w2.origin_[index2];
        surface_[index1] = w2.surface_[index2];
        //normal_[index1] = w2.normal_[index2];

        return true;
    }

    const scalar diff = distSqr_[index1] - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if
    (
        (diff < SMALL)
     || ((distSqr_[index1] > SMALL) && (diff/distSqr_[index1] < tol))
    )
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_[index1] = dist2;
        origin_[index1] = w2.origin_[index2];
        surface_[index1] = w2.surface_[index2];
        //normal_[index1] = w2.normal_[index2];

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallPoints::wallPoints()
:
    origin_(0),
    distSqr_(0),
    surface_(0)
    //normal_(0)
{}


inline Foam::wallPoints::wallPoints
(
    const UList<point>& origin,
    const UList<scalar>& distSqr,
    const UList<FixedList<label, 3>>& surface
    //const UList<vector>& normal
)
:
    origin_(origin),
    distSqr_(distSqr),
    surface_(surface)
    //normal_(normal)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallPoints::valid(TrackingData& td) const
{
    return origin_.size();
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::wallPoints::sameGeometry
(
    const polyMesh&,
    const wallPoints&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::wallPoints::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{
    for (auto& o : origin_)
    {
        o -= faceCentre;
    }
}


// No geometric data.
template<class TrackingData>
inline void Foam::wallPoints::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData&
)
{
    for (auto& o : origin_)
    {
        o = Foam::transform(rotTensor, o);
    }
}


// No geometric data.
template<class TrackingData>
inline void Foam::wallPoints::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{
    // back to absolute form
    for (auto& o : origin_)
    {
        o += faceCentre;
    }
}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::wallPoints::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const wallPoints& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    const point& cc = mesh.cellCentres()[thisCelli];

    bool hasChanged = false;

    forAll(neighbourInfo.surface_, i)
    {
        const FixedList<label, 3>& nbrSurface = neighbourInfo.surface_[i];

        const scalar d2 = magSqr(cc-neighbourInfo.origin_[i]);

        // Optionally check against surface-based block size
        bool propagate = false;
        if (nbrSurface[0] == labelMax)
        {
            // nbrSurface has special value to ignore regionToBlockSize
            propagate = true;
        }
        else
        {
            // Check distance from nbr origin to cc against max walking distance
            const scalar blockSize =
                td.regionToBlockSize_[nbrSurface[0]][nbrSurface[1]];

            propagate = (d2 < Foam::sqr(3*blockSize));
        }

        if (propagate)
        {
            // Real distance less than max gap distance. Note that it should
            // be at least 2 * blockSize (since gap is two cells across).
            // Should be
            // a little bit more to account for castellated path walk. Bit
            // heuristic - try 3. Should be as low as possible to kill off
            // unnecessary waves asap.

            // Find in my surfaces
            label index = surface_.find(nbrSurface);
            if (index == -1)
            {
                // Append
                origin_.append(neighbourInfo.origin_[i]);
                distSqr_.append(d2);
                surface_.append(nbrSurface);
                //normal_.append(neighbourInfo.normal_[i]);
                hasChanged = true;
            }
            else
            {
                hasChanged =
                    update(cc, index, neighbourInfo, i, tol, td)
                 || hasChanged;
            }
        }
        else
        {
            // Real distance more than gap distance so ignore
            //Pout<< "at cell:" << cc << " ignoring nbr info:"
            //    << neighbourInfo.origin_[i]
            //    << " from surface:" << nbrSurface[0]
            //    << " from region:" << nbrSurface[1]
            //    << " bloxkSize:" << blockSize
            //    << " distance:" << Foam::sqrt(d2)
            //    << endl;
        }
    }

    return hasChanged;
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::wallPoints::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const wallPoints& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From cell to its faces.
    bool hasChanged = false;

    if (!td.isBlockedFace_[thisFacei])
    {
        const point& fc = mesh.faceCentres()[thisFacei];

        forAll(neighbourInfo.surface_, i)
        {
            const FixedList<label, 3>& nbrSurface = neighbourInfo.surface_[i];

            const scalar d2 = magSqr(fc-neighbourInfo.origin_[i]);

            // Optionally check against surface-based block size
            bool propagate = false;
            if (nbrSurface[0] == labelMax)
            {
                // nbrSurface has special value to ignore regionToBlockSize
                propagate = true;
            }
            else
            {
                // Check distance from nbr origin to cc against max walking
                // distance
                const scalar blockSize =
                    td.regionToBlockSize_[nbrSurface[0]][nbrSurface[1]];

                propagate = (d2 < Foam::sqr(3*blockSize));
            }

            if (propagate)
            {
                // Real distance less than max gap distance

                // Find in my surfaces
                label index = surface_.find(nbrSurface);
                if (index == -1)
                {
                    // Append
                    origin_.append(neighbourInfo.origin_[i]);
                    distSqr_.append(d2);
                    surface_.append(nbrSurface);
                    //normal_.append(neighbourInfo.normal_[i]);
                    hasChanged = true;
                }
                else
                {
                    hasChanged =
                        update(fc, index, neighbourInfo, i, tol, td)
                     || hasChanged;
                }
            }
            else
            {
                // Real distance more than gap distance so ignore
                //Pout<< "at face:" << fc << " ignoring nbr info:"
                //    << neighbourInfo.origin_[i]
                //    << " from surface:" << nbrSurface[0]
                //    << " from region:" << nbrSurface[1]
                //    << " bloxkSize:" << blockSize
                //    << " distance:" << Foam::sqrt(d2)
                //    << endl;
            }
        }
    }

    return hasChanged;
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::wallPoints::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const wallPoints& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From face to face (e.g. coupled faces)
    bool hasChanged = false;

    if (!td.isBlockedFace_[thisFacei])
    {
        const point& fc = mesh.faceCentres()[thisFacei];

        forAll(neighbourInfo.surface_, i)
        {
            const FixedList<label, 3>& nbrSurface = neighbourInfo.surface_[i];

            const scalar d2 = magSqr(fc-neighbourInfo.origin_[i]);

            // Optionally check against surface-based block size
            bool propagate = false;
            if (nbrSurface[0] == labelMax)
            {
                // nbrSurface has special value to ignore regionToBlockSize
                propagate = true;
            }
            else
            {
                // Check distance from nbr origin to cc against max walking
                // distance
                const scalar blockSize =
                    td.regionToBlockSize_[nbrSurface[0]][nbrSurface[1]];

                propagate = (d2 < Foam::sqr(3*blockSize));
            }

            if (propagate)
            {
                // Real distance less than max gap distance

                // Find in my surfaces
                const label index = surface_.find(nbrSurface);
                if (index == -1)
                {
                    // Append
                    origin_.append(neighbourInfo.origin_[i]);
                    distSqr_.append(d2);
                    surface_.append(nbrSurface);
                    //normal_.append(neighbourInfo.normal_[i]);
                    hasChanged = true;
                }
                else
                {
                    hasChanged =
                        update(fc, index, neighbourInfo, i, tol, td)
                     || hasChanged;
                }
            }
            else
            {
                // Real distance more than gap distance so ignore
                //Pout<< "at face:" << fc << " ignoring nbr info:"
                //    << neighbourInfo.origin_[i]
                //    << " from surface:" << nbrSurface[0]
                //    << " from region:" << nbrSurface[1]
                //    << " bloxkSize:" << blockSize
                //    << " distance:" << Foam::sqrt(d2)
                //    << endl;
            }
        }
    }

    return hasChanged;
}


template<class TrackingData>
inline bool Foam::wallPoints::equal
(
    const wallPoints& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallPoints::operator==
(
    const wallPoints& rhs
) const
{
    return
        surface_ == rhs.surface_
     && distSqr_ == rhs.distSqr_
     && origin_ == rhs.origin_;
     //&& normal_ == rhs.normal_;
}


inline bool Foam::wallPoints::operator!=
(
    const wallPoints& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
