/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020,2022 OpenCFD Ltd.
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

Class
    Foam::wallPoints

Description
    For use with FaceCellWave. Determines topological distance to starting faces

SourceFiles
    wallPointsI.H
    wallPoints.C

\*---------------------------------------------------------------------------*/

#ifndef wallPoints_H
#define wallPoints_H

#include "point.H"
#include "tensor.H"
#include "DynamicList.H"
#include "labelList.H"
#include "scalarList.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class wallPoints;
Istream& operator>>(Istream&, wallPoints&);
Ostream& operator<<(Ostream&, const wallPoints&);

/*---------------------------------------------------------------------------*\
                           Class wallPoints Declaration
\*---------------------------------------------------------------------------*/

class wallPoints
{
public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Per face whether the face should not be walked through
        const bitSet& isBlockedFace_;

        //- Per surface, per region the blockSize
        const List<scalarList>& regionToBlockSize_;

        trackData
        (
            const bitSet& isBlockedFace,
            const List<scalarList>& regionToBlockSize
        )
        :
            isBlockedFace_(isBlockedFace),
            regionToBlockSize_(regionToBlockSize)
        {}
    };


protected:

    // Protected Data

        //- Starting points
        DynamicList<point> origin_;

        //- Distance (squared) from cellcenter to origin
        DynamicList<scalar> distSqr_;

        //- Originating surface,region and topological region
        DynamicList<FixedList<label, 3>> surface_;

        //- Originating normal
        //DynamicList<vector> normal_;


    // Protected Member Functions

        //- Evaluate distance to point.
        //  Update distSqr, origin from whomever is nearer pt.
        //  \return true if w2 is closer to point, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const label index1,
            const wallPoints& w2,
            const label index2,

            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Default construct
        inline wallPoints();

        //- Construct from count
        inline wallPoints
        (
            const UList<point>& origin,
            const UList<scalar>& distSqr,
            const UList<FixedList<label, 3>>& surface
            //const UList<vector>& normal
        );


    // Member Functions

        // Access

            const List<point>& origin() const
            {
                return origin_;
            }
            const List<scalar>& distSqr() const
            {
                return distSqr_;
            }

            const List<FixedList<label, 3>>& surface() const
            {
                return surface_;
            }

            //const List<vector>& normal() const
            //{
            //    return normal_;
            //}


        // Needed by FaceCellWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const wallPoints&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const wallPoints& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const wallPoints& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const wallPoints& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const wallPoints&, TrackingData&) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const wallPoints&) const;

        //- Test for inequality
        inline bool operator!=(const wallPoints&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallPoints&);
        friend Istream& operator>>(Istream&, wallPoints&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
