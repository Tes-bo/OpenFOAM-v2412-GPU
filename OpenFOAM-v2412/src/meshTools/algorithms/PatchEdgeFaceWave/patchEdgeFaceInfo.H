/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::patchEdgeFaceInfo

Description

SourceFiles
    patchEdgeFaceInfoI.H
    patchEdgeFaceInfo.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceInfo_H
#define patchEdgeFaceInfo_H

#include "label.H"
#include "scalar.H"
#include "point.H"
#include "tensor.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class patchEdgeFaceInfo;

Istream& operator>>(Istream&, patchEdgeFaceInfo&);
Ostream& operator<<(Ostream&, const patchEdgeFaceInfo&);

/*---------------------------------------------------------------------------*\
                           Class patchEdgeFaceInfo Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceInfo
{
    // Private Data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from point to origin
        scalar distSqr_;


    // Private Member Functions

        //- Evaluate distance to point.
        //  Update distSqr, origin from whomever is nearer pt.
        //  \return true if w2 is closer to point, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const patchEdgeFaceInfo& w2,
            const scalar tol,
            TrackingData& td
        );

        //- Combine current with w2. Update distSqr, origin if w2 has smaller
        //  quantities and returns true.
        template<class TrackingData>
        inline bool update
        (
            const patchEdgeFaceInfo& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Default construct
        inline patchEdgeFaceInfo();

        //- Construct from origin, distance squared
        inline patchEdgeFaceInfo(const point& origin, const scalar distSqr);


    // Member Functions

        // Access

            const point& origin() const
            {
                return origin_;
            }

            scalar distSqr() const
            {
                return distSqr_;
            }


        // Needed by MeshWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label edgeI,
                const label facei,
                const patchEdgeFaceInfo& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const patchEdgeFaceInfo& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label facei,
                const label edgeI,
                const patchEdgeFaceInfo& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceInfo&, TrackingData& td) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const patchEdgeFaceInfo&) const;

        //- Test for inequality
        inline bool operator!=(const patchEdgeFaceInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceInfo&);
        friend Istream& operator>>(Istream&, patchEdgeFaceInfo&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for patchEdgeFaceInfo
template<> struct is_contiguous<patchEdgeFaceInfo> : std::true_type {};

//- Contiguous scalar data for patchEdgeFaceInfo
template<> struct is_contiguous_scalar<patchEdgeFaceInfo> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFaceInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
