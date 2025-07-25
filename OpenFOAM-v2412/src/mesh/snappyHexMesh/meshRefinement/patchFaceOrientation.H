/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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
    Foam::patchFaceOrientation

Description
    Transport of orientation for use in PatchEdgeFaceWave.

SourceFiles
    patchFaceOrientationI.H
    patchFaceOrientation.C

\*---------------------------------------------------------------------------*/

#ifndef patchFaceOrientation_H
#define patchFaceOrientation_H

#include "tensor.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class patchFaceOrientation;

Istream& operator>>(Istream&, patchFaceOrientation&);
Ostream& operator<<(Ostream&, const patchFaceOrientation&);

/*---------------------------------------------------------------------------*\
                    Class patchFaceOrientation Declaration
\*---------------------------------------------------------------------------*/

class patchFaceOrientation
{
    // Private Data

        //- Flip status
        label flipStatus_;

public:

    // Constructors

        //- Default construct
        inline patchFaceOrientation();

        //- Construct from components
        inline patchFaceOrientation(const label);


    // Member Functions

    // Access

        //- Orientation
        label flipStatus() const
        {
            return flipStatus_;
        }

        //- Reverse the orientation
        inline void flip();


        // Needed by MeshWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label edgeI,
                const label facei,
                const patchFaceOrientation& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const patchFaceOrientation& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label facei,
                const label edgeI,
                const patchFaceOrientation& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const patchFaceOrientation&, TrackingData&) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const patchFaceOrientation&) const;

        //- Test for equality
        inline bool operator!=(const patchFaceOrientation&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchFaceOrientation&);
        friend Istream& operator>>(Istream&, patchFaceOrientation&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for patchFaceOrientation
template<> struct is_contiguous<patchFaceOrientation> : std::true_type {};

//- Contiguous label data for patchFaceOrientation
template<> struct is_contiguous_label<patchFaceOrientation> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchFaceOrientationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
