/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::patchEdgeFaceRegions

Description
    Transport of regions for use in PatchEdgeFaceWave.

    Set element to -1 to denote blocked.

SourceFiles
    patchEdgeFaceRegionsI.H
    patchEdgeFaceRegions.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceRegions_H
#define patchEdgeFaceRegions_H

#include "labelList.H"
#include "scalar.H"
#include "tensor.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class patchEdgeFaceRegions;
Istream& operator>>(Istream&, patchEdgeFaceRegions&);
Ostream& operator<<(Ostream&, const patchEdgeFaceRegions&);

/*---------------------------------------------------------------------------*\
                    Class patchEdgeFaceRegions Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceRegions
{
    // Private Data

        //- Region per point
        labelList regions_;

public:

    // Constructors

        //- Default construct
        inline patchEdgeFaceRegions();

        //- Construct from regions
        inline patchEdgeFaceRegions(const labelList& regions);

        //- Construct from regions (on edge)
        inline patchEdgeFaceRegions(const labelPair& regions);


    // Member Functions

        // Access

            const labelList& regions() const
            {
                return regions_;
            }
            labelList& regions()
            {
                return regions_;
            }


        // Needed by MeshWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class Patch, class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const Patch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class Patch, class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const Patch& patch,
                const label edgeI,
                const label facei,
                const patchEdgeFaceRegions& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class Patch, class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const Patch& patch,
                const patchEdgeFaceRegions& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class Patch, class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const Patch& patch,
                const label facei,
                const label edgeI,
                const patchEdgeFaceRegions& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceRegions&, TrackingData&) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const patchEdgeFaceRegions&) const;

        //- Test for inequality
        inline bool operator!=(const patchEdgeFaceRegions&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceRegions&);
        friend Istream& operator>>(Istream&, patchEdgeFaceRegions&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchEdgeFaceRegionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
