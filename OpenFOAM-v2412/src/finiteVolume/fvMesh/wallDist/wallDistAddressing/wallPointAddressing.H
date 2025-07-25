/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::wallPointAddressing

Description
    Holds information (coordinate and origin) regarding nearest wall point.

    Is like wallPoint but transfer extra (passive) data.
    Used e.g. in wall distance calculation with wall reflection vectors.

SourceFiles
    wallPointAddressingI.H
    wallPointAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef wallPointAddressing_H
#define wallPointAddressing_H

#include "wallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalIndexAndTransform;
class wallPointAddressing;

Ostream& operator<<(Ostream&, const wallPointAddressing&);
Istream& operator>>(Istream&, wallPointAddressing&);


/*---------------------------------------------------------------------------*\
                     Class wallPointAddressing Declaration
\*---------------------------------------------------------------------------*/

class wallPointAddressing
:
    public wallPoint
{
    // Private Data

        //- Data at nearest wall center
        labelPair data_;


    // Private Member Functions

        //- Initialise data so it uses the nulltransformindex
        inline void setNull(const globalIndexAndTransform& gt);

        //- Evaluate distance to point.
        //  Update distSqr, origin from whomever is nearer pt.
        //  \return true if w2 is closer to point, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const wallPointAddressing& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Default construct
        inline wallPointAddressing();

        //- Construct from origin, transforms, distance
        inline wallPointAddressing
        (
            const point& origin,
            const labelPair& data,
            const scalar distSqr
        );


    // Member Functions

        // Access

            const labelPair& data() const
            {
                return data_;
            }
            labelPair& data()
            {
                return data_;
            }


        // Needed by MeshWave

            //- Convert any absolute coordinates into relative to (patch)face
            //- centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCelli,
                const label neighbourFacei,
                const wallPointAddressing& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFacei,
                const label neighbourCelli,
                const wallPointAddressing& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const label thisFacei,
                const wallPointAddressing& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallPointAddressing&);
        friend Istream& operator>>(Istream&, wallPointAddressing&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if data type is contiguous
template<>
struct is_contiguous<wallPointAddressing>
: is_contiguous<wallPoint>
{};

template<>
struct is_contiguous_scalar<wallPointAddressing>
: is_contiguous_scalar<wallPoint>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
