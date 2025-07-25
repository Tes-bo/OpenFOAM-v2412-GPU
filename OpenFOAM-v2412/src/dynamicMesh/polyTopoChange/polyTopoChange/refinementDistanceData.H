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
    Foam::refinementDistanceData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in FaceCellWave.

SourceFiles
    refinementDistanceDataI.H
    refinementDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef refinementDistanceData_H
#define refinementDistanceData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class refinementDistanceData;

Istream& operator>>(Istream&, refinementDistanceData&);
Ostream& operator<<(Ostream&, const refinementDistanceData&);


/*---------------------------------------------------------------------------*\
                           Class refinementDistanceData Declaration
\*---------------------------------------------------------------------------*/

class refinementDistanceData
{
    // Private Data

        //- Unrefined (level0) buffer size (nBufferLayers*level0Size)
        scalar level0Size_;

        //- Nearest point with highest level
        point origin_;
        label originLevel_;


    // Private Member Functions

        //- Updates with neighbouring data.
        //  \return true if something changed.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const refinementDistanceData& neighbourInfo,
            const scalar tol,
            TrackingData&
        );

public:

    // Constructors

        //- Default construct
        inline refinementDistanceData();

        //- Construct from count
        inline refinementDistanceData
        (
            const scalar level0Size,
            const point& origin,
            const label level
        );


    // Member Functions

        // Access

            scalar level0Size() const
            {
                return level0Size_;
            }
            scalar& level0Size()
            {
                return level0Size_;
            }

            const point& origin() const
            {
                return origin_;
            }
            point& origin()
            {
                return origin_;
            }

            label originLevel() const
            {
                return originLevel_;
            }
            label& originLevel()
            {
                return originLevel_;
            }


        // Other

            //- Calculates the wanted level at a given point.
            //  Walks out from the origin.
            inline label wantedLevel(const point& pt) const;


        // Needed by FaceCellWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData&) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const refinementDistanceData&,
                const scalar,
                TrackingData&
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
                TrackingData&
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData&
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData&
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const refinementDistanceData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const refinementDistanceData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const refinementDistanceData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const refinementDistanceData&, TrackingData&)
            const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const refinementDistanceData&) const;

        //- Test for inequality
        inline bool operator!=(const refinementDistanceData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refinementDistanceData&);
        friend Istream& operator>>(Istream&, refinementDistanceData&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for refinementDistanceData
template<> struct is_contiguous<refinementDistanceData> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refinementDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
