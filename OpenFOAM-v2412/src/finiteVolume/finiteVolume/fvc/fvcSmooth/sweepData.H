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
    Foam::sweepData

Description
    Helper class used by fvc::sweep function.

SourceFiles
    sweepData.H
    sweepDataI.H

\*---------------------------------------------------------------------------*/

#ifndef sweepData_H
#define sweepData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class sweepData;

Ostream& operator<<(Ostream&, const sweepData&);
Istream& operator>>(Istream&, sweepData&);

/*---------------------------------------------------------------------------*\
                       Class sweepData Declaration
\*---------------------------------------------------------------------------*/

class sweepData
{
    // Private Data

        scalar value_;

        point origin_;


    // Private Member Functions

        //- Update gets information from neighbouring face/cell and
        //- uses this to update itself (if necessary) and return true
        template<class TrackingData>
        inline bool update
        (
            const sweepData& svf,
            const point& position,
            const scalar tol,
            TrackingData& td
        );


public:

    // Generated Methods

        //- Copy construct
        sweepData(const sweepData&) = default;

        //- Copy assignment
        sweepData& operator=(const sweepData&) = default;


    // Constructors

        //- Default construct
        inline sweepData();

        //- Construct from components
        inline sweepData(const scalar value, const point& origin);


    // Member Functions

        // Access

            //- Return value
            scalar value() const
            {
                return value_;
            }

            //- The value
            scalar& value()
            {
                return value_;
            }

            //- Return origin
            const point& origin() const
            {
                return origin_;
            }


        // Needed by FaceCellWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const sweepData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to
            //  (patch)face centre
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

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const sweepData&, TrackingData& td) const;


    // Member Operators

        //- Assign new value
        void operator=(const scalar value)
        {
            value_ = value;
        }

        //- Test for equality
        inline bool operator==(const sweepData&) const;

        //- Test for inequality
        inline bool operator!=(const sweepData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const sweepData& rhs)
        {
            return os << rhs.value_ << rhs.origin_;
        }

        friend Istream& operator>>(Istream& is, sweepData& rhs)
        {
            return is >> rhs.value_ >> rhs.origin_;
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for sweepData
template<> struct is_contiguous<sweepData> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sweepDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
