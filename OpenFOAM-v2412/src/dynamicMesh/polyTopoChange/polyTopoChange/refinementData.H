/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020,2024 OpenCFD Ltd.
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
    Foam::refinementData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in FaceCellWave.

SourceFiles
    refinementDataI.H
    refinementData.C

\*---------------------------------------------------------------------------*/

#ifndef refinementData_H
#define refinementData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class refinementData;

Istream& operator>>(Istream&, refinementData&);
Ostream& operator<<(Ostream&, const refinementData&);

/*---------------------------------------------------------------------------*\
                           Class refinementData Declaration
\*---------------------------------------------------------------------------*/

class refinementData
{
    // Private Data

        //- Count which triggers refinement
        label refinementCount_;

        //- Refinement level
        label count_;

public:

    // Constructors

        //- Default construct
        inline refinementData();

        //- Construct from count
        inline refinementData(const label refinementCount, const label count);


    // Member Functions

        // Access

            label refinementCount() const
            {
                return refinementCount_;
            }
            label& refinementCount()
            {
                return refinementCount_;
            }

            label count() const
            {
                return count_;
            }
            label& count()
            {
                return count_;
            }

            bool isRefined() const
            {
                return count_ >= refinementCount_;
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
                const refinementData&,
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
                const refinementData& neighbourInfo,
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
                const refinementData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const refinementData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const refinementData&, TrackingData& td) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const refinementData&) const;

        //- Test for inequality
        inline bool operator!=(const refinementData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refinementData&);
        friend Istream& operator>>(Istream&, refinementData&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- No-op rotational transform for base types
refinementData transform(const tensor&, const refinementData val);

//- Contiguous data for refinementData
template<> struct is_contiguous<refinementData> : std::true_type {};

//- Contiguous label data for refinementData
template<> struct is_contiguous_label<refinementData> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refinementDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
