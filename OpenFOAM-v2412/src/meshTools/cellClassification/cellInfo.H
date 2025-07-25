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
    Foam::cellInfo

Description
    Holds information regarding type of cell. Used in inside/outside
    determination in cellClassification.

SourceFiles
    cellInfoI.H
    cellInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cellInfo_H
#define cellInfo_H

#include "point.H"
#include "label.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class cellInfo;

Istream& operator>>(Istream&, cellInfo&);
Ostream& operator<<(Ostream&, const cellInfo&);


/*---------------------------------------------------------------------------*\
                           Class cellInfo Declaration
\*---------------------------------------------------------------------------*/

class cellInfo
{
    // Private Data

        label type_;


    // Private Member Functions

        //- Update current cell/face type with neighbouring type.
        //  \return true if information needs to propagate, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const cellInfo& w2,
            const label thisFacei,
            const label thisCelli,
            const label neighbourFacei,
            const label neighbourCelli,
            TrackingData& td
        );


public:

    // Constructors

        //- Default construct - as cellClassification::NOTSET
        inline cellInfo();

        //- Construct from cellClassification type
        inline explicit cellInfo(const label ctype);


    // Member Functions

        // Access

            label type() const
            {
                return type_;
            }
            label& type()
            {
                return type_;
            }


        // Needed by MeshWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const cellInfo&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch& patch,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch& patch,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor& rotTensor,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const cellInfo& neighbourInfo,
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
                const cellInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const cellInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const cellInfo&, TrackingData& td) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const cellInfo&) const;

        //- Test for inequality
        inline bool operator!=(const cellInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const cellInfo&);
        friend Istream& operator>>(Istream&, cellInfo&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for cellInfo
template<> struct is_contiguous<cellInfo> : std::true_type {};

//- Contiguous label data for cellInfo
template<> struct is_contiguous_label<cellInfo> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
