/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::wallNormalInfo

Description
    Holds information regarding nearest wall point.
    Used in wall refinement.

SourceFiles
    wallNormalInfoI.H
    wallNormalInfo.C

\*---------------------------------------------------------------------------*/

#ifndef wallNormalInfo_H
#define wallNormalInfo_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
class wallNormalInfo;

Istream& operator>>(Istream&, wallNormalInfo&);
Ostream& operator<<(Ostream&, const wallNormalInfo&);

/*---------------------------------------------------------------------------*\
                           Class wallNormalInfo Declaration
\*---------------------------------------------------------------------------*/

class wallNormalInfo
{
    // Private Data

        //- Normal at nearest wall point
        vector normal_;


    // Private Member Functions

        //- Evaluate distance to point and update normal_
        template<class TrackingData>
        inline bool update(const wallNormalInfo& w2, TrackingData& td);


public:

    // Constructors

        //- Default construct
        inline wallNormalInfo();

        //- Construct from normal
        inline wallNormalInfo(const vector& normal);



    // Member Functions

        // Access

            const vector& normal() const
            {
                return normal_;
            }
            vector& normal()
            {
                return normal_;
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
                const wallNormalInfo&,
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
                const wallNormalInfo& neighbourInfo,
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
                const wallNormalInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const wallNormalInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const wallNormalInfo&, TrackingData& td) const;


    // Member Operators

        //- Test for equality
        inline bool operator==(const wallNormalInfo&) const;

        //- Test for inequality
        inline bool operator!=(const wallNormalInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallNormalInfo&);
        friend Istream& operator>>(Istream&, wallNormalInfo&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for wallNormalInfo
template<> struct is_contiguous<wallNormalInfo> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallNormalInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
