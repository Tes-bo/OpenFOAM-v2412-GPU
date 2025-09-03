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
    Foam::topoDistanceData

Description
    For use with FaceCellWave. Determines topological distance to starting
    faces. Templated on passive transported data.

SourceFiles
    topoDistanceDataI.H
    topoDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef topoDistanceData_H
#define topoDistanceData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
template<class Type> class topoDistanceData;

template<class Type>
Istream& operator>>(Istream&, topoDistanceData<Type>&);
template<class Type>
Ostream& operator<<(Ostream&, const topoDistanceData<Type>&);

/*---------------------------------------------------------------------------*\
                      Class topoDistanceData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class topoDistanceData
{
    // Private Data

        //- Distance
        label distance_;

        //- Starting data
        Type data_;


public:

    // Constructors

        //- Construct null with invalid (-1) for distance, null constructor
        //  for data
        inline topoDistanceData();

        //- Construct from components
        inline topoDistanceData
        (
            const label distance,
            const Type& data
        );


    // Member Functions

        // Access

            inline label distance() const
            {
                return distance_;
            }

            inline const Type& data() const
            {
                return data_;
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
                const topoDistanceData<Type>&,
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
                const topoDistanceData<Type>& neighbourInfo,
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
                const topoDistanceData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const topoDistanceData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal
            (
                const topoDistanceData<Type>&,
                TrackingData& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const topoDistanceData<Type>&) const;

        inline bool operator!=(const topoDistanceData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const topoDistanceData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            topoDistanceData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if data type is contiguous
template<class Type>
struct is_contiguous<topoDistanceData<Type>> : is_contiguous<Type> {};

//- Data are contiguous label if data type is label
template<class Type>
struct is_contiguous_label<topoDistanceData<Type>> :
    is_contiguous_label<Type> {};

//- Data are contiguous scalar if data type is scalar
template<class Type>
struct is_contiguous_scalar<topoDistanceData<Type>> :
    is_contiguous_scalar<Type>{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "topoDistanceData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "topoDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
