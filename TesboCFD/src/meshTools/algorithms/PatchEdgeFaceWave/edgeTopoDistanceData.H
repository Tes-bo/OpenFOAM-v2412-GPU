/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019,2020 OpenCFD Ltd.
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
    Foam::edgeTopoDistanceData

Description
    For use with PatchEdgeFaceWave. Determines topological distance to
    starting edges. Templated on passive transported data.

SourceFiles
    edgeTopoDistanceDataI.H
    edgeTopoDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef edgeTopoDistanceData_H
#define edgeTopoDistanceData_H

#include "point.H"
#include "tensor.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
template<class Type, class PrimitivePatchType>
class edgeTopoDistanceData;

template<class Type, class PrimitivePatchType>
Istream& operator>>
(
    Istream&,
    edgeTopoDistanceData<Type, PrimitivePatchType>&
);
template<class Type, class PrimitivePatchType>
Ostream& operator<<
(
    Ostream&,
    const edgeTopoDistanceData<Type, PrimitivePatchType>&
);


/*---------------------------------------------------------------------------*\
                    Class edgeTopoDistanceData Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class PrimitivePatchType = indirectPrimitivePatch>
class edgeTopoDistanceData
{
protected:

    // Protected data

        //- Distance
        label distance_;

        //- Starting data
        Type data_;


public:

    typedef Type dataType;

    // Constructors

        //- Construct null with invalid (-1) for distance, null constructor
        //  for data
        inline edgeTopoDistanceData();

        //- Construct from distance, data
        inline edgeTopoDistanceData
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


        // Needed by PatchEdgeFaceWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const PrimitivePatchType& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const PrimitivePatchType& patch,
                const label edgeI,
                const label facei,
                const edgeTopoDistanceData<Type, PrimitivePatchType>& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const PrimitivePatchType& patch,
                const edgeTopoDistanceData<Type, PrimitivePatchType>& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const PrimitivePatchType& patch,
                const label facei,
                const label edgeI,
                const edgeTopoDistanceData<Type, PrimitivePatchType>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const edgeTopoDistanceData<Type, PrimitivePatchType>&,
                TrackingData&
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==
        (
            const edgeTopoDistanceData<Type, PrimitivePatchType>&
        ) const;
        inline bool operator!=
        (
            const edgeTopoDistanceData<Type, PrimitivePatchType>&
        ) const;


    // IOstream Operators

        friend Ostream& operator<< <Type, PrimitivePatchType>
        (
            Ostream&,
            const edgeTopoDistanceData<Type, PrimitivePatchType>&
        );
        friend Istream& operator>> <Type, PrimitivePatchType>
        (
            Istream&,
            edgeTopoDistanceData<Type, PrimitivePatchType>&
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //


//- Data are contiguous if data type is contiguous
template<class Type, class PrimitivePatchType>
struct is_contiguous<edgeTopoDistanceData<Type, PrimitivePatchType>> :
    is_contiguous<Type> {};

//- Data are contiguous label if data type is label
template<class Type, class PrimitivePatchType>
struct is_contiguous_label<edgeTopoDistanceData<Type, PrimitivePatchType>> :
    is_contiguous_label<Type> {};

//- Data are contiguous scalar if data type is scalar
template<class Type, class PrimitivePatchType>
struct is_contiguous_scalar<edgeTopoDistanceData<Type, PrimitivePatchType>> :
    is_contiguous_scalar<Type>{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "edgeTopoDistanceData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeTopoDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
