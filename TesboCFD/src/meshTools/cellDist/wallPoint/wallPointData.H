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
    Foam::wallPointData

Description
    Holds information (coordinate and normal) regarding nearest wall point.

    Is like wallPoint but transfer extra (passive) data.
    Used e.g. in wall distance calculation with wall reflection vectors.

SourceFiles
    wallPointDataI.H
    wallPointData.C

\*---------------------------------------------------------------------------*/

#ifndef wallPointData_H
#define wallPointData_H

#include "wallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class wallPointData;

template<class Type> Istream& operator>>(Istream&, wallPointData<Type>&);
template<class Type> Ostream& operator<<(Ostream&, const wallPointData<Type>&);


/*---------------------------------------------------------------------------*\
                        Class wallPointData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wallPointData
:
    public wallPoint
{
    // Private Data

        //- Data at nearest wall center
        Type data_;


    // Private Member Functions

        //- Evaluate distance to point.
        //  Update distSqr, origin from whomever is nearer pt.
        //  \return true if w2 is closer to point, false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const wallPointData<Type>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Public Typedefs

        //- Type of additional data carried
        typedef Type dataType;


    // Constructors

        //- Default construct
        inline wallPointData();

        //- Construct from origin, normal, distance
        inline wallPointData
        (
            const point& origin,
            const Type& data,
            const scalar distSqr
        );


    // Member Functions

        // Access

            const Type& data() const
            {
                return data_;
            }
            Type& data()
            {
                return data_;
            }


        // Needed by MeshWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCelli,
                const label neighbourFacei,
                const wallPointData<Type>& neighbourWallInfo,
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
                const wallPointData<Type>& neighbourWallInfo,
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
                const wallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <Type>(Ostream&, const wallPointData<Type>&);
        friend Istream& operator>> <Type>(Istream&, wallPointData<Type>&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if data type is contiguous
template<class Type>
struct is_contiguous<wallPointData<Type>> : is_contiguous<Type> {};

//- Data are contiguous label if data type is label
template<class Type>
struct is_contiguous_label<wallPointData<Type>> : is_contiguous_label<Type> {};

//- Data are contiguous scalar if data type is scalar
template<class Type>
struct is_contiguous_scalar<wallPointData<Type>> : is_contiguous_scalar<Type>{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wallPointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
