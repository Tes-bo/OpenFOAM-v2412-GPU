/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::directionalWallPointData

Description
    Holds information (coordinate and normal) regarding the nearest wall point.

    Variant of \c wallPointData that ignores the specified normal component
    before comparing. This is used e.g. to find the distance to the wall
    in the z-direction only.

SourceFiles
    directionalWallPointDataI.H
    directionalWallPointData.C

\*---------------------------------------------------------------------------*/

#ifndef directionalWallPointData_H
#define directionalWallPointData_H

#include "wallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class directionalWallPointData;

// Forward declaration of friend functions and operators

template<class Type> Istream&
operator>>(Istream&, directionalWallPointData<Type>&);
template<class Type> Ostream&
operator<<(Ostream&, const directionalWallPointData<Type>&);


/*---------------------------------------------------------------------------*\
                  Class directionalWallPointData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class directionalWallPointData
:
    public wallPointData<Type>
{
    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const directionalWallPointData<Type>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    typedef Type dataType;


    // Constructors

        //- Construct null
        inline directionalWallPointData();

        //- Construct from origin, normal, distance
        inline directionalWallPointData
        (
            const point& origin,
            const Type& data,
            const scalar distSqr
        );


    // Member Functions

        // Needed by meshWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh& mesh,
                const label thisCelli,
                const label neighbourFacei,
                const directionalWallPointData<Type>& neighbourWallInfo,
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
                const directionalWallPointData<Type>& neighbourWallInfo,
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
                const directionalWallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

    // Member Operators

    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const directionalWallPointData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            directionalWallPointData<Type>&
        );
};


//  Data associated with directionalWallPointData type are contiguous. List the
//  usual ones.

template<> struct is_contiguous<directionalWallPointData<bool>> :
is_contiguous<wallPoint> {};

template<> struct is_contiguous<directionalWallPointData<label>> :
is_contiguous_label<wallPoint> {};

template<> struct is_contiguous<directionalWallPointData<scalar>> :
is_contiguous_scalar<wallPoint> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "directionalWallPointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "directionalWallPointDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
