/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::tetIndices

Description
    Storage and named access for the indices of a tet which is part of
    the decomposition of a cell.

    Tets are designated by
    - cell (of course)
    - face on cell
    - three points on face (faceBasePt, facePtA, facePtB)
      When constructing from a mesh and index in the face (tetPtI):
        - faceBasePt is the mesh.tetBasePtIs() base point
        - facePtA is tetPtI away from faceBasePt
        - facePtB is next one after/before facePtA
        e.g.:

            +---+
            |2 /|
            | / |
            |/ 1|  <- tetPt (so 1 for first triangle, 2 for second)
            +---+
            ^
           faceBasePt

SourceFiles
    tetIndicesI.H
    tetIndices.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tetIndices_H
#define Foam_tetIndices_H

#include "label.H"
#include "tetrahedron.H"
#include "triangle.H"
#include "polyMesh.H"
#include "triFace.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class tetIndices;

Istream& operator>>(Istream&, tetIndices&);
Ostream& operator<<(Ostream&, const tetIndices&);

/*---------------------------------------------------------------------------*\
                         Class tetIndices Declaration
\*---------------------------------------------------------------------------*/

class tetIndices
{
    // Private Data

        //- Cell that this is a decomposed tet of
        label celli_;

        //- Face that holds this decomposed tet
        label facei_;

        //- Point on the face, *relative to the base point*,
        //- which characterises this tet on the face.
        label tetPti_;


    // Private Static Data

        //- Maximum number of bad tet warnings
        static constexpr int maxNWarnings = 100;

        //- Current number of bad tet warnings.
        //- Warnings stop when this reaches the maximum number.
        static int nWarnings_;


public:

    // Constructors

        //- Default construct, with invalid labels (-1)
        inline constexpr tetIndices() noexcept;

        //- Construct from components
        inline constexpr tetIndices
        (
            label celli,
            label facei,
            label tetPointi
        ) noexcept;


    //- Destructor
    ~tetIndices() = default;


    // Member Functions

    // Access

        //- Return the cell index
        label cell() const noexcept { return celli_; }

        //- Non-const access to the cell index
        label& cell() noexcept { return celli_; }

        //- Return the face index
        label face() const noexcept { return facei_; }

        //- Non-const access to the face index
        label& face() noexcept { return facei_; }

        //- Return the characterising tet point index
        label tetPt() const noexcept { return tetPti_; }

        //- Non-const access to the characterising tet point index
        label& tetPt() noexcept { return tetPti_; }


    // Mesh-related

        //- The indices corresponding to the tri on the face for this tet.
        //- The normal of the tri points out of the cell.
        inline triFace faceTriIs
        (
            const polyMesh& mesh,
            const bool warn = true
        ) const;

        //- Local indices corresponding to the tri on the face for this tet.
        //- The normal of the tri points out of the cell.
        inline triFace triIs
        (
            const polyMesh& mesh,
            const bool warn = true
        ) const;


        //- The tet geometry for this tet,
        //- where point0 is the cell centre.
        inline tetPointRef tet(const polyMesh& mesh) const;

        //- The tet geometry for this tet (using old positions),
        //- where point0 is the cell centre.
        inline tetPointRef oldTet(const polyMesh& mesh) const;

        //- The triangle geometry for the face for this tet.
        //- The normal of the tri points out of the cell
        inline triPointRef faceTri(const polyMesh& mesh) const;

        //- The triangle geometry for the face for this tet
        //- (using old positions)
        inline triPointRef oldFaceTri(const polyMesh& mesh) const;

        //- The x/y/z position for given barycentric coordinates
        //- (where point0 is the cell centre).
        inline point barycentricToPoint
        (
            const polyMesh& mesh,
            const barycentric& bary
        ) const;


    // Other

        //- Compare tetIndices for equality.
        //- Compares cell, face, tetPt elements in order, stopping at the
        //- first inequality.
        //
        //  \returns negative/zero/positive from the last element compared
        static inline label compare
        (
            const tetIndices& a,
            const tetIndices& b
        ) noexcept;


    // IOstream Operators

        friend Istream& operator>>(Istream&, tetIndices&);
        friend Ostream& operator<<(Ostream&, const tetIndices&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for tetIndices
template<> struct is_contiguous<tetIndices> : std::true_type {};

//- Contiguous label data for tetIndices
template<> struct is_contiguous_label<tetIndices> : std::true_type {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const tetIndices& a, const tetIndices& b) noexcept;
inline bool operator!=(const tetIndices& a, const tetIndices& b) noexcept;


inline bool operator<(const tetIndices& a, const tetIndices& b) noexcept
{
    return (tetIndices::compare(a, b) < 0);
}

inline bool operator<=(const tetIndices& a, const tetIndices& b) noexcept
{
    return (tetIndices::compare(a, b) <= 0);
}

inline bool operator>(const tetIndices& a, const tetIndices& b) noexcept
{
    return (tetIndices::compare(a, b) > 0);
}

inline bool operator>=(const tetIndices& a, const tetIndices& b) noexcept
{
    return (tetIndices::compare(a, b) >= 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetIndicesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
