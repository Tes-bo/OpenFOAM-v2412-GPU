/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::ensightFaMesh

Description
    Encapsulation of area meshes for writing in ensight format.

Note
    Currently restricted to a single faMesh representation.
    The face elements are created from a specified subset of polyMesh
    faces. The original ordering of these faces is retained in the
    ensightFaces faceOrder().

SourceFiles
    ensightFaMesh.C
    ensightFaMeshI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_ensightMesh_H
#define Foam_fa_ensightMesh_H

#include "ensightFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;
class ensightGeoFile;
class ensightFaMesh;

/*---------------------------------------------------------------------------*\
                         Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightFaMesh
{
    // Private Data

        //- Reference to the finite-area mesh
        const faMesh& mesh_;

        //- Face elements for the area mesh (currently only one)
        ensightFaces areaPart_;

        //- Track if it needs an update
        mutable bool needsUpdate_;

        //- Output verbosity level
        int verbose_;


    // Private Member Functions

        //- Clear all storage
        void clear();

        //- Enforce consistent index/part numbering
        void renumber();


public:

    // Generated Methods

        //- No copy construct
        ensightFaMesh(const ensightFaMesh&) = delete;

        //- No copy assignment
        void operator=(const ensightFaMesh&) = delete;


    // Constructors

        //- Construct from mesh with all default options
        explicit ensightFaMesh(const faMesh& mesh);


    // Member Functions

        //- Output verbosity level
        int verbose() const noexcept;

        //- Change the output verbosity level.
        //  \return old level
        int verbose(const int level) noexcept;


    // Access

        //- Reference to the underlying faMesh
        const faMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Face elements for finite-area
        const ensightFaces& areaPart() const noexcept
        {
            return areaPart_;
        }


    // Other

        //- Does the content need an update?
        bool needsUpdate() const noexcept
        {
            return needsUpdate_;
        }

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if already marked as expired.
        inline bool expire();

        //- Update for new mesh
        void correct();


    // Output

        //- Write geometry to file (normally in parallel).
        //- Adds beginGeometry() marker.
        void write
        (
            ensightGeoFile& os,
            bool parallel = UPstream::parRun()
        ) const;

        //- Write geometry to file (normally in parallel).
        //- Adds beginGeometry() marker.
        FOAM_DEPRECATED_FOR(2024-05, "write(ensightGeoFile&, ...")
        inline void write
        (
            autoPtr<ensightGeoFile>& os,
            bool parallel = UPstream::parRun()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightFaMeshI.H"

#endif

// ************************************************************************* //
