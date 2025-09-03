/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::cyclicPeriodicAMIPolyPatch

Description
    Cyclic patch for periodic Arbitrary Mesh Interface (AMI)

    Uses automated matching to fill all of the faces on owner and neighbour.
    Can be used for translational periodic matching as well. Uses the
    transforms from a specified periodic patch to transform the geometry.

SourceFiles
    cyclicPeriodicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIPolyPatch_H
#define cyclicPeriodicAMIPolyPatch_H

#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OBJstream;

/*---------------------------------------------------------------------------*\
                    Class cyclicPeriodicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMIPolyPatch
:
    public cyclicAMIPolyPatch
{
private:

    // Private data

        //- Current number of transformations (+ve forward, -ve backward).
        //- Used internally to stores the last state that gave a starting
        //- match. Speeds up the automatic matching.
        mutable label nTransforms_;

        //- Number of sectors in a rotationally periodic geometry (optional)
        const label nSectors_;

        //- Maximum number of attempts to match the AMI geometry
        const label maxIter_;


    // Private Member Functions

        //- Synchronise the periodic transformations
        void syncTransforms() const;

        //- Debug: write obj files of patch (collected on master)
        void writeOBJ(const primitivePatch& p, OBJstream& str) const;

        //- Reset the AMI interpolator
        virtual void resetAMI() const;


public:

    //- Runtime type information
    TypeName("cyclicPeriodicAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicPeriodicAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        cyclicPeriodicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPeriodicAMIPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPeriodicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPeriodicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicPeriodicAMIPolyPatch();


    // Member Functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
