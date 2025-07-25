/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::ijkAddressing

Description
    A simple i-j-k (row-major order) to linear addressing.

SourceFiles
    ijkAddressingI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ijkAddressing_H
#define Foam_ijkAddressing_H

#include "labelVector.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ijkAddressing Declaration
\*---------------------------------------------------------------------------*/

class ijkAddressing
{
    // Private Data

        //- The number of items in the i,j,k directions.
        labelVector sizes_;


public:

    // Constructors

        //- Construct zero-size addressing
        inline ijkAddressing();

        //- Construct with addressing
        inline explicit ijkAddressing(const labelVector& ijk);

        //- Construct with addressing components
        inline ijkAddressing(const label ni, const label nj, const label nk);


    // Member Functions

    // Access

        //- Addressing is considered empty if any component is zero
        inline bool empty() const;

        //- The (i,j,k) addressing dimensions
        inline const labelVector& sizes() const;

        //- Return the (i,j,k) dimensions for modification
        inline labelVector& sizes();

        //- Return the total i*j*k size
        inline label size() const;

        //- The addressing dimension in the given direction
        inline const label& size(const vector::components cmpt) const;

        //- Reset to (0,0,0) sizing
        inline void clear();

        //- Change the sizing parameters
        inline void reset(const label ni, const label nj, const label nk);

        //- Change the sizing parameters
        inline void reset(const labelVector& newSizes);

        //- Linear addressing index (offset) for an (i,j,k) position.
        inline label index(const label i, const label j, const label k) const;

        //- Linear addressing index (offset) for an (i,j,k) position.
        inline label index(const labelVector& ijk) const;

        //- The (i,j,k) indexing from linear addressing.
        inline labelVector index(const label idx) const;


    // Checks

        //- Check indices are within ni,nj,nk range.
        //  Optionally allow an extra index for point addressing
        inline void checkIndex
        (
            const label i,
            const label j,
            const label k,
            const bool allowExtra = false
        ) const;

        //- Check indices are within ni,nj,nk range.
        //  Optionally allow an extra index for point addressing
        inline void checkIndex
        (
            const labelVector& ijk,
            const bool allowExtra = false
        ) const;

        //- Check that all components of sizes() are non-negative
        inline void checkSizes() const;

        //- Check that all components of sizes() match
        inline void checkSizes(const labelVector& other) const;

        //- Check that the total size matches
        inline void checkSizes(const label nTotal) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ijkAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
