/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::hexMatcher

Description
    A cellMatcher for hex cells (cellModel::HEX).

See also
    cellMatcher, cellModel

SourceFiles
    hexMatcher.C

\*---------------------------------------------------------------------------*/

#ifndef hexMatcher_H
#define hexMatcher_H

#include "cellMatcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hexMatcher Declaration
\*---------------------------------------------------------------------------*/

class hexMatcher
:
    public cellMatcher
{
    // Static Data Members

        // Constants for this shape
        static constexpr label vertPerCell = 8;
        static constexpr label facePerCell = 6;
        static constexpr label maxVertPerFace = 4;


    // Private Member Functions

        //- No copy construct
        hexMatcher(const hexMatcher&) = delete;

        //- No copy assignment
        hexMatcher& operator=(const hexMatcher&) = delete;


public:

    // Constructors

        //- Default construct
        hexMatcher();


    //- Destructor
    ~hexMatcher() = default;


    // Static Functions

        //- Test if given list of faces satisfies criteria for HEX.
        //- (6 quad)
        static bool test(const UList<face>& faces);

        //- Test if given cell satisfies criteria for HEX.
        //- (6 quad)
        static bool test(const primitiveMesh& mesh, const label celli);


    // Member Functions

        virtual label nVertPerCell() const
        {
            return vertPerCell;
        }

        virtual label nFacePerCell() const
        {
            return facePerCell;
        }

        virtual label nMaxVertPerFace() const
        {
            return maxVertPerFace;
        }

        virtual label faceHashValue() const;

        virtual bool faceSizeMatch(const faceList&, const labelList&) const;

        virtual bool matchShape
        (
            const bool checkOnly,
            const faceList& faces,
            const labelList& faceOwner,
            const label celli,
            const labelList& myFaces
        );

        virtual bool isA(const faceList& faces)
        {
            return hexMatcher::test(faces);
        }

        virtual bool isA(const primitiveMesh& mesh, const label celli)
        {
            return hexMatcher::test(mesh, celli);
        }

        virtual bool matches
        (
            const primitiveMesh& mesh,
            const label celli,
            cellShape& shape
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
