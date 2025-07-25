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
    Foam::referredWallFace

Description
    Storage for referred wall faces.  Stores patch index, face and
    associated points

SourceFiles
    referredWallFaceI.H
    referredWallFace.C

\*---------------------------------------------------------------------------*/

#ifndef referredWallFace_H
#define referredWallFace_H

#include "face.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class referredWallFace;
Istream& operator>>(Istream&, referredWallFace&);
Ostream& operator<<(Ostream&, const referredWallFace&);

/*---------------------------------------------------------------------------*\
                        Class referredWallFace Declaration
\*---------------------------------------------------------------------------*/

class referredWallFace
:
    public face
{
    // Private Data

        //- Points of face
        pointField pts_;

        //- Index of originating patch
        label patchi_;


public:

    // Generated Methods

        //- Default construct
        referredWallFace() = default;


    // Constructors

        //- Construct from components
        referredWallFace
        (
            const face& f,
            const pointField& pts,
            label patchi
        );


    // Member Functions

        // Access

            //- Return access to the stored points
            inline const pointField& points() const;

            //- Return non-const access to the stored points
            inline pointField& points();

            //- Return access to the patch index
            inline label patchIndex() const;

            //- Return non-const access to the patch index
            inline label& patchIndex();


    // Member Operators

        bool operator==(const referredWallFace&) const;
        bool operator!=(const referredWallFace&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, referredWallFace&);
        friend Ostream& operator<<(Ostream&, const referredWallFace&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredWallFaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
