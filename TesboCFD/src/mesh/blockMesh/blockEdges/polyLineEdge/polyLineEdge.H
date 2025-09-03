/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::blockEdges::polyLineEdge

Description
    A blockEdge defined in terms of a series of straight line segments.

SourceFiles
    polyLineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_polyLineEdge_H
#define blockEdges_polyLineEdge_H

#include "blockEdge.H"
#include "polyLine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                        Class polyLineEdge Declaration
\*---------------------------------------------------------------------------*/

class polyLineEdge
:
    public blockEdge,
    public polyLine
{
    // Private Member Functions

        //- No copy construct
        polyLineEdge(const polyLineEdge&) = delete;

        //- No copy assignment
        void operator=(const polyLineEdge&) = delete;


public:

    // Static data members
    TypeName("polyLine");


    // Constructors

        //- Construct from components
        polyLineEdge
        (
            const pointField& points,       //!< Referenced point field
            const edge& fromTo,             //!< Start/end in point field
            const pointField& intermediate  //!< The intermediate points
        );

        //- Construct from components
        polyLineEdge
        (
            const pointField& points,       //!< Referenced point field
            const label from,               //!< Start point in point field
            const label to,                 //!< End point in point field
            const pointField& intermediate  //!< The intermediate points
        );

        //- Construct from Istream and point field.
        polyLineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces&  /*unused*/,
            const pointField& points,  //!< Referenced point field
            Istream& is
        );


    //- Destructor
    virtual ~polyLineEdge() = default;


    // Member Functions

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar lambda) const;

        //- The length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
