/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::blockEdges::projectEdge

Description
    Defines the edge from the projection onto a surface (single surface)
    or intersection of two surfaces.

SourceFiles
    projectEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_projectEdge_H
#define blockEdges_projectEdge_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointConstraint;

namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                         Class projectEdge Declaration
\*---------------------------------------------------------------------------*/

class projectEdge
:
    public blockEdge
{
    // Private Data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;


    // Private Member Functions

        //- Single point find nearest
        void findNearest(const point&, point& near, pointConstraint&) const;

        //- No copy construct
        projectEdge(const projectEdge&) = delete;

        //- No copy assignment
        void operator=(const projectEdge&) = delete;


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream and point field.
        projectEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,   //!< Referenced point field
            Istream& is
        );


    //- Destructor
    virtual ~projectEdge() = default;


    // Member Functions

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- The point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual tmp<pointField> position(const scalarList&) const;

        //- The length of the edge
        //  \note NotImplemented
        virtual scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
