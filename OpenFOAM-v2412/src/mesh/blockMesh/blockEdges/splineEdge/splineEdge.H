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
    Foam::blockEdges::splineEdge

Description
    A blockEdge interface for Catmull-Rom splines.

See also
    BSpline
    CatmullRomSpline

SourceFiles
    splineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_splineEdge_H
#define blockEdges_splineEdge_H

#include "blockEdge.H"
#include "CatmullRomSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                      Class splineEdge Declaration
\*---------------------------------------------------------------------------*/

class splineEdge
:
    public blockEdge,
    public CatmullRomSpline
{
    // Private Member Functions

        //- No copy construct
        splineEdge(const splineEdge&) = delete;

        //- No copy assignment
        void operator=(const splineEdge&) = delete;


public:

    //- Runtime type information
    TypeName("spline");


    // Constructors

        //- Construct from components
        splineEdge
        (
            const pointField& points,   //!< Referenced point field
            const edge& fromTo,         //!< Start/end in point field
            const pointField& internalPoints
        );

        //- Construct from components
        splineEdge
        (
            const pointField& points,   //!< Referenced point field
            const label from,           //!< Start point in point field
            const label to,             //!< End point in point field
            const pointField& internalPoints
        );

        //- Construct from Istream and point field.
        splineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces&  /*unused*/,
            const pointField& points,  //!< Referenced point field
            Istream&
        );


    //- Destructor
    virtual ~splineEdge() = default;


    // Member Functions

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- The length of the spline curve
        //  \note NotImplemented
        virtual scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
