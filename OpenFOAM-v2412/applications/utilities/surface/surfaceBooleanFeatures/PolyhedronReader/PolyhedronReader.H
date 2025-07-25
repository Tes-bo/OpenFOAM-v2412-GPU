/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::PolyhedronReader

Description

SourceFiles
    PolyhedronReader.C

\*---------------------------------------------------------------------------*/

#ifndef PolyhedronReader_H
#define PolyhedronReader_H

#include "CGALIndexedPolyhedron.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PolyhedronReader Declaration
\*---------------------------------------------------------------------------*/

class PolyhedronReader
{
    // Private Classes

        template<class HDS>
        class Build_triangle
        :
            public CGAL::Modifier_base<HDS>
        {
            const triSurface& s_;

        public:

            Build_triangle(const triSurface& s)
            :
                s_(s)
            {}

            void operator()(HDS& hds);
        };


public:

    // Constructors

        //- Populate polyhedron from surface
        PolyhedronReader(const triSurface& s, Polyhedron& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
