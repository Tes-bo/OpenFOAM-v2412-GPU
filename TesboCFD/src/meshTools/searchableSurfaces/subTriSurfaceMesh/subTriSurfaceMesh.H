/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::subTriSurfaceMesh

Description
    On-the-fly subset of triSurfaceMesh (based on selected patches)

Note

    Example usage:
    \verbatim
    spheres_solid1
    {
        type subTriSurfaceMesh;

        // Surface to operate on. Avoid duplicate loading
        surface spheres.stl;

        // Regions to operate on (regular expressions allowed)
        patches (solid1);
    }
    \endverbatim

    \heading Dictionary parameters
    \table
        Property | Description                              | Required | Default
        type     | subTriSurfaceMesh                        | selector |
        surface  | Name of the underlying surface           | yes   |
        patches  | List of surface region names or regexs   | yes   |
    \endtable

SourceFiles
    subTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef subTriSurfaceMesh_H
#define subTriSurfaceMesh_H

#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class subTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class subTriSurfaceMesh
:
    public triSurfaceMesh
{
    // Private member functions

        //- Extract patch names of triSurface
        static wordList patchNames(const triSurface& s);

        //- Select regions by name
        static labelList selectedRegions
        (
            const triSurface& s,
            const wordRes& regionNameMatcher
        );

        //- Subset triSurface based on regions
        static triSurface subset(const IOobject&, const dictionary&);


public:

    //- Runtime type information
    TypeName("subTriSurfaceMesh");


    // Constructors

        //- Construct from IO and dictionary.
        subTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~subTriSurfaceMesh() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
