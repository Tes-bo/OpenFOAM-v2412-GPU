/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::cellPointConnectivity

Description
    This class provides ordered connectivity for each point of each cell. Lists
    are available of the points and faces surrounding each point of a cell. The
    lists are ordered so that the connected points describe a polygonal cone.
    For a convex cell, any three sequantial cone edges form a positive basis.

SourceFiles
    cellPointConnectivity.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointConnectivity_H
#define cellPointConnectivity_H

#include "polyMesh.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cellPointConnectivity Declaration
\*---------------------------------------------------------------------------*/

class cellPointConnectivity
:
    public MoveableMeshObject<polyMesh>
{
    // Private data

        //- Reference to the polyMesh
        const polyMesh& mesh_;

        //- Lists of point-point connections
        labelListListList cellPointPoints_;

        //- Lists of point-face connections
        labelListListList cellPointFaces_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellPointConnectivity(const cellPointConnectivity&);

        //- Disallow default bitwise assignment
        void operator=(const cellPointConnectivity&);

        //- Generate the connectivity
        void generateCellPointConnectivity(label cellI);


public:

    // Run-time type information
    TypeName("cellPointConnectivity");


    // Constructors

        //- Construct an IOobject and a polymesh
        cellPointConnectivity(const polyMesh&);


    //- Destructor
    ~cellPointConnectivity();


    // Member Functions

        // Access

            //- Access the point-point connections
            const labelListListList& cellPointPoints() const
            {
                return cellPointPoints_;
            }

            //- Access the point-face connections
            const labelListListList& cellPointFaces() const
            {
                return cellPointFaces_;
            }

        // Edit

            //- No action required on move points
            virtual bool movePoints()
            {
                return true;
            }

            //- Dummy write for regIOobject
            virtual bool writeData(Ostream&) const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
