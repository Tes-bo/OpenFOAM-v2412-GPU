/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::meshWriter

Description
    Write OpenFOAM meshes and/or results to another CFD format.

\par Files

    "constant/boundaryRegion" is an IOMap<dictionary> that contains
    the boundary type and names. eg,
    \verbatim
        (
            0
            {
                BoundaryType    wall;
                Label           Default_Boundary_Region;
            }

            1
            {
                BoundaryType    inlet;
                Label           inlet_1;
            }

            ...

            4
            {
                BoundaryType    pressure;
                Label           outlet;
            }
        )
    \endverbatim


SourceFiles
    meshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef meshWriter_H
#define meshWriter_H

#include "polyMesh.H"
#include "boundaryRegion.H"
#include "cellTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshWriter Declaration
\*---------------------------------------------------------------------------*/

class meshWriter
{
    // Private Member Functions

        //- No copy construct
        meshWriter(const meshWriter&) = delete;

        //- No copy assignment
        void operator=(const meshWriter&) = delete;


protected:

    // Protected data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Scaling factor for points (eg, [m] -> [mm])
        scalar scaleFactor_;

        //- boundaryRegion persistent data saved as a dictionary
        boundaryRegion boundaryRegion_;

        //- cellTable persistent data saved as a dictionary
        cellTable cellTable_;

        //- cellTable IDs for each cell
        labelList cellTableId_;

public:

    // Static data members

        //- Specify a default mesh name
        static string defaultMeshName;


    // Constructors

        //- Create a writer object with given output scaling.
        //  Treats a zero or negative scale factor as unity scaling.
        meshWriter
        (
            const polyMesh&,
            const scalar scaling = 1.0
        );


    //- Destructor
    virtual ~meshWriter();


    // Member Functions

    // Write

        //- Write volume mesh. Subclass must supply this method
        virtual bool write
        (
            const fileName& timeName = fileName::null
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
