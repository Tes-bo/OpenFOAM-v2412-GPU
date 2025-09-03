/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::geometrySurface

Description
    Read and visualize surface geometry files.

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | The type: geometry                    | yes |
        files       | The files to read                     | yes |
    \endtable

    Standard file types (vtk, vtp, obj, stl, stlb) are read with the
    VTK-native readers. Other file types use the OpenFOAM surfMesh
    readers and convert to VTK.

SourceFiles
    geometrySurface.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_geometrySurface_H
#define functionObjects_runTimePostPro_geometrySurface_H

#include "surface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkPolyData;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                       Class geometrySurface Declaration
\*---------------------------------------------------------------------------*/

class geometrySurface
:
    public surface
{
protected:

    // Protected Data

        //- File names
        List<fileName> fileNames_;


    // Protected Member Functions

        //- Add surface (file) to scene
        void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer,
            const fileName& fName
        ) const;

        //- No copy construct
        geometrySurface(const geometrySurface&) = delete;

        //- No copy assignment
        void operator=(const geometrySurface&) = delete;


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("geometry");


    // Constructors

        //- Construct from dictionary
        geometrySurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );

        //- Construct from components
        geometrySurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours,
            const List<fileName>& fileNames
        );


    //- Destructor
    virtual ~geometrySurface() = default;


    // Member Functions

        //- Add geometry surface(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- Clear files used to create the object(s) - no-op
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
