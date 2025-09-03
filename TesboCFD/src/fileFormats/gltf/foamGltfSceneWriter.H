/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::glTF::sceneWriter

Description
    Wrapper for glTF scene for file output

SourceFiles
    foamGltfSceneWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gltf_sceneWriter_H
#define Foam_gltf_sceneWriter_H

#include "autoPtr.H"
#include "fileName.H"
#include "foamGltfScene.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OFstream;

namespace glTF
{

/*---------------------------------------------------------------------------*\
                     Class glTF::sceneWriter Declaration
\*---------------------------------------------------------------------------*/

class sceneWriter
{
    // Private Data

        //- The backend output stream (json)
        autoPtr<OFstream> ofile_;

        //- The scene to output
        autoPtr<glTF::scene> scene_;


public:

    // Generated Methods

        //- No copy construct
        sceneWriter(const sceneWriter&) = delete;

        //- No copy assignment
        void operator=(const sceneWriter&) = delete;


    // Constructors

        //- Default construct
        sceneWriter() = default;

        //- Construct and open with given file name
        explicit sceneWriter(const fileName& outputFile);


    //- Destructor - calls close()
    ~sceneWriter();


    // Member Functions

        //- True if output file and scene exist
        bool good() const noexcept;

        //- The json file name. Empty with !good()
        const fileName& path() const;

        //- Const access to the scene. Error if good() is not true!
        const scene& getScene() const;

        //- Non-const access to the scene. Error if good() is not true!
        scene& getScene();


        //- Flush, output and open a new file for output
        void open(const fileName& outputFile);

        //- Write scene and close file
        void close();

        //- Same as good()
        bool valid() const noexcept { return good(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
