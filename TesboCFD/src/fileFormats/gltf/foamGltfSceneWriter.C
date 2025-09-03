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

\*---------------------------------------------------------------------------*/

#include "foamGltfSceneWriter.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::sceneWriter::sceneWriter(const fileName& outputFile)
:
    ofile_(nullptr),
    scene_(nullptr)
{
    open(outputFile);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::glTF::sceneWriter::~sceneWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::glTF::sceneWriter::good() const noexcept
{
    return (ofile_ && scene_);
}


const Foam::fileName& Foam::glTF::sceneWriter::path() const
{
    return (ofile_ ? ofile_->name() : fileName::null);
}


const Foam::glTF::scene& Foam::glTF::sceneWriter::getScene() const
{
    return *scene_;
}


Foam::glTF::scene& Foam::glTF::sceneWriter::getScene()
{
    return *scene_;
}


void Foam::glTF::sceneWriter::open(const fileName& outputFile)
{
    close();

    fileName jsonFile(outputFile);
    jsonFile.replace_ext("gltf");

    // Note: called on master only
    if (!isDir(jsonFile.path()))
    {
        mkDir(jsonFile.path());
    }

    ofile_.reset(new OFstream(jsonFile));
    scene_.reset(new glTF::scene());
}


void Foam::glTF::sceneWriter::close()
{
    if (ofile_ && scene_)
    {
        scene_->write(*ofile_);
    }
    ofile_.reset(nullptr);
    scene_.reset(nullptr);
}


// ************************************************************************* //
