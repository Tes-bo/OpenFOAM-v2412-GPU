/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "foamGltfScene.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::scene::scene()
:
    objects_(),
    meshes_(),
    bufferViews_(),
    accessors_(),
    animations_(),
    bytes_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::glTF::mesh& Foam::glTF::scene::getMesh(label meshi)
{
    const label lastMeshi = (meshes_.size() - 1);

    if (meshi < 0)
    {
        meshi = (lastMeshi < 0 ? static_cast<label>(0) : lastMeshi);
    }

    if (meshi > lastMeshi)
    {
        FatalErrorInFunction
            << "Mesh " << meshi << " out of range: " << lastMeshi
            << abort(FatalError);
    }

    return meshes_[meshi];
}


Foam::label Foam::glTF::scene::addColourToMesh
(
    const vectorField& fld,
    const word& name,
    const label meshi,
    const scalarField& alpha
)
{
    auto& gmesh = getMesh(meshi);

    auto& bv = bufferViews_.create(name);
    bv.byteOffset() = bytes_;
    bv.byteLength() = fld.size()*3*sizeof(float); // 3 components
    bv.target() = key(targetTypes::ARRAY_BUFFER);
    bytes_ += bv.byteLength();

    auto& acc = accessors_.create(name);
    acc.bufferViewId() = bv.id();
    acc.set(fld, false); // no min-max

    auto& obj = objects_.create(name);

    if (alpha.empty())
    {
        obj.addData(fld);
    }
    else
    {
        bv.byteLength() += fld.size()*sizeof(float);
        bytes_ += fld.size()*sizeof(float);

        acc.type() = "VEC4";

        // Support uniform alpha vs full alpha field
        tmp<scalarField> talpha(alpha);

        if (alpha.size() == 1 && alpha.size() < fld.size())
        {
            talpha = tmp<scalarField>::New(fld.size(), alpha[0]);
        }

        obj.addData(fld, talpha());
    }

    gmesh.addColour(acc.id());

    return acc.id();
}


Foam::label Foam::glTF::scene::createAnimation(const word& name)
{
    animations_.create(name);
    return animations_.size() - 1;
}


void Foam::glTF::scene::addToAnimation
(
    const label animationi,
    const label inputId,
    const label outputId,
    const label meshId,
    const string& interpolation
)
{
    if (animationi > animations_.size() - 1)
    {
        FatalErrorInFunction
            << "Animation " << animationi << " out of range "
            << (animations_.size() - 1)
            << abort(FatalError);
    }

    const label nodeId = meshId + 1; // offset by 1 for parent node

    // Note
    // using 1 mesh per node +1 parent node => meshes_.size() nodes in total
    if (nodeId > meshes_.size())
    {
        FatalErrorInFunction
            << "Node " << nodeId << " out of range " << meshes_.size()
            << abort(FatalError);
    }

    animations_[animationi].addTranslation
    (
        inputId,
        outputId,
        nodeId,
        interpolation
    );
}


void Foam::glTF::scene::write(const fileName& outputFile)
{
    fileName jsonFile(outputFile);
    jsonFile.replace_ext("gltf");

    // Note: called on master only

    if (!isDir(jsonFile.path()))
    {
        mkDir(jsonFile.path());
    }

    OFstream os(jsonFile);
    write(os);
}


void Foam::glTF::scene::write(Ostream& os)
{
    fileName binFile(os.name());
    binFile.replace_ext("bin");

    // Write binary file
    // Note: using stdStream
    OFstream bin(binFile, IOstreamOption::BINARY);
    auto& osbin = bin.stdStream();

    label totalBytes = 0;
    for (const auto& object : objects_.data())
    {
        for (const auto& data : object.data())
        {
            osbin.write
            (
                reinterpret_cast<const char*>(&data),
                sizeof(float)
            );

            totalBytes += sizeof(float);
        }
    }

    // Write json file
    os << "{" << nl << incrIndent;

    os  << indent << "\"asset\" : {" << nl << incrIndent
        << indent << "\"generator\" : \"OpenFOAM - www.openfoam.com\"," << nl
        << indent << "\"version\" : \"2.0\"" << nl << decrIndent
        << indent << "}," << nl;

    os  << indent << "\"extras\" : {" << nl << incrIndent
        /* << content */
        << decrIndent
        << indent << "}," << nl;

    os  << indent << "\"scene\": 0," << nl;

    os  << indent << "\"scenes\": [{" << nl << incrIndent
        << indent << "\"nodes\" : [0]" << nl << decrIndent
        << indent << "}]," << nl;

    os  << indent << "\"buffers\" : [{" << nl << incrIndent
        << indent << "\"uri\" : " << string(fileName::name(binFile))
        << "," << nl
        << indent << "\"byteLength\" : " << totalBytes << nl << decrIndent
        << indent << "}]," << nl;

    os  << indent << "\"nodes\" : [" << nl << incrIndent
        << indent << "{" << nl << incrIndent
        << indent << "\"children\" : [" << nl << incrIndent;

    // List of child node indices
    os  << indent;
    forAll(meshes_, meshi)
    {
        const label nodeId = meshi + 1;

        os  <<  nodeId;

        if (meshi != meshes_.size() - 1) os  << ", ";

        if ((meshi+1) % 10 == 0) os << nl << indent;
    }

    os  << decrIndent << nl << indent << "]," << nl
        << indent << "\"name\" : \"parent\"" << nl << decrIndent
        << indent << "}," << nl;

    // List of child meshes
    forAll(meshes_, meshi)
    {
        os  << indent << "{" << nl << incrIndent
            << indent << "\"mesh\" : " << meshi << nl << decrIndent
            << indent << "}";

        if (meshi != meshes_.size() - 1) os  << ",";

        os  << nl;
    }

    os  << decrIndent << indent << "]";

    meshes_.write(os, "meshes");

    bufferViews_.write(os, "bufferViews");

    accessors_.write(os, "accessors");

    animations_.write(os, "animations");

    os  << nl;

    os  << decrIndent << "}" << endl;
}


// ************************************************************************* //
