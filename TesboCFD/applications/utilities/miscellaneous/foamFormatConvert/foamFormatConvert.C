/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Application
    foamFormatConvert

Group
    grpMiscUtilities

Description
    Converts all IOobjects associated with a case into the format specified
    in the controlDict.

    Mainly used to convert binary mesh/field files to ASCII.

    Problem: any zero-size List written binary gets written as '0'. When
    reading the file as a dictionary this is interpreted as a label. This
    is (usually) not a problem when doing patch fields since these get the
    'uniform', 'nonuniform' prefix. However zone contents are labelLists
    not labelFields and these go wrong. For now hacked a solution where
    we detect the keywords in zones and redo the dictionary entries
    to be labelLists.

Usage

    - foamFormatConvert [OPTION]

    \param -noConstant \n
    Exclude the constant/ directory from the times list

    \param -enableFunctionEntries \n
    By default all dictionary preprocessing of fields is disabled

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "cellIOList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "cloud.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"
#include "passiveParticleCloud.H"
#include "fieldDictionary.H"

#include "writeMeshObject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// Hack to do zones which have Lists in them. See above.
bool writeZones
(
    const word& name,
    const fileName& meshDir,
    Time& runTime,
    IOstreamOption::compressionType compression
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    bool writeOk = false;

    if (io.typeHeaderOk<cellZoneMesh>(false))
    {
        Info<< "        Reading " << io.headerClassName()
            << " : " << name << endl;

        // Switch off type checking (for reading e.g. faceZones as
        // generic list of dictionaries).
        const word oldTypeName = IOPtrList<entry>::typeName;
        const_cast<word&>(IOPtrList<entry>::typeName) = word::null;

        IOPtrList<entry> meshObject(io);

        for (entry& e : meshObject)
        {
            if (e.isDict())
            {
                dictionary& d = e.dict();

                if (d.found("faceLabels"))
                {
                    d.set("faceLabels", labelList(d.lookup("faceLabels")));
                }

                if (d.found("flipMap"))
                {
                    d.set("flipMap", boolList(d.lookup("flipMap")));
                }

                if (d.found("cellLabels"))
                {
                    d.set("cellLabels", labelList(d.lookup("cellLabels")));
                }

                if (d.found("pointLabels"))
                {
                    d.set("pointLabels", labelList(d.lookup("pointLabels")));
                }
            }
        }

        const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
        // Fake type back to what was in field
        const_cast<word&>(meshObject.type()) = io.headerClassName();

        Info<< "        Writing " << name << endl;

        // Force writing as ASCII
        writeOk = meshObject.regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII, compression),
            true
        );
    }

    return writeOk;
}


// Reduction for non-empty strings.
template<class StringType>
struct uniqueEqOp
{
    void operator()(List<StringType>& x, const List<StringType>& y) const
    {
        List<StringType> newX(x.size()+y.size());
        label n = 0;
        forAll(x, i)
        {
            if (!x[i].empty())
            {
                newX[n++] = x[i];
            }
        }
        forAll(y, i)
        {
            if (!y[i].empty() && !x.found(y[i]))
            {
                newX[n++] = y[i];
            }
        }
        newX.setSize(n);
        x.transfer(newX);
    }
};


template<class Type>
bool writeOptionalMeshObject
(
    const word& name,
    const fileName& meshDir,
    Time& runTime,
    bool writeOnProc
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    // Check if available and the correct type.
    // Done as typeHeaderOk<regIOobject> + isHeaderClass to ensure
    // local-only reading and circumvent is_globalIOobject<Type> check
    // in typeHeaderOk<Type>

    bool readOnProc =
    (
        io.typeHeaderOk<regIOobject>(false)
     && io.isHeaderClass<Type>()
    );

    bool writeOk = false;

    if (returnReduceOr(readOnProc))
    {
        Info<< "        Reading " << Type::typeName << " : " << name << endl;
        Type meshObject(io, readOnProc && writeOnProc);

        Info<< "        Writing " << name << endl;
        writeOk = meshObject.regIOobject::write(readOnProc && writeOnProc);
    }

    return writeOk;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Converts all IOobjects associated with a case into the format"
        " specified in the controlDict"
    );
    timeSelector::addOptions();
    argList::addBoolOption
    (
        "noConstant",
        "Exclude the 'constant/' dir in the times list"
    );
    argList::addBoolOption
    (
        "enableFunctionEntries",
        "Enable expansion of dictionary directives - #include, #codeStream etc"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"

    // enable noConstant by switching
    if (!args.found("noConstant"))
    {
        args.setOption("constant");
    }
    else
    {
        args.unsetOption("constant");
        Info<< "Excluding the constant directory." << nl << endl;
    }


    #include "createTime.H"
    // Optional mesh (used to read Clouds)
    autoPtr<polyMesh> meshPtr;

    const bool enableEntries = args.found("enableFunctionEntries");
    if (enableEntries)
    {
        Info<< "Allowing dictionary preprocessing ('#include', '#codeStream')."
            << endl;
    }

    const int oldFlag = entry::disableFunctionEntries;
    if (!enableEntries)
    {
        // By default disable dictionary expansion for fields
        entry::disableFunctionEntries = 1;
    }

    // Make sure we do not use the master-only reading since we read
    // fields (different per processor) as dictionaries.
    IOobject::fileModificationChecking = IOobject::timeStamp;


    // Specified region or default region
    #include "getRegionOption.H"
    if (!polyMesh::regionName(regionName).empty())
    {
        Info<< "Using region " << regionName << nl << endl;
    }

    const fileName meshDir
    (
        polyMesh::meshDir(regionName)
    );


    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Convert all the standard mesh files
        writeMeshObject<cellCompactIOList, cellIOList>
        (
            "cells",
            meshDir,
            runTime,
            false   // do not check typeName since varies between binary/ascii
        );
        writeMeshObject<labelIOList>("owner", meshDir, runTime);
        writeMeshObject<labelIOList>("neighbour", meshDir, runTime);
        writeMeshObject<faceCompactIOList, faceIOList>
        (
            "faces",
            meshDir,
            runTime,
            false   // do not check typeName since varies between binary/ascii
        );
        writeMeshObject<pointIOField>("points", meshDir, runTime);
        // Write boundary in ascii. This is only needed for fileHandler to
        // kick in. Should not give problems since always writing ascii.
        writeZones("boundary", meshDir, runTime, IOstreamOption::UNCOMPRESSED);
        writeMeshObject<labelIOList>("pointProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>("faceProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>("cellProcAddressing", meshDir, runTime);
        writeMeshObject<labelIOList>
        (
            "boundaryProcAddressing",
            meshDir,
            runTime
        );

        // foamyHexMesh vertices
        writeMeshObject<pointIOField>
        (
            "internalDelaunayVertices",
            polyMesh::regionName(regionName),
            runTime
        );

        if (runTime.writeFormat() == IOstreamOption::ASCII)
        {
            // Only do zones when converting from binary to ascii
            // The other way gives problems since working on dictionary level.
            const IOstreamOption::compressionType compress =
                runTime.writeCompression();
            writeZones("cellZones", meshDir, runTime, compress);
            writeZones("faceZones", meshDir, runTime, compress);
            writeZones("pointZones", meshDir, runTime, compress);
        }

        // Get list of objects from the database
        IOobjectList objects
        (
            runTime,
            runTime.timeName(),
            polyMesh::regionName(regionName)
        );

        forAllConstIters(objects, iter)
        {
            const IOobject& io = *(iter.val());

            if
            (
                io.isHeaderClass<volScalarField>()
             || io.isHeaderClass<volVectorField>()
             || io.isHeaderClass<volSphericalTensorField>()
             || io.isHeaderClass<volSymmTensorField>()
             || io.isHeaderClass<volTensorField>()

             || io.isHeaderClass<surfaceScalarField>()
             || io.isHeaderClass<surfaceVectorField>()
             || io.isHeaderClass<surfaceSphericalTensorField>()
             || io.isHeaderClass<surfaceSymmTensorField>()
             || io.isHeaderClass<surfaceTensorField>()

             || io.isHeaderClass<pointScalarField>()
             || io.isHeaderClass<pointVectorField>()
             || io.isHeaderClass<pointSphericalTensorField>()
             || io.isHeaderClass<pointSymmTensorField>()
             || io.isHeaderClass<pointTensorField>()

             || io.isHeaderClass<volScalarField::Internal>()
             || io.isHeaderClass<volVectorField::Internal>()
             || io.isHeaderClass<volSphericalTensorField::Internal>()
             || io.isHeaderClass<volSymmTensorField::Internal>()
             || io.isHeaderClass<volTensorField::Internal>()
            )
            {
                Info<< "        Reading " << io.headerClassName()
                    << " : " << io.name() << endl;

                fieldDictionary fDict(io, io.headerClassName());

                Info<< "        Writing " << io.name() << endl;
                fDict.regIOobject::write();
            }
        }


        // Check for lagrangian
        fileNameList lagrangianDirs
        (
            1,
            fileHandler().filePath
            (
                runTime.timePath()
              / polyMesh::regionName(regionName)
              / cloud::prefix
            )
        );

        Pstream::combineReduce(lagrangianDirs, uniqueEqOp<fileName>());

        if (!lagrangianDirs.empty())
        {
            if (meshPtr)
            {
                meshPtr->readUpdate();
            }
            else
            {
                Info<< "        Create polyMesh for time = "
                    << runTime.timeName() << endl;

                meshPtr.reset
                (
                    new polyMesh
                    (
                        IOobject
                        (
                            polyMesh::defaultRegion,
                            runTime.timeName(),
                            runTime,
                            IOobject::MUST_READ
                        )
                    )
                );
            }

            const auto& mesh = meshPtr();

            fileNameList cloudDirs
            (
                fileHandler().readDir
                (
                    lagrangianDirs[0],
                    fileName::DIRECTORY
                )
            );

            Pstream::combineReduce(cloudDirs, uniqueEqOp<fileName>());

            for (const auto& cloudDir : cloudDirs)
            {
                fileName dir(cloud::prefix/cloudDir);

                // Read with origProcId,origId fields
                passiveParticleCloud parcels(mesh, cloudDir, true);

                const bool writeOnProc = parcels.size();

                parcels.writeObject
                (
                    IOstreamOption
                    (
                        runTime.writeFormat(),
                        runTime.writeCompression()
                    ),
                    writeOnProc
                );


                // Do local scan for valid cloud objects
                wordList cloudFields
                (
                    IOobjectList(runTime, runTime.timeName(), dir).sortedNames()
                );

                // Combine with all other cloud objects
                Pstream::combineReduce(cloudFields, uniqueEqOp<word>());

                for (const word& name : cloudFields)
                {
                    // These ones already done by cloud itself
                    if
                    (
                        name == "positions"
                     || name == "coordinates"
                     || name == "origProcId"
                     || name == "origId"
                    )
                    {
                        continue;
                    }

                    #undef  doLocalCode
                    #define doLocalCode(Type)                                 \
                    if                                                        \
                    (                                                         \
                        writeOptionalMeshObject<Type>                         \
                        (                                                     \
                            name, dir, runTime, writeOnProc                   \
                        )                                                     \
                    )                                                         \
                    {                                                         \
                        continue;                                             \
                    }

                    doLocalCode(labelIOField);
                    doLocalCode(scalarIOField);
                    doLocalCode(vectorIOField);
                    doLocalCode(sphericalTensorIOField);
                    doLocalCode(symmTensorIOField);
                    doLocalCode(tensorIOField);

                    doLocalCode(labelFieldIOField);
                    doLocalCode(vectorFieldIOField);

                    #undef doLocalCode

                    Info<< "        Failed converting " << name << endl;
                }
            }
        }

        Info<< endl;
    }

    entry::disableFunctionEntries = oldFlag;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
