/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    vtkUnstructuredToFoam

Group
    grpMeshConversionUtilities

Description
    Convert legacy VTK file (ascii) containing an unstructured grid
    to an OpenFOAM mesh without boundary information.

Usage
    \b vtkUnstructuredToFoam \<XXX.vtk\>

    Options:
      - \par -no-fields
        Do not attempt to recreate volFields

Note
    The .vtk format does not contain any boundary information.
    It is purely a description of the internal mesh. This also limits the
    usefulness of reconstructing the volFields.

    Not extensively tested.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "vtkUnstructuredReader.H"

#include "columnFvMesh.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void constructVolFields(fvMesh& mesh, const vtkUnstructuredReader& reader)
{
    const auto fields(reader.cellData().csorted<IOField<Type>>());
    for (const auto& field : fields)
    {
        Info<< "Constructing volField " << field.name() << endl;

        // field is
        //  - cell data followed by
        //  - boundary face data


        auto tfld = GeometricField<Type, fvPatchField, volMesh>::New
        (
            field.name(),
            mesh,
            dimless
        );
        auto& fld = tfld.ref();
        fld.instance() = mesh.time().timeName();
        fld.writeOpt() = IOobject::AUTO_WRITE;

        // Fill cell values
        fld.internalFieldRef().field() =
            UIndirectList<Type>(field, reader.cellMap());

        // Fill boundary values
        const auto& map = reader.faceMap();
        if (map.size())
        {
            for (auto& pfld : fld.boundaryFieldRef())
            {
                const auto& pp = pfld.patch();

                forAll(pfld, i)
                {
                    const label bFacei = pp.patch().offset()+i;
                    pfld[i] = field[map[bFacei]];
                }
            }
        }

        regIOobject::store(std::move(tfld));
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert legacy VTK file (ascii) containing an unstructured grid"
        " to an OpenFOAM mesh without boundary information"
    );

    argList::noParallel();
    argList::addOptionCompat("no-fields", {"noFields", 2106});
    argList::addArgument("vtk-file", "The input legacy ascii vtk file");

    #include "setRootCase.H"
    #include "createTime.H"

    const bool doFields = !args.found("no-fields");

    IFstream mshStream(args.get<fileName>(1));

    vtkUnstructuredReader reader(runTime, mshStream);

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        std::move(reader.points()),
        reader.cells(),
        faceListList(),
        wordList(),
        wordList(),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing mesh ..." << endl;

    mesh.removeFiles();
    mesh.write();


    if (doFields)
    {
        // Re-read mesh as fvMesh so we can have fields
        Info<< "Re-reading mesh ..." << endl;
        #include "createMesh.H"

        constructVolFields<scalar>(mesh, reader);
        constructVolFields<vector>(mesh, reader);
        constructVolFields<sphericalTensor>(mesh, reader);
        constructVolFields<symmTensor>(mesh, reader);
        constructVolFields<tensor>(mesh, reader);

        // No need to write the mesh, only fields
        mesh.thisDb().write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
