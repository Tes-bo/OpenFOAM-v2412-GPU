/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    patchSummary

Group
    grpMiscUtilities

Description
    Write field and boundary condition info for each patch at each requested
    time instance.

    Default action is to write a single entry for patches/patchGroups with the
    same boundary conditions. Use the -expand option to print every patch
    separately. In case of multiple groups matching it will print only the
    first one.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "pointFields.H"
#include "IOobjectList.H"
#include "patchSummaryTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write field and boundary condition info for each patch"
        " at each requested time instance"
    );

    timeSelector::addOptions();

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "expand",
        "Do not combine patches"
    );
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool expand = args.found("expand");


    #include "createNamedMesh.H"
    const polyBoundaryMesh& bm = mesh.boundaryMesh();


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update the mesh if changed
        if (mesh.readUpdate() == polyMesh::TOPO_PATCH_CHANGE)
        {
            Info<< "Detected changed patches. Recreating patch group table."
                << endl;
        }

        const IOobjectList objects(mesh, runTime.timeName());

        Info<< "Reading fields:" << endl;

        // Read fields
        #undef  createFields
        #define createFields(FieldType, Variable)                   \
        PtrList<FieldType> Variable                                 \
        (                                                           \
            readFields<FieldType>(objects, mesh)                    \
        );

        createFields(volScalarField, vsf);
        createFields(volVectorField, vvf);
        createFields(volSphericalTensorField, vsptf);
        createFields(volSymmTensorField, vsytf);
        createFields(volTensorField, vtf);

        // Point fields
        const pointMesh& pMesh = pointMesh::New(mesh);

        #undef  createFields
        #define createFields(FieldType, Variable)                   \
        PtrList<FieldType> Variable                                 \
        (                                                           \
            readFields<FieldType>(objects, pMesh)                   \
        );

        createFields(pointScalarField, psf);
        createFields(pointVectorField, pvf);
        createFields(pointSphericalTensorField, psptf);
        createFields(pointSymmTensorField, psytf);
        createFields(pointTensorField, ptf);

        #undef createFields

        Info<< endl;


        if (expand)
        {
            // Print each patch separately

            forAll(bm, patchi)
            {
                Info<< bm[patchi].type() << "\t: " << bm[patchi].name() << nl;

                outputFieldList(vsf, patchi);
                outputFieldList(vvf, patchi);
                outputFieldList(vsptf, patchi);
                outputFieldList(vsytf, patchi);
                outputFieldList(vtf, patchi);

                outputFieldList(psf, patchi);
                outputFieldList(pvf, patchi);
                outputFieldList(psptf, patchi);
                outputFieldList(psytf, patchi);
                outputFieldList(ptf, patchi);
                Info<< endl;
            }
        }
        else
        {
            // Collect for each patch the bc type per field. Merge similar
            // patches.

            // Per 'group', the map from fieldname to patchfield type
            DynamicList<HashTable<word>> fieldToTypes(bm.size());
            // Per 'group' the patches
            DynamicList<DynamicList<label>> groupToPatches(bm.size());

            forAll(bm, patchi)
            {
                HashTable<word> fieldToType;
                collectFieldList(vsf, patchi, fieldToType);
                collectFieldList(vvf, patchi, fieldToType);
                collectFieldList(vsptf, patchi, fieldToType);
                collectFieldList(vsytf, patchi, fieldToType);
                collectFieldList(vtf, patchi, fieldToType);

                collectFieldList(psf, patchi, fieldToType);
                collectFieldList(pvf, patchi, fieldToType);
                collectFieldList(psptf, patchi, fieldToType);
                collectFieldList(psytf, patchi, fieldToType);
                collectFieldList(ptf, patchi, fieldToType);

                label groupI = fieldToTypes.find(fieldToType);
                if (groupI == -1)
                {
                    DynamicList<label> group(1);
                    group.append(patchi);
                    groupToPatches.append(group);
                    fieldToTypes.append(fieldToType);
                }
                else
                {
                    groupToPatches[groupI].append(patchi);
                }
            }


            forAll(groupToPatches, groupI)
            {
                const DynamicList<label>& patchIDs = groupToPatches[groupI];

                if (patchIDs.size() > 1)
                {
                    // Check if part of a group
                    wordList groups;
                    labelHashSet nonGroupPatches;
                    bm.matchGroups(patchIDs, groups, nonGroupPatches);

                    for (const label patchi : nonGroupPatches.sortedToc())
                    {
                        Info<< bm[patchi].type()
                            << "\t: " << bm[patchi].name() << nl;
                    }
                    for (const word& groupName : groups)
                    {
                        Info<< "group\t: " << groupName << nl;
                    }

                    const label patchi = patchIDs[0];

                    outputFieldList(vsf, patchi);
                    outputFieldList(vvf, patchi);
                    outputFieldList(vsptf, patchi);
                    outputFieldList(vsytf, patchi);
                    outputFieldList(vtf, patchi);

                    outputFieldList(psf, patchi);
                    outputFieldList(pvf, patchi);
                    outputFieldList(psptf, patchi);
                    outputFieldList(psytf, patchi);
                    outputFieldList(ptf, patchi);
                    Info<< endl;
                }
                else
                {
                    // No group.
                    for (const label patchi : patchIDs)
                    {
                        Info<< bm[patchi].type()
                            << "\t: " << bm[patchi].name() << nl;

                        outputFieldList(vsf, patchi);
                        outputFieldList(vvf, patchi);
                        outputFieldList(vsptf, patchi);
                        outputFieldList(vsytf, patchi);
                        outputFieldList(vtf, patchi);

                        outputFieldList(psf, patchi);
                        outputFieldList(pvf, patchi);
                        outputFieldList(psptf, patchi);
                        outputFieldList(psytf, patchi);
                        outputFieldList(ptf, patchi);
                        Info<< endl;
                    }
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
