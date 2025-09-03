/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    createZeroDirectory

Group
    grpPreProcessingUtilities

Description
    Creates a zero directory with fields appropriate for the chosen solver and
    turbulence model. Operates on both single and multi-region cases.

Usage
    The set-up is configured using a 'caseProperties' dictionary, located under
    the $FOAM_CASE/system (or system/regionName if multi-region) directory.
    This consists of a lists of initial and boundary conditions, e.g.

    \verbatim
    initialConditions
    {
        U           uniform (0 0 0);
        p           uniform 0;
    }

    boundaryConditions
    {
        topWall
        {
            category        wall;
            patches         (movingWall);
            type            noSlip;
            options
            {
                wallFunction    highReynolds;
                motion          moving;
            };
            values
            {
                U           uniform (1 0 0);
            }
        }

        walls
        {
            category        wall;
            patches         (fixedWalls);
            type            noSlip;
            options
            {
                wallFunction    highReynolds;
                motion          stationary;
            };
        }
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "caseInfo.H"
#include "boundaryTemplates.H"
#include "solverTemplate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const word getClassType(const word& pType)
{
    if (pType == pTraits<scalar>::typeName)
    {
        return GeometricField<scalar, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<vector>::typeName)
    {
        return GeometricField<vector, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<sphericalTensor>::typeName)
    {
        return GeometricField<sphericalTensor, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<symmTensor>::typeName)
    {
        return GeometricField<symmTensor, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<tensor>::typeName)
    {
        return GeometricField<tensor, fvPatchField, volMesh>::typeName;
    }
    else
    {
        // Error
        return word::null;
    }
}


void createFieldFiles
(
    const Time& runTime,
    const word& regionName,
    const word& regionPrefix,
    const wordList& fieldNames,
    const wordList& fieldTypes,
    const PtrList<dimensionSet>& fieldDimensions,
    const caseInfo& cInfo,
    const boundaryTemplates& bcTemplates
)
{
    Info<< "    Generating field files" << nl << endl;

    // Create files
    forAll(fieldNames, i)
    {
        const_cast<word&>(IOdictionary::typeName) =
            getClassType(fieldTypes[i]);

        dictionary field;

        field.add
        (
            "#include",
            "<system>"/regionName/"caseProperties"
        );

        field.add("dimensions", fieldDimensions[i]);

        string iField("${:initialConditions." + fieldNames[i] + '}');
        field.add("internalField", iField.c_str());

        dictionary boundaryField =
            cInfo.generateBoundaryField
            (
                regionPrefix,
                fieldNames[i],
                bcTemplates
            );

        field.add("boundaryField", boundaryField);

        // Expand all of the dictionary redirections and remove unnecessary
        // entries

        dictionary field2;

        {
            OCharStream os;
            os  << field;
            ISpanStream is(os.view());

            entry::disableFunctionEntries = 0;
            is >> field2;
            entry::disableFunctionEntries = 1;

            field2.remove("#include");
            field2.remove("initialConditions");
            field2.remove("boundaryConditions");
        }

        // Construct and write field dictionary. Note use of localIOdictionary
        localIOdictionary fieldOut
        (
            IOobject
            (
                fieldNames[i],
                "0",
                regionName,
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            field2
        );

        fieldOut.regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII),
            true
        );
    }
}


// Main program:
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a 0/ directory with fields appropriate for the chosen"
        " solver and turbulence model."
    );

    argList::addOption
    (
        "templateDir",
        "dir",
        "Read case set-up templates from specified location"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Reading controlDict" << nl << endl;

    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // Template directory: default is from PROJECT/etc directory
    //
    // Can use
    // - foamEtcDir("caseDicts/createZeroDirectory", 0007);
    // - expand "<etc:o>/caseDicts/createZeroDirectory"
    // - expand "${WM_PROJECT_DIR}/etc/caseDicts/createZeroDirectory"
    //
    // Use "${WM_PROJECT_DIR}/" version for nicer error message

    fileName baseDir
    (
        args.getOrDefault<fileName>
        (
            "templateDir",
            "${WM_PROJECT_DIR}/etc/caseDicts/createZeroDirectoryTemplates"
        )
    );

    baseDir.expand();
    baseDir.toAbsolute();

    if (!Foam::isDir(baseDir))
    {
        FatalErrorInFunction
            << "templateDir " << baseDir << nl
            << "Does not point to a folder containing case set-up templates"
            << nl
            << exit(FatalError);
    }

    // Keep variable substitutions - delay until after creation of controlDict
    // to allow #include files to be processed
    entry::disableFunctionEntries = 1;

    // Read the solver
    const word solverName(controlDict.get<word>("application"));

    // Generate solver template
    const solverTemplate solver(baseDir, runTime, solverName);

    // Read the boundary condition templates
    const boundaryTemplates bcTemplates
    (
        baseDir,
        runTime,
        solver.type()
    );

    Info<< endl;

    const label nRegion = solver.nRegion();
    for (label regionI = 0; regionI < nRegion; regionI++)
    {
        const word& regionName = solver.regionName(regionI);

        if (regionName.empty())
        {
            Info<< "Region: " << polyMesh::defaultRegion << " (default)"
                << endl;
        }
        else
        {
            Info<< "Region: " << regionName << endl;
        }

        // Read the case set-up info for the current region
        const caseInfo cInfo(runTime, regionName);

        cInfo.checkPatches(solver.regionType(regionI), bcTemplates);

        createFieldFiles
        (
            runTime,
            regionName,
            solver.regionType(regionI),
            solver.fieldNames(regionI),
            solver.fieldTypes(regionI),
            solver.fieldDimensions(regionI),
            cInfo,
            bcTemplates
        );
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
