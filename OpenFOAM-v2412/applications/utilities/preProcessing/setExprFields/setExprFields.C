/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    setExprFields

Group
    grpPreProcessingUtilities

Description
    Set values on a selected set of cells/patch-faces via a dictionary.

Note
    Based on funkySetFields from
    Bernhard Gschaider <bgschaid@hfd-research.com>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "exprOps.H"
#include "volumeExprDriver.H"
#include "timeSelector.H"
#include "readFields.H"


using namespace Foam;

using FieldAssociation = expressions::FieldAssociation;

word fieldGeoType(const FieldAssociation geoType)
{
    switch (geoType)
    {
        case FieldAssociation::POINT_DATA : return "points"; break;
        case FieldAssociation::FACE_DATA : return "faces"; break;
        case FieldAssociation::VOLUME_DATA : return "cells"; break;
        default: break;
    }
    return "unknown";
}


//- Simple control structure to with collected switches to simplify passing
struct setExprFieldsControl
{
    bool dryRun;
    bool debugParsing;
    bool cacheVariables;
    bool hasDimensions;
    bool createNew;
    bool keepPatches;
    bool correctPatches;
    bool correctBCs;
    IOstreamOption streamOpt;
};


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    if (correctBCs)
    {
        Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, pointPatchField, pointMesh>& field
)
{
    if (correctBCs)
    {
        Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{}


template<class GeoField>
bool setField
(
    const word& fieldName,
    const GeoField& evaluated,
    const boolField& fieldMask,
    const dimensionSet& dims,
    const wordList& valuePatches,

    const setExprFieldsControl& ctrl
)
{
    Info<< "setField(" << fieldName << "): "
        << pTraits<GeoField>::typeName << endl;

    const auto& mesh = evaluated.mesh();

    tmp<GeoField> toutput;

    if (ctrl.createNew)
    {
        // Create with zero
        toutput = GeoField::New
        (
            fieldName,
            mesh,
            Foam::zero{},  // value
            dims
        );
    }
    else
    {
        // Read
        toutput = tmp<GeoField>::New
        (
            IOobject
            (
                fieldName,
                mesh.thisDb().time().timeName(),
                mesh.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh
        );
    }

    auto& output = toutput.ref();

    label numValuesChanged = 0;

    // Internal field
    if (fieldMask.empty())
    {
        // No field-mask - set entire internal field
        numValuesChanged = output.size();

        output.primitiveFieldRef() = evaluated;
    }
    else
    {
        auto& internal = output.primitiveFieldRef();

        forAll(internal, idx)
        {
            if (fieldMask[idx])
            {
                internal[idx] = evaluated[idx];
                ++numValuesChanged;
            }
        }
    }

    // Boundary fields
    forAll(evaluated.boundaryField(), patchi)
    {
        auto& pf = output.boundaryFieldRef()[patchi];

        if (pf.patch().coupled())
        {
            pf == evaluated.boundaryField()[patchi];
        }
    }

    doCorrectBoundaryConditions(ctrl.correctBCs, output);

    const label numTotal = returnReduce(output.size(), sumOp<label>());
    reduce(numValuesChanged, sumOp<label>());

    if (numValuesChanged == numTotal)
    {
        Info<< "Set all ";
    }
    else
    {
        Info<< "Set " << numValuesChanged << " of ";
    }
    Info<< numTotal << " values" << endl;

    if (ctrl.hasDimensions)
    {
        Info<< "Setting dimensions to " << dims << endl;
        output.dimensions().reset(dims);
    }

    if (ctrl.dryRun)
    {
        Info<< "(dry-run): Writing to " << output.name() << nl;
    }
    else
    {
        Info<< "Writing to " << output.name() << nl;
        output.writeObject(ctrl.streamOpt, true);
    }

    return true;
}


void evaluate
(
    const fvMesh& mesh,
    const word& fieldName,
    const expressions::exprString& valueExpr_,
    const expressions::exprString& maskExpr_,
    const dictionary& dict,
    const dimensionSet& dims,
    const wordList& valuePatches,

    const setExprFieldsControl& ctrl
)
{
    word oldFieldType;

    if (ctrl.createNew)
    {
        Info<< "Set new field: " << fieldName;
    }
    else
    {
        IOobject io
        (
            fieldName,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!io.typeHeaderOk<regIOobject>(false))
        {
            FatalErrorInFunction
                << "Field '" << fieldName
                << "' seems to be missing. Use 'create'" << nl
                << exit(FatalError);
        }

        oldFieldType = io.headerClassName();

        Info<< "Modify field: " << fieldName
            << " (type: " << oldFieldType << ')';
    }


    Info<< " time=" << mesh.thisDb().time().timeName() << nl
        << "Expression:" << nl
        << ">>>>" << nl
        << valueExpr_.c_str() << nl
        << "<<<<" << nl;

    bool evalFieldMask =
        (maskExpr_.size() && maskExpr_ != "true" && maskExpr_ != "1");

    if (evalFieldMask)
    {
        Info<< "field-mask:" << nl
            << ">>>>" << nl
            << maskExpr_.c_str() << nl
            << "<<<<" << nl;
    }

    if (ctrl.keepPatches)
    {
        Info<< "Keeping patches unaltered" << endl;
    }
    else if (!valuePatches.empty())
    {
        Info<< "Setting patches " << flatOutput(valuePatches)
            << " to fixed value" << endl;
    }

    Info<< endl;

    expressions::volumeExprDriver driver(mesh);

    driver.setCaching(ctrl.cacheVariables);

    driver.readDict(dict);

    if (ctrl.debugParsing)
    {
        Info<< "Parsing expression: " << valueExpr_ << "\nand field-mask "
            << maskExpr_ << nl << endl;
        driver.setDebugging(true, true);
    }


    driver.clearVariables();


    // Handle "field-mask" evaluation

    boolField fieldMask;
    FieldAssociation maskFieldAssoc(FieldAssociation::NO_DATA);

    if (evalFieldMask)
    {
        if (ctrl.debugParsing)
        {
            Info<< "Parsing field-mask:" << maskExpr_ << endl;
        }

        driver.parse(maskExpr_);
        if (ctrl.debugParsing)
        {
            Info<< "Parsed field-mask" << endl;
        }

        if (driver.isLogical())
        {
            auto& result = driver.result();
            if (result.is_bool())
            {
                fieldMask = result.getResult<bool>();
                maskFieldAssoc = driver.fieldAssociation();
            }
        }

        // Slightly pedantic...
        driver.clearField();
        driver.clearResult();

        evalFieldMask = (maskFieldAssoc != FieldAssociation::NO_DATA);

        if (!evalFieldMask)
        {
            FatalErrorInFunction
                << " mask: " << maskExpr_
                << " does not evaluate to a logical expression: "
                << driver.resultType() << nl
                #ifdef FULLDEBUG
                << "contents: " << fieldMask
                #endif
                << exit(FatalError);
        }

        if (ctrl.debugParsing)
        {
            Info<< "Field-mask evaluates to "
                << fieldMask << nl;
        }
    }

    if (ctrl.debugParsing)
    {
        Info<< "Parsing expression:" << valueExpr_ << endl;
    }

    driver.parse(valueExpr_);

    if (ctrl.debugParsing)
    {
        Info<< "Parsed expression" << endl;
    }

    if (evalFieldMask && maskFieldAssoc != driver.fieldAssociation())
    {
        FatalErrorInFunction
            << "Mismatch between field-mask geometric type ("
            << fieldGeoType(maskFieldAssoc) << ") and" << nl
            << "expression geometric type ("
            << fieldGeoType(driver.fieldAssociation()) << ')' << nl
            << nl
            << "expression: " << valueExpr_ << nl
            << "field-mask: " << maskExpr_ << nl
            << nl
            << exit(FatalError);
    }

    if (!oldFieldType.empty() && driver.resultType() != oldFieldType)
    {
        FatalErrorInFunction
            << "Inconsistent types: " << fieldName << " is  "
            << oldFieldType
            << " but the expression evaluates to "
            << driver.resultType()
            << exit(FatalError);
    }

    Info<< "Dispatch ... " << driver.resultType() << nl;


    bool applied = false;
    switch (driver.fieldAssociation())
    {
        #undef  doLocalCode
        #define doLocalCode(GeoField)                                        \
        {                                                                    \
            const auto* ptr = driver.isResultType<GeoField>();               \
            if (ptr)                                                         \
            {                                                                \
                applied = setField                                           \
                (                                                            \
                    fieldName,                                               \
                    *ptr,                                                    \
                    fieldMask,                                               \
                    dims,                                                    \
                    valuePatches,                                            \
                    ctrl                                                     \
                );                                                           \
                break;                                                       \
            }                                                                \
        }

        case FieldAssociation::VOLUME_DATA:
        {
            doLocalCode(volScalarField);
            doLocalCode(volVectorField);
            doLocalCode(volTensorField);
            doLocalCode(volSymmTensorField);
            doLocalCode(volSphericalTensorField);
            break;
        }
        case FieldAssociation::FACE_DATA:
        {
            doLocalCode(surfaceScalarField);
            doLocalCode(surfaceVectorField);
            doLocalCode(surfaceTensorField);
            doLocalCode(surfaceSymmTensorField);
            doLocalCode(surfaceSphericalTensorField);
            break;
        }
        case FieldAssociation::POINT_DATA:
        {
            doLocalCode(pointScalarField);
            doLocalCode(pointVectorField);
            doLocalCode(pointTensorField);
            doLocalCode(pointSymmTensorField);
            doLocalCode(pointSphericalTensorField);
            break;
        }

        default: break;
        #undef doLocalCode
    }

    if (!applied)
    {
        FatalErrorInFunction
            << "Expression evaluates to an unsupported type: "
            << driver.resultType() << nl << nl
            << "Expression " << valueExpr_ << nl << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Normally without functionObjects, with -withFunctionObjects to enable
    argList::noFunctionObjects(true);

    // No -constant, no special treatment for 0/
    timeSelector::addOptions(false);

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of the controlDict setting"
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for setExprFieldsDict"
    );
    argList::addDryRunOption
    (
        "Evaluate but do not write"
    );
    argList::addVerboseOption
    (
        "Additional verbosity",
        true // Advanced option
    );
    argList::addOption
    (
        "load-fields",
        "wordList",
        "Specify field or fields to preload. Eg, 'T' or '(p T U)'",
        true // Advanced option
    );
    argList::addOption
    (
        "field",
        "name",
        "The field to create/overwrite"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "expression",
        "expr",
        "The expression to evaluate"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "field-mask",
        "logic",
        "The field mask (logical condition) when to apply the expression"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOptionCompat("field-mask", {"condition", 2106});
    argList::addOption
    (
        "dimensions",
        "dims",
        "The dimensions for created fields"
        " (command-line operation)",
        true // Advanced option
    );

    argList::addBoolOption
    (
        "debug-parser",
        "Additional debugging information",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "no-variable-cache",
        "Disable caching of expression variables",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "create",
        "Create a new field"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "keepPatches",
        "Leave patches unaltered"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "value-patches",
        "(patches)",
        "A list of patches that receive a fixed value"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "dummy-phi",
        "Provide a zero phi field"
        " (command-line operation)",
        true // Advanced option
    );

    // Future?
    #if 0
    argList::addBoolOption
    (
        "noCorrectPatches",
        ""
    );
    argList::addBoolOption
    (
        "correctResultBoundaryFields",
        "",
        true
    );
    #endif

    #include "addRegionOption.H"
    #include "setRootCase.H"

    #include "createTime.H"

    const word dictName("setExprFieldsDict");

    instantList times = timeSelector::select0(runTime, args);

    if (times.empty())
    {
        FatalErrorInFunction
            << "No times selected." << nl
            << exit(FatalError);
    }

    // Disable dimension checking during operations
    dimensionSet::checking(false);

    #include "createNamedMesh.H"

    autoPtr<surfaceScalarField> dummyPhi;

    autoPtr<IOdictionary> exprDictPtr;

    // Sort out conflicts

    const bool useCommandArgs = args.found("field");

    if (useCommandArgs)
    {
        bool fatalCombination = false;

        if (args.found("dict"))
        {
            fatalCombination = true;
            FatalErrorInFunction
                << "Cannot specify both dictionary and command-line arguments"
                << nl << endl;
        }

        if (args.found("create") && args.found("keepPatches"))
        {
            fatalCombination = true;
            FatalErrorInFunction
                << "Cannot specify both 'create' and 'keepPatches'" << nl
                << endl;
        }

        if (!args.found("expression"))
        {
            fatalCombination = true;
            FatalErrorInFunction
                << "Missing mandatory 'expression' option'" << nl
                << endl;
        }
        if (fatalCombination)
        {
            FatalError
                << exit(FatalError);
        }
    }
    else
    {
        // Carp about inapplicable options (non-fatal)

        wordHashSet badOptions
        ({
            "create", "keepPatches", "value-patches",
            "field-mask", "expression", "dimensions"
        });
        badOptions.retain(args.options());

        if (!badOptions.empty())
        {
            // Non-fatal (warning)
            FatalErrorInFunction
                << "Using a dictionary. Cannot specify these options:" << nl
                << flatOutput(badOptions.sortedToc()) << nl
                << endl;
        }

        #include "setSystemMeshDictionaryIO.H"
        exprDictPtr.reset(new IOdictionary(dictIO));
    }

    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "\nTime = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // preload fields specified on command-line
        if (timei == 0)
        {
            wordList preloadFields;
            args.readListIfPresent("load-fields", preloadFields);
            readFieldsHandler(mesh).execute(preloadFields);
        }

        if (args.found("dummy-phi") && !dummyPhi)
        {
            Info<< "Adding a dummy phi" << endl;
            dummyPhi.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phi",
                        mesh.thisDb().time().constant(),
                        mesh.thisDb(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(Zero)
                )
            );
        }

        if (args.allowFunctionObjects())
        {
            runTime.functionObjects().start();
        }

        if (useCommandArgs)
        {
            const word fieldName(args.get<word>("field"));

            Info<< "Using command-line options for "
                << fieldName << nl << endl;

            setExprFieldsControl ctrl;

            ctrl.dryRun = args.dryRun();
            ctrl.debugParsing = args.found("debug-parser");
            ctrl.cacheVariables = !args.found("no-variable-caching");

            ctrl.createNew = args.found("create");
            ctrl.keepPatches = args.found("keepPatches");
            ctrl.correctPatches = !args.found("noCorrectPatches");
            ctrl.correctBCs = args.found("correctResultBoundaryFields");
            ctrl.hasDimensions = args.found("dimensions");
            ctrl.streamOpt.format(runTime.writeFormat());
            if (args.found("ascii"))
            {
                ctrl.streamOpt.format(IOstreamOption::ASCII);
            }

            expressions::exprString valueExpr_(args["expression"]);

            expressions::exprString maskExpr_;
            args.readIfPresent("field-mask", maskExpr_);

            dimensionSet dims;
            if (ctrl.hasDimensions)
            {
                ITstream is(args.lookup("dimensions"));
                is >> dims;
            }

            evaluate
            (
                mesh,
                fieldName,
                valueExpr_,
                maskExpr_,
                dictionary::null,
                dims,
                args.getList<word>("value-patches", false),

                ctrl
            );
        }
        else if (exprDictPtr)
        {
            const dictionary& exprDict = *exprDictPtr;

            // preload fields specified in dictionary
            {
                wordList preloadFields;
                exprDict.readIfPresent("readFields", preloadFields);
                readFieldsHandler(mesh).execute(preloadFields);
            }

            // Read set construct info from dictionary
            PtrList<entry> actions(exprDict.lookup("expressions"));

            for (const entry& dEntry : actions)
            {
                if (!dEntry.isDict())
                {
                    Info<< "Ignore non-dictionary entry: "
                        << dEntry.keyword() << nl;
                    continue;
                }

                const dictionary& dict = dEntry.dict();

                setExprFieldsControl ctrl;

                ctrl.dryRun = args.dryRun();
                ctrl.debugParsing = args.found("debug-parser");
                ctrl.cacheVariables = !args.found("no-variable-caching");

                ctrl.createNew = dict.getOrDefault("create", false);
                ctrl.keepPatches = dict.getOrDefault("keepPatches", false);

                ctrl.correctPatches = !args.found("noCorrectPatches");
                ctrl.correctBCs = args.found("correctResultBoundaryFields");
                ctrl.streamOpt.format(runTime.writeFormat());
                if (args.found("ascii"))
                {
                    ctrl.streamOpt.format(IOstreamOption::ASCII);
                }

                if (ctrl.createNew && ctrl.keepPatches)
                {
                    FatalIOErrorInFunction(dict)
                        << "Cannot specify both 'create' and 'keepPatches'"
                        << nl << endl
                        << exit(FatalIOError);
                }

                // Local override
                dict.readIfPresent
                (
                    "correctResultBoundaryFields",
                    ctrl.correctBCs
                );


                const word fieldName(dict.get<word>("field"));

                expressions::exprString valueExpr_("expression", dict);

                expressions::exprString maskExpr_;
                {
                    const entry* eptr = dict.findCompat
                    (
                        "fieldMask", {{"condition", 2106}},
                        keyType::LITERAL
                    );

                    if (eptr)
                    {
                        maskExpr_.readEntry(eptr->keyword(), dict);
                        maskExpr_.trim();
                    }
                }

                // Optional: "dimensions"
                dimensionSet dims;
                ctrl.hasDimensions = dims.readIfPresent("dimensions", dict);

                if (args.verbose() && !timei)
                {
                    // Report once
                    Info<< "Processing" << dict << nl;
                }

                evaluate
                (
                    mesh,
                    fieldName,
                    valueExpr_,
                    maskExpr_,
                    dict,
                    dims,
                    dict.getOrDefault<wordList>("valuePatches", wordList()),

                    ctrl
                );
            }
        }
        else
        {
            FatalErrorInFunction
                << "No command-line or dictionary??" << nl << endl
                << exit(FatalError);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
