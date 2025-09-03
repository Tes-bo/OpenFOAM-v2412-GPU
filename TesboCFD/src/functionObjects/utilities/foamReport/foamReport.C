/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "foamReport.H"
#include "addToRunTimeSelectionTable.H"
#include "argList.H"
#include "clock.H"
#include "cloud.H"
#include "foamVersion.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "stringOps.H"
#include "substitutionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(foamReport, 0);
    addToRunTimeSelectionTable(functionObject, foamReport, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::foamReport::setStaticBuiltins()
{
    substitutionModel::addBuiltinStr("OF_HOST", Foam::hostName());
    substitutionModel::addBuiltinStr
    (
        "OF_PROC_ZERO_DIR",
        Pstream::parRun() ? "processor0" : ""
    );

    substitutionModel::addBuiltin("OF_API", foamVersion::api);
    substitutionModel::addBuiltinStr("OF_PATCH", foamVersion::patch);
    substitutionModel::addBuiltinStr("OF_BUILD", foamVersion::build);
    substitutionModel::addBuiltinStr("OF_BUILD_ARCH", foamVersion::buildArch);
    substitutionModel::addBuiltinStr("OF_VERSION", foamVersion::version);

    substitutionModel::addBuiltinStr("OF_DATE_START", clock::date());
    substitutionModel::addBuiltinStr("OF_CLOCK_START", clock::clockTime());

    substitutionModel::addBuiltinStr("OF_EXECUTABLE", argList::envExecutable());
    substitutionModel::addBuiltinStr("OF_CASE_PATH", argList::envGlobalPath());
    substitutionModel::addBuiltinStr("OF_CASE_NAME", time().globalCaseName());

    substitutionModel::addBuiltin("OF_NPROCS", Pstream::nProcs());

    // Set mesh builtins when there is only 1 mesh
    const auto meshes = time_.lookupClass<fvMesh>();
    if (meshes.size() == 1)
    {
        const auto& mesh = *(meshes.begin().val());
        substitutionModel::addBuiltin("OF_MESH_NCELLS", mesh.nCells());
        substitutionModel::addBuiltin("OF_MESH_NFACES", mesh.nFaces());
        substitutionModel::addBuiltin("OF_MESH_NEDGES", mesh.nEdges());
        substitutionModel::addBuiltin("OF_MESH_NPOINTS", mesh.nPoints());
        substitutionModel::addBuiltin
        (
            "OF_MESH_NINTERNALFACES",
            mesh.nInternalFaces()
        );
        substitutionModel::addBuiltin
        (
            "OF_MESH_NBOUNDARYFACES",
            mesh.nBoundaryFaces()
        );
        substitutionModel::addBuiltin
        (
            "OF_MESH_NPATCHES",
            mesh.boundaryMesh().nNonProcessor()
        );
        substitutionModel::addBuiltin
        (
            "OF_MESH_BOUNDS_MIN",
            mesh.bounds().min()
        );
        substitutionModel::addBuiltin
        (
            "OF_MESH_BOUNDS_MAX",
            mesh.bounds().max()
        );
    }
}


void Foam::functionObjects::foamReport::setDynamicBuiltins()
{
    // Overwrite existing entries
    substitutionModel::setBuiltinStr("OF_TIME", time().timeName());
    substitutionModel::setBuiltin("OF_NTIMES", time().times().size());
    substitutionModel::setBuiltin("OF_TIME_INDEX", time().timeIndex());
    substitutionModel::setBuiltin("OF_TIME_DELTAT", time().deltaTValue());

    substitutionModel::setBuiltinStr("OF_DATE_NOW", clock::date());
    substitutionModel::setBuiltinStr("OF_CLOCK_NOW", clock::clockTime());

    substitutionModel::setBuiltin("OF_NREGIONS", time().names<fvMesh>().size());
    substitutionModel::setBuiltin("OF_NCLOUDS", time().names<cloud>().size());
}


bool Foam::functionObjects::foamReport::parseTemplate(const fileName& fName)
{
    Info<< "    Reading template from " << fName << endl;

    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Unable to open file " << fName << endl;
    }

    DynamicList<string> contents;
    string buffer;

    label lineNo = 0;
    while (is.good())
    {
        is.getLine(buffer);

        // Collect keys for this line and clean the buffer
        const wordList keys(substitutionModel::getKeys(buffer));

        Tuple2<label, DynamicList<label>> nullValue(-1, DynamicList<label>());

        // Assemble table of keyword and lines where the keyword appears
        for (const word& key : keys)
        {
            if (modelKeys_.insert(key, nullValue))
            {
                // Set substitution model responsible for this keyword
                label modeli = -1;
                forAll(substitutions_, i)
                {
                    if (substitutions_[i].valid(key))
                    {
                        modeli = i;
                        break;
                    }
                }

                // Note: cannot check that key/model is set here
                // - dynamic builtins not ready yet...

                modelKeys_[key].first() = modeli;
            }

            DynamicList<label>& lineNos = modelKeys_[key].second();
            lineNos.push_back(lineNo);
        }

        contents.push_back(buffer);

        ++lineNo;
    }

    templateContents_.transfer(contents);

    return templateContents_.size() > 0;
}


bool Foam::functionObjects::foamReport::apply(Ostream& os) const
{
    List<string> out(templateContents_);

    forAllConstIters(modelKeys_, iter)
    {
        const word& key = iter.key();
        const label modeli = iter.val().first();
        const DynamicList<label>& lineNos = iter.val().second();

        DebugInfo<< "key:" << key << endl;

        for (const label linei : lineNos)
        {
            if (modeli == -1)
            {
                if (!substitutionModel::replaceBuiltin(key, out[linei]))
                {
                    WarningInFunction
                        << "Unable to find substitution for " << key
                        << " on line " << linei << endl;
                }
            }
            else
            {
                substitutions_[modeli].apply(key, out[linei]);
            }
        }
    }

    for (const auto& line : out)
    {
        os  << line.c_str() << nl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::foamReport::foamReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    writeFile(runTime, name, typeName, dict),
    templateFile_(),
    modelKeys_(),
    substitutions_(),
    debugKeys_(dict.getOrDefault<bool>("debugKeys", false))
{
    read(dict);

    setStaticBuiltins();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::foamReport::read(const dictionary& dict)
{
    if (stateFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;

        dict.readEntry("template", templateFile_);

        Info<< "    Template: " << templateFile_ << endl;

        const word ext = templateFile_.ext();

        if (ext.size())
        {
            setExt("." + ext);
        }
        else
        {
            setExt(ext);
        }

        Info<< "    Reading substitutions" << endl;

        const dictionary& subsDict = dict.subDict("substitutions");

        substitutions_.resize(subsDict.size());

        label i = 0;
        for (const entry& e : subsDict)
        {
            if (!e.isDict())
            {
                FatalIOErrorInFunction(subsDict)
                    << "Substitution models must be provided in dictionary "
                    << "format"
                    << exit(FatalIOError);
            }

            substitutions_.set(i++, substitutionModel::New(e.dict(), time()));
        }

        parseTemplate(templateFile_.expand());

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::foamReport::execute()
{
    for (auto& sub : substitutions_)
    {
        sub.update();
    }

    return true;
}


bool Foam::functionObjects::foamReport::write()
{
    if (!Pstream::master()) return true;

    setDynamicBuiltins();

    auto filePtr = newFileAtTime(name(), time().value());
    auto& os = filePtr();

    // Reset stream width (by default assumes fixed width tabular output)
    os.width(0);

    // Perform the substitutions
    apply(os);

    if (debugKeys_)
    {
        os  << "Model keys:" << nl;
        for (const auto& model : substitutions_)
        {
            os  << model.type() << ":" << model.keys() << nl;
        }

        os  << "Builtins:" << nl;
        substitutionModel::writeBuiltins(os);
    }

    return true;
}


// ************************************************************************* //
