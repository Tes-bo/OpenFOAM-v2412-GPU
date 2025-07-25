/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2017-2020 German Aerospace Center (DLR)
    Copyright (C) 2020 Johan Roenby
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
    setAlphaField

Description
    Uses cutCellIso to create a volume fraction field from either a cylinder,
    a sphere or a plane.

    Original code supplied by Johan Roenby, DHI (2016)
    Modification Henning Scheufler, DLR (2019)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"

#include "triSurface.H"
#include "triSurfaceTools.H"

#include "implicitFunction.H"

#include "cutCellIso.H"
#include "cutFaceIso.H"
#include "OBJstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void isoFacesToFile
(
    const DynamicList<List<point>>& faces,
    const word& fileName
)
{
    // Writing isofaces to OBJ file for inspection in paraview
    OBJstream os(fileName + ".obj");

    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            Info<< "Writing file: " << fileName << endl;

            for (const DynamicList<List<point>>& procFaces : allProcFaces)
            {
                for (const List<point>& facePts : procFaces)
                {
                    os.writeFace(facePts);
                }
            }
        }
    }
    else
    {
        Info<< "Writing file: " << fileName << endl;

        for (const List<point>& facePts : faces)
        {
            os.writeFace(facePts);
        }
    }
}

void setAlpha
(
    volScalarField& alpha1,
    DynamicList<List<point>>& facePts,
    scalarField& f,
    const bool writeOBJ
)
{
    const fvMesh& mesh = alpha1.mesh();
    cutCellIso cutCell(mesh, f);
    cutFaceIso cutFace(mesh, f);

    forAll(alpha1, cellI)
    {
        cutCell.calcSubCell(cellI, 0.0);

        alpha1[cellI] = clamp(cutCell.VolumeOfFluid(), zero_one{});

        if (writeOBJ && (mag(cutCell.faceArea()) >= 1e-14))
        {
            facePts.append(cutCell.facePoints());
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].size() > 0)
        {
            const label start = mesh.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                cutFace.calcSubFace(facei, 0.0);
                alphap[patchFacei] =
                    mag(cutFace.subFaceArea())/magSfp[patchFacei];
            }
        }
    }
}


int main(int argc, char *argv[])
{
    argList::noFunctionObjects();           // Never use function objects

    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addNote
    (
        "Uses cutCellIso to create a volume fraction field from an "
        "implicit function."
    );

    argList::addOption
    (
        "dict",
        "file",
        "Alternative setAlphaFieldDict dictionary"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Set time from specified time options, or no-op
    timeSelector::setTimeIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const word dictName("setAlphaFieldDict");
    #include "setSystemMeshDictionaryIO.H"

    IOdictionary setAlphaFieldDict(dictIO);

    Info<< "Reading " << setAlphaFieldDict.name() << endl;

    const word fieldName = setAlphaFieldDict.get<word>("field");
    const bool invert = setAlphaFieldDict.getOrDefault("invert", false);
    const bool writeOBJ = setAlphaFieldDict.getOrDefault("writeOBJ", false);

    Info<< "Reading field " << fieldName << nl << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    autoPtr<implicitFunction> func = implicitFunction::New
    (
        setAlphaFieldDict.get<word>("type"),
        setAlphaFieldDict
    );

    scalarField f(mesh.nPoints(), Zero);

    forAll(f, pi)
    {
        f[pi] = func->value(mesh.points()[pi]);
    };

    DynamicList<List<point>> facePts;

    setAlpha(alpha1, facePts, f, writeOBJ);

    if (writeOBJ)
    {
        isoFacesToFile(facePts, fieldName + "0");
    }

    ISstream::defaultPrecision(18);

    if (invert)
    {
        alpha1 = scalar(1) - alpha1;
    }

    Info<< "Writing new alpha field " << alpha1.name() << endl;
    alpha1.write();

    const scalarField& alpha = alpha1.internalField();

    Info<< "sum(alpha*V):" << gSum(mesh.V()*alpha)
        << ", 1-max(alpha1): " << 1 - gMax(alpha)
        << " min(alpha1): " << gMin(alpha) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
