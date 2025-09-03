/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    surfaceFeatureExtract

Group
    grpSurfaceUtilities

Description
    Extracts and writes surface features to file. All but the basic feature
    extraction is a work-in-progress.

    The extraction process is driven by the \a system/surfaceFeatureExtractDict
    dictionary, but the \a -dict option can be used to define an alternative
    location.

    The \a system/surfaceFeatureExtractDict dictionary contains entries
    for each extraction process.
    The name of the individual dictionary is used to load the input surface
    (found under \a constant/triSurface) and also as the basename for the
    output.

    If the \c surfaces entry is present in a sub-dictionary, it has absolute
    precedence over a surface name deduced from the dictionary name.
    If the dictionary name itself does not have an extension, the \c surfaces
    entry becomes mandatory since in this case the dictionary name cannot
    represent an input surface file (ie, there is no file extension).
    The \c surfaces entry is a wordRe list, which allows loading and
    combining of multiple surfaces. Any exactly specified surface names must
    exist, but surfaces selected via regular expressions need not exist.
    The selection mechanism preserves order and is without duplicates.
    For example,
    \verbatim
    dictName
    {
        surfaces    (surface1.stl "other.*" othersurf.obj);
        ...
    }
    \endverbatim

    When loading surfaces, the points/faces/regions of each surface are
    normally offset to create an aggregated surface. No merging of points
    or faces is done. The optional entry \c loadingOption can be used to
    adjust the treatment of the regions when loading single or multiple files,
    with selections according to the Foam::triSurfaceLoader::loadingOption
    enumeration.
    \verbatim
    dictName
    {
        // Optional treatment of surface regions when loading
        // (single, file, offset, merge)
        loadingOption   file;
        ...
    }
    \endverbatim
    The \c loadingOption is primarily used in combination with the
    \c intersectionMethod (specifically its \c region option).
    The default \c loadingOption is normally \c offset,
    but this changes to \c file if the \c intersectionMethod
    \c region is being used.

    Once surfaces have been loaded, the first stage is to extract
    features according to the specified \c extractionMethod with values
    as per the following table:
    \table
        extractionMethod   | Description
        none               | No feature extraction
        extractFromFile    | Load features from the file named in featureEdgeFile
        extractFromSurface | Extract features from surface geometry
    \endtable

    There are a few entries that influence the extraction behaviour:
    \verbatim
        // File to use for extractFromFile input
        featureEdgeFile     "FileName"

        // Mark edges whose adjacent surface normals are at an angle less
        // than includedAngle as features
        // - 0  : selects no edges
        // - 180: selects all edges
        includedAngle       120;

        // Do not mark region edges
        geometricTestOnly   yes;
    \endverbatim

    This initial set of edges can be trimmed:
    \verbatim
        trimFeatures
        {
            // Remove features with fewer than the specified number of edges
            minElem         0;

            // Remove features shorter than the specified cumulative length
            minLen          0.0;
        }
    \endverbatim

    and subsetted
    \verbatim
    subsetFeatures
    {
        // Use a plane to select feature edges (normal)(basePoint)
        // Only keep edges that intersect the plane
        plane           (1 0 0)(0 0 0);

        // Select feature edges using a box // (minPt)(maxPt)
        // Only keep edges inside the box:
        insideBox       (0 0 0)(1 1 1);

        // Only keep edges outside the box:
        outsideBox      (0 0 0)(1 1 1);

        // Keep nonManifold edges (edges with >2 connected faces where
        // the faces form more than two different normal planes)
        nonManifoldEdges yes;

        // Keep open edges (edges with 1 connected face)
        openEdges       yes;
    }
    \endverbatim

    Subsequently, additional features can be added from another file:
    \verbatim
        addFeatures
        {
            // Add (without merging) another extendedFeatureEdgeMesh
            name        axZ.extendedFeatureEdgeMesh;
        }
    \endverbatim

    The intersectionMethod provides a final means of adding additional
    features. These are loosely termed "self-intersection", since it
    detects the face/face intersections of the loaded surface or surfaces.

    \table
        intersectionMethod | Description
        none    | Do nothing
        self    | All face/face intersections
        region  | Limit face/face intersections to those between different regions.
    \endtable
    The optional \c tolerance tuning parameter is available for handling
    the face/face intersections, but should normally not be touched.

    As well as the normal extendedFeatureEdgeMesh written,
    other items can be selected with boolean switches:

    \table
        Output option | Description
        closeness | Output the closeness of surface elements to other surface elements.
        curvature | Output surface curvature
        featureProximity | Output the proximity of feature points and edges to another
        writeObj  | Write features to OBJ format for postprocessing
        writeVTK  | Write closeness/curvature/proximity fields as VTK for postprocessing
    \endtable

Note
   Although surfaceFeatureExtract can do many things, there are still a fair
   number of corner cases where it may not produce the desired result.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "edgeMeshTools.H"
#include "surfaceFeaturesExtraction.H"
#include "surfaceIntersection.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "triSurfaceMesh.H"
#include "foamVtkSurfaceWriter.H"
#include "unitConversion.H"
#include "plane.H"
#include "point.H"
#include "triSurfaceLoader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract and write surface feature lines to file.\n"
        "Feature line extraction only valid on closed manifold surfaces."
    );

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    argList::addOption
    (
        "dict",
        "file",
        "Read surfaceFeatureExtractDict from specified location"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< nl
        << "Note: "
        << "Feature line extraction only valid on closed manifold surfaces"
        << nl << nl;

    const word dictName("surfaceFeatureExtractDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;
    const IOdictionary dict(dictIO);

    // Loader for available triSurface surface files
    triSurfaceLoader loader(runTime);

    // Where to write VTK output files
    const fileName vtkOutputDir = runTime.constantPath()/"triSurface";

    for (const entry& dEntry : dict)
    {
        if (!dEntry.isDict())  // dictionary entries only
        {
            continue;
        }

        // Note - do not check for shell meta-characters
        // - user responsibility

        const word& dictName = dEntry.keyword();
        const dictionary& surfaceDict = dEntry.dict();

        if (!surfaceDict.found("extractionMethod"))
        {
            // Insist on an extractionMethod
            continue;
        }

        // The output name based in dictionary name (without extensions)
        const word outputName = dictName.lessExt();

        autoPtr<surfaceFeaturesExtraction::method> extractor =
            surfaceFeaturesExtraction::method::New
            (
                surfaceDict
            );

        // We don't needs the intersectionMethod yet, but can use it
        // for setting a reasonable loading option
        const surfaceIntersection::intersectionType selfIntersect =
            surfaceIntersection::selfIntersectionNames.getOrDefault
            (
                "intersectionMethod",
                surfaceDict,
                surfaceIntersection::NONE
            );

        const Switch writeObj("writeObj", surfaceDict, Switch::OFF);

        const Switch writeVTK("writeVTK", surfaceDict, Switch::OFF);

        // The "surfaces" entry is normally optional, but make it mandatory
        // if the dictionary name doesn't have an extension
        // (ie, probably not a surface filename at all).
        // If it is missing, this will fail nicely with an appropriate error
        // message.
        if (surfaceDict.found("surfaces") || !dictName.has_ext())
        {
            loader.select(surfaceDict.get<wordRes>("surfaces"));
        }
        else
        {
            loader.select(dictName);
        }

        // DebugVar(loader.available());
        // DebugVar(outputName);

        if (loader.selected().empty())
        {
            FatalErrorInFunction
                << "No surfaces specified/found for entry: "
                << dictName << exit(FatalError);
        }

        Info<< "Surfaces     : ";
        if (loader.selected().size() == 1)
        {
            Info<< loader.selected().first() << nl;
        }
        else
        {
            Info<< flatOutput(loader.selected()) << nl;
        }
        Info<< "Output       : " << outputName << nl;

        // Loading option - default depends on context
        triSurfaceLoader::loadingOption loadingOption =
            triSurfaceLoader::loadingOptionNames.getOrDefault
            (
                "loadingOption",
                surfaceDict,
                (
                    selfIntersect == surfaceIntersection::SELF_REGION
                  ? triSurfaceLoader::FILE_REGION
                  : triSurfaceLoader::OFFSET_REGION
                )
            );

        Info<<"Load options : "
            << triSurfaceLoader::loadingOptionNames[loadingOption] << nl
            << "Write options:"
            << " writeObj=" << writeObj
            << " writeVTK=" << writeVTK << nl;

        scalar scaleFactor = -1;
        // Allow rescaling of the surface points (eg, mm -> m)
        if (surfaceDict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
        {
            Info<<"Scaling : " << scaleFactor << nl;
        }

        // Load a single file, or load and combine multiple selected files
        autoPtr<triSurface> surfPtr = loader.load(loadingOption, scaleFactor);
        if (!surfPtr || surfPtr->empty())
        {
            FatalErrorInFunction
                << "Problem loading surface(s) for entry: "
                << dictName << exit(FatalError);
        }

        triSurface surf = *surfPtr;

        Info<< nl
            << "Statistics:" << nl;
        surf.writeStats(Info);

        // Need a copy as plain faces if outputting VTK format
        faceList faces;
        if (writeVTK)
        {
            faces.setSize(surf.size());
            forAll(surf, fi)
            {
                faces[fi] = surf[fi];
            }
        }


        const dictionary* subDictPtr = nullptr;

        //
        // Extract features using the preferred extraction method
        //
        autoPtr<surfaceFeatures> features = extractor().features(surf);

        // Trim set
        // ~~~~~~~~

        // Option: "trimFeatures" (dictionary)
        if ((subDictPtr = surfaceDict.findDict("trimFeatures")) != nullptr)
        {
            const dictionary& trimDict = *subDictPtr;

            const scalar minLen =
                trimDict.getOrDefault<scalar>("minLen", 0);
            const label minElem =
                trimDict.getOrDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minLen > 0 || minElem > 0)
            {
                if (minLen > 0)
                {
                    Info<< "Removing features of length < "
                        << minLen << endl;
                }
                if (minElem > 0)
                {
                    Info<< "Removing features with number of edges < "
                        << minElem << endl;
                }

                features().trimFeatures
                (
                    minLen, minElem, extractor().includedAngle()
                );
            }
        }

        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(features().toStatus());

        // Option: "subsetFeatures" (dictionary)
        if ((subDictPtr = surfaceDict.findDict("subsetFeatures")) != nullptr)
        {
            const dictionary& subsetDict = *subDictPtr;

            treeBoundBox bb;

            // Suboption: "insideBox"
            if (subsetDict.readIfPresent("insideBox", bb))
            {
                Info<< "Subset edges inside box " << bb << endl;
                features().subsetBox(edgeStat, bb);

                {
                    OBJstream os("subsetBox.obj");

                    Info<< "Dumping bounding box " << bb
                        << " as lines to obj file "
                        << os.name() << endl;

                    os.write(bb);
                }
            }
            // Suboption: "outsideBox"
            else if (subsetDict.readIfPresent("outsideBox", bb))
            {
                Info<< "Exclude edges outside box " << bb << endl;
                features().excludeBox(edgeStat, bb);

                {
                    OBJstream os("deleteBox.obj");

                    Info<< "Dumping bounding box " << bb
                        << " as lines to obj file "
                        << os.name() << endl;

                    os.write(bb);
                }
            }

            // Suboption: "nonManifoldEdges" (false: remove non-manifold edges)
            if (!subsetDict.getOrDefault("nonManifoldEdges", true))
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                features().checkFlatRegionEdge
                (
                    edgeStat,
                    1e-5,   // tol
                    extractor().includedAngle()
                );
            }

            if (!subsetDict.getOrDefault("strictNonManifoldEdges", true))
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces)"
                    << endl;

                features().excludeNonManifold(edgeStat);
            }

            // Suboption: "openEdges" (false: remove open edges)
            if (!subsetDict.getOrDefault("openEdges", true))
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                features().excludeOpen(edgeStat);
            }

            // Suboption: "plane"
            if (subsetDict.found("plane"))
            {
                plane cutPlane(subsetDict.lookup("plane"));

                Info<< "Only include feature edges that intersect the plane"
                    << " with normal " << cutPlane.normal()
                    << " and origin " << cutPlane.origin() << endl;

                features().subsetPlane(edgeStat, cutPlane);
            }
        }

        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat, extractor().includedAngle());

        Info<< nl << "Initial ";
        newSet.writeStats(Info);

        boolList surfBaffleRegions(surf.patches().size(), false);
        if (surfaceDict.found("baffles"))
        {
            wordRes baffleSelect(surfaceDict.get<wordRes>("baffles"));

            wordList patchNames(surf.patches().size());
            forAll(surf.patches(), patchi)
            {
                patchNames[patchi] = surf.patches()[patchi].name();
            }

            labelList indices(baffleSelect.matching(patchNames));

            for (const label patchId : indices)
            {
                surfBaffleRegions[patchId] = true;
            }

            if (indices.size())
            {
                Info<< "Adding " << indices.size() << " baffle regions: (";

                forAll(surfBaffleRegions, patchi)
                {
                    if (surfBaffleRegions[patchi])
                    {
                        Info<< ' ' << patchNames[patchi];
                    }
                }
                Info<< " )" << nl << nl;
            }
        }

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            outputName + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        );


        if ((subDictPtr = surfaceDict.findDict("addFeatures")) != nullptr)
        {
            const dictionary& addFeaturesDict = *subDictPtr;

            const word addFeName = addFeaturesDict.get<word>("name");

            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            edgeMeshTools::writeStats(Info, addFeMesh);

            feMesh.add(addFeMesh);
        }

        if (selfIntersect != surfaceIntersection::NONE)
        {
            triSurfaceSearch query(surf);
            surfaceIntersection intersect(query, surfaceDict);

            // Remove rounding noise. For consistency could use 1e-6,
            // as per extractFromFile implementation

            intersect.mergePoints(10*SMALL);

            labelPair sizeInfo
            (
                intersect.cutPoints().size(),
                intersect.cutEdges().size()
            );

            if (intersect.cutEdges().size())
            {
                extendedEdgeMesh addMesh
                (
                    intersect.cutPoints(),
                    intersect.cutEdges()
                );

                feMesh.add(addMesh);

                sizeInfo[0] = addMesh.points().size();
                sizeInfo[1] = addMesh.edges().size();
            }
            Info<< nl
                << "intersection: "
                << surfaceIntersection::selfIntersectionNames[selfIntersect]
                << nl
                << "    points : " << sizeInfo[0] << nl
                << "    edges  : " << sizeInfo[1] << nl;
        }

        Info<< nl << "Final ";
        edgeMeshTools::writeStats(Info, feMesh);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj(feMesh.path()/outputName);
        }

        feMesh.write();

        // Write a featureEdgeMesh (.eMesh) for backwards compatibility
        // Used by snappyHexMesh (JUN-2017)
        if (true)
        {
            featureEdgeMesh bfeMesh
            (
                IOobject
                (
                    outputName + ".eMesh",      // name
                    runTime.constant(),         // instance
                    "triSurface",
                    runTime,                    // registry
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                feMesh.points(),
                feMesh.edges()
            );

            Info<< nl << "Writing featureEdgeMesh to "
                << bfeMesh.objectPath() << endl;

            bfeMesh.regIOobject::write();
        }


        // Output information

        const bool optCloseness =
            surfaceDict.getOrDefault("closeness", false);

        const bool optProximity =
            surfaceDict.getOrDefault("featureProximity", false);

        const bool optCurvature =
            surfaceDict.getOrDefault("curvature", false);


        // For VTK legacy format, we would need an a priori count of
        // CellData and PointData fields.
        // For convenience, we therefore only use the XML formats

        autoPtr<vtk::surfaceWriter> vtkWriter;

        if (optCloseness || optProximity || optCurvature)
        {
            if (writeVTK)
            {
                vtkWriter.reset
                (
                    new vtk::surfaceWriter
                    (
                        surf.points(),
                        faces,
                        (vtkOutputDir / outputName),
                        false  // serial only
                    )
                );

                vtkWriter->writeGeometry();

                Info<< "Writing VTK to "
                    << runTime.relativePath(vtkWriter->output()) << nl;
            }
        }
        else
        {
            continue;  // Nothing to output
        }


        // Option: "closeness"
        if (optCloseness)
        {
            Pair<tmp<scalarField>> tcloseness =
                triSurfaceTools::writeCloseness
                (
                    runTime,
                    outputName,
                    surf,
                    45,  // internalAngleTolerance
                    10   // externalAngleTolerance
                );

            if (vtkWriter)
            {
                vtkWriter->beginCellData();
                vtkWriter->write("internalCloseness", tcloseness[0]());
                vtkWriter->write("externalCloseness", tcloseness[1]());
            }
        }

        // Option: "featureProximity"
        if (optCloseness)
        {
            const scalar maxProximity =
                surfaceDict.getOrDefault<scalar>("maxFeatureProximity", 1);

            tmp<scalarField> tproximity =
                edgeMeshTools::writeFeatureProximity
                (
                    runTime,
                    outputName,
                    feMesh,
                    surf,
                    maxProximity
                );

            if (vtkWriter)
            {
                vtkWriter->beginCellData();
                vtkWriter->write("featureProximity", tproximity());
            }
        }

        // Option: "curvature"
        if (optCurvature)
        {
            tmp<scalarField> tcurvature =
                triSurfaceTools::writeCurvature
                (
                    runTime,
                    outputName,
                    surf
                );

            if (vtkWriter)
            {
                vtkWriter->beginPointData();
                vtkWriter->write("curvature", tcurvature());
            }
        }

        Info<< endl;
    }

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
