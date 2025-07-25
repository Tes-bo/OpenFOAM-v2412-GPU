#include "writeFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "polyMeshTools.H"
#include "syncTools.H"
#include "tetrahedron.H"
#include "regionSplit.H"
#include "wallDist.H"
#include "cellAspectRatio.H"

using namespace Foam;

void maxFaceToCell
(
    const scalarField& faceData,
    volScalarField& cellData
)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = -GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = max(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void minFaceToCell
(
    const scalarField& faceData,
    volScalarField& cellData
)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = min(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void minFaceToCell
(
    const surfaceScalarField& faceData,
    volScalarField& cellData,
    const bool correctBoundaryConditions
)
{
    scalarField& cellFld = cellData.ref();

    cellFld = GREAT;

    const labelUList& own = cellData.mesh().owner();
    const labelUList& nei = cellData.mesh().neighbour();

    // Internal faces
    forAll(own, facei)
    {
        cellFld[own[facei]] = min(cellFld[own[facei]], faceData[facei]);
        cellFld[nei[facei]] = min(cellFld[nei[facei]], faceData[facei]);
    }

    // Patch faces
    forAll(faceData.boundaryField(), patchi)
    {
        const fvsPatchScalarField& fvp = faceData.boundaryField()[patchi];
        const labelUList& fc = fvp.patch().faceCells();

        forAll(fc, i)
        {
            cellFld[fc[i]] = min(cellFld[fc[i]], fvp[i]);
        }
    }

    volScalarField::Boundary& bfld = cellData.boundaryFieldRef();

    forAll(bfld, patchi)
    {
        bfld[patchi] = faceData.boundaryField()[patchi];
    }
    if (correctBoundaryConditions)
    {
        cellData.correctBoundaryConditions();
    }
}


void writeSurfaceField
(
    const fvMesh& mesh,
    const fileName& fName,
    const scalarField& faceData
)
{
    // Write single surfaceScalarField

    surfaceScalarField fld
    (
        IOobject
        (
            fName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    fld.primitiveFieldRef() = faceData;
    //fld.correctBoundaryConditions();
    Info<< "    Writing face data to " << fName << endl;
    fld.write();
}


void Foam::writeFields
(
    const fvMesh& mesh,
    const wordHashSet& selectedFields,
    const bool writeFaceFields
)
{
    if (selectedFields.empty())
    {
        return;
    }

    Info<< "Writing fields with mesh quality parameters" << endl;

    if (selectedFields.found("nonOrthoAngle"))
    {
        //- Face based orthogonality
        const scalarField faceOrthogonality
        (
            polyMeshTools::faceOrthogonality
            (
                mesh,
                mesh.faceAreas(),
                mesh.cellCentres()
            )
        );

        //- Face based angle
        const scalarField nonOrthoAngle
        (
            radToDeg
            (
                Foam::acos(min(scalar(1), max(scalar(-1), faceOrthogonality)))
            )
        );

        //- Cell field - max of either face
        volScalarField cellNonOrthoAngle
        (
            IOobject
            (
                "nonOrthoAngle",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        //- Take max
        maxFaceToCell(nonOrthoAngle, cellNonOrthoAngle);
        Info<< "    Writing non-orthogonality (angle) to "
            << cellNonOrthoAngle.name() << endl;
        cellNonOrthoAngle.write();

        if (writeFaceFields)
        {
            writeSurfaceField
            (
                mesh,
                "face_nonOrthoAngle",
                SubField<scalar>(nonOrthoAngle, mesh.nInternalFaces())
            );
        }
    }

    if (selectedFields.found("faceWeight"))
    {
        volScalarField cellWeights
        (
            IOobject
            (
                "faceWeight",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            wordList                                // wanted bc types
            (
                mesh.boundary().size(),
                fvPatchFieldBase::calculatedType()
            ),
            mesh.weights().boundaryField().types()  // current bc types
        );
        //- Take min
        minFaceToCell(mesh.weights(), cellWeights, false);
        Info<< "    Writing face interpolation weights (0..0.5) to "
            << cellWeights.name() << endl;
        cellWeights.write();

        if (writeFaceFields)
        {
            writeSurfaceField(mesh, "face_faceWeight", mesh.weights());
        }
    }


    // Skewness
    // ~~~~~~~~

    if (selectedFields.found("skewness"))
    {
        //- Face based skewness
        const scalarField faceSkewness
        (
            polyMeshTools::faceSkewness
            (
                mesh,
                mesh.points(),
                mesh.faceCentres(),
                mesh.faceAreas(),
                mesh.cellCentres()
            )
        );

        //- Cell field - max of either face
        volScalarField cellSkewness
        (
            IOobject
            (
                "skewness",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        //- Take max
        maxFaceToCell(faceSkewness, cellSkewness);
        Info<< "    Writing face skewness to " << cellSkewness.name() << endl;
        cellSkewness.write();

        if (writeFaceFields)
        {
            writeSurfaceField
            (
                mesh,
                "face_skewness",
                SubField<scalar>(faceSkewness, mesh.nInternalFaces())
            );
        }
    }


    // cellDeterminant
    // ~~~~~~~~~~~~~~~

    if (selectedFields.found("cellDeterminant"))
    {
        volScalarField cellDeterminant
        (
            IOobject
            (
                "cellDeterminant",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        cellDeterminant.primitiveFieldRef() =
            primitiveMeshTools::cellDeterminant
            (
                mesh,
                mesh.geometricD(),
                mesh.faceAreas(),
                syncTools::getInternalOrCoupledFaces(mesh)
            );
        cellDeterminant.correctBoundaryConditions();
        Info<< "    Writing cell determinant to "
            << cellDeterminant.name() << endl;
        cellDeterminant.write();
    }


    // Aspect ratio
    // ~~~~~~~~~~~~

    if (selectedFields.found("aspectRatio"))
    {
        volScalarField aspectRatio
        (
            IOobject
            (
                "aspectRatio",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );


        scalarField cellOpenness;
        polyMeshTools::cellClosedness
        (
            mesh,
            mesh.geometricD(),
            mesh.faceAreas(),
            mesh.cellVolumes(),
            cellOpenness,
            aspectRatio.ref()
        );

        aspectRatio.correctBoundaryConditions();
        Info<< "    Writing aspect ratio to " << aspectRatio.name() << endl;
        aspectRatio.write();
    }

    if (selectedFields.found("cellAspectRatio"))
    {
        volScalarField aspectRatio
        (
            IOobject
            (
                "cellAspectRatio",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        aspectRatio.ref().field() = cellAspectRatio(mesh);

        aspectRatio.correctBoundaryConditions();
        Info<< "    Writing approximate aspect ratio to "
            << aspectRatio.name() << endl;
        aspectRatio.write();
    }


    // cell type
    // ~~~~~~~~~

    if (selectedFields.found("cellShapes"))
    {
        volScalarField shape
        (
            IOobject
            (
                "cellShapes",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        const cellShapeList& cellShapes = mesh.cellShapes();
        forAll(cellShapes, cellI)
        {
            const cellModel& model = cellShapes[cellI].model();
            shape[cellI] = model.index();
        }
        shape.correctBoundaryConditions();
        Info<< "    Writing cell shape (hex, tet etc.) to " << shape.name()
            << endl;
        shape.write();
    }

    if (selectedFields.found("cellVolume"))
    {
        volScalarField V
        (
            IOobject
            (
                "cellVolume",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimVolume, Zero)
        );
        V.ref() = mesh.V();
        Info<< "    Writing cell volume to " << V.name() << endl;
        V.write();
    }

    if (selectedFields.found("cellVolumeRatio"))
    {
        const scalarField faceVolumeRatio
        (
            polyMeshTools::volRatio
            (
                mesh,
                mesh.V()
            )
        );

        volScalarField cellVolumeRatio
        (
            IOobject
            (
                "cellVolumeRatio",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        //- Take min
        minFaceToCell(faceVolumeRatio, cellVolumeRatio);
        Info<< "    Writing cell volume ratio to "
            << cellVolumeRatio.name() << endl;
        cellVolumeRatio.write();

        if (writeFaceFields)
        {
            writeSurfaceField
            (
                mesh,
                "face_cellVolumeRatio",
                SubField<scalar>(faceVolumeRatio, mesh.nInternalFaces())
            );
        }
    }

    // minTetVolume
    if (selectedFields.found("minTetVolume"))
    {
        volScalarField minTetVolume
        (
            IOobject
            (
                "minTetVolume",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar("minTetVolume", dimless, GREAT),
            fvPatchFieldBase::zeroGradientType()
        );


        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();
        const pointField& p = mesh.points();
        forAll(own, facei)
        {
            const face& f = mesh.faces()[facei];
            const point& fc = mesh.faceCentres()[facei];

            {
                const point& ownCc = mesh.cellCentres()[own[facei]];
                scalar& ownVol = minTetVolume[own[facei]];
                forAll(f, fp)
                {
                    scalar tetQual = tetPointRef
                    (
                        p[f[fp]],
                        p[f.nextLabel(fp)],
                        ownCc,
                        fc
                    ).quality();
                    ownVol = min(ownVol, tetQual);
                }
            }
            if (mesh.isInternalFace(facei))
            {
                const point& neiCc = mesh.cellCentres()[nei[facei]];
                scalar& neiVol = minTetVolume[nei[facei]];
                forAll(f, fp)
                {
                    scalar tetQual = tetPointRef
                    (
                        p[f[fp]],
                        p[f.nextLabel(fp)],
                        fc,
                        neiCc
                    ).quality();
                    neiVol = min(neiVol, tetQual);
                }
            }
        }

        minTetVolume.correctBoundaryConditions();
        Info<< "    Writing minTetVolume to " << minTetVolume.name() << endl;
        minTetVolume.write();
    }

    // minPyrVolume
    if (selectedFields.found("minPyrVolume"))
    {
        volScalarField minPyrVolume
        (
            IOobject
            (
                "minPyrVolume",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar("minPyrVolume", dimless, GREAT),
            fvPatchFieldBase::zeroGradientType()
        );

        // Get owner and neighbour pyr volumes
        scalarField ownPyrVol(mesh.nFaces());
        scalarField neiPyrVol(mesh.nInternalFaces());
        primitiveMeshTools::facePyramidVolume
        (
            mesh,
            mesh.points(),
            mesh.cellCentres(),

            ownPyrVol,
            neiPyrVol
        );

        // Get min pyr vol per cell
        scalarField& cellFld = minPyrVolume.ref();
        cellFld = GREAT;

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        // Internal faces
        forAll(own, facei)
        {
            cellFld[own[facei]] = min(cellFld[own[facei]], ownPyrVol[facei]);
            cellFld[nei[facei]] = min(cellFld[nei[facei]], neiPyrVol[facei]);
        }

        // Patch faces
        for (const auto& fvp : minPyrVolume.boundaryField())
        {
            const labelUList& fc = fvp.patch().faceCells();

            forAll(fc, i)
            {
                const label meshFacei = fvp.patch().start();
                cellFld[fc[i]] = min(cellFld[fc[i]], ownPyrVol[meshFacei]);
            }
        }

        minPyrVolume.correctBoundaryConditions();
        Info<< "    Writing minPyrVolume to " << minPyrVolume.name() << endl;
        minPyrVolume.write();

        if (writeFaceFields)
        {
            scalarField minFacePyrVol(neiPyrVol);
            minFacePyrVol = min
            (
                minFacePyrVol,
                SubField<scalar>(ownPyrVol, mesh.nInternalFaces())
            );
            writeSurfaceField(mesh, "face_minPyrVolume", minFacePyrVol);
        }
    }

    if (selectedFields.found("cellRegion"))
    {
        volScalarField cellRegion
        (
            IOobject
            (
                "cellRegion",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        regionSplit rs(mesh);
        forAll(rs, celli)
        {
            cellRegion[celli] = rs[celli];
        }
        cellRegion.correctBoundaryConditions();
        Info<< "    Writing cell region to " << cellRegion.name() << endl;
        cellRegion.write();
    }

    if (selectedFields.found("wallDistance"))
    {
        // See if wallDist.method entry in fvSchemes before calling factory
        // method of wallDist. Have 'failing' version of wallDist::New instead?
        const dictionary& schemesDict =
            static_cast<const fvSchemes&>(mesh).schemesDict();
        if (schemesDict.found("wallDist"))
        {
            if (schemesDict.subDict("wallDist").found("method"))
            {
                // Wall distance
                volScalarField y("wallDistance", wallDist::New(mesh).y());
                Info<< "    Writing wall distance to " << y.name() << endl;
                y.write();

                // Wall-reflection vectors
                //const volVectorField& n = wallDist::New(mesh).n();
                //Info<< "    Writing wall normal to " << n.name() << endl;
                //n.write();
            }
        }
    }

    if (selectedFields.found("cellZone"))
    {
        volScalarField cellZone
        (
            IOobject
            (
                "cellZone",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(word::null, dimless, -1)
        );

        const cellZoneMesh& czs = mesh.cellZones();
        for (const auto& zone : czs)
        {
            UIndirectList<scalar>(cellZone, zone) = zone.index();
        }

        cellZone.correctBoundaryConditions();
        Info<< "    Writing cell zoning to " << cellZone.name() << endl;
        cellZone.write();
    }
    if (selectedFields.found("faceZone"))
    {
        // Determine for each face the zone index (scalar for ease of
        // manipulation)
        scalarField zoneID(mesh.nFaces(), -1);
        const faceZoneMesh& czs = mesh.faceZones();
        for (const auto& zone : czs)
        {
            UIndirectList<scalar>(zoneID, zone) = zone.index();
        }


        // Split into internal and boundary values
        surfaceScalarField faceZone
        (
            IOobject
            (
                "faceZone",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(word::null, dimless, -1)
        );

        faceZone.primitiveFieldRef() =
            SubField<scalar>(zoneID, mesh.nInternalFaces());
        surfaceScalarField::Boundary& bfld = faceZone.boundaryFieldRef();
        for (auto& pfld : bfld)
        {
            const fvPatch& fvp = pfld.patch();
            pfld == SubField<scalar>(zoneID, fvp.size(), fvp.start());
        }

        //faceZone.correctBoundaryConditions();
        Info<< "    Writing face zoning to " << faceZone.name() << endl;
        faceZone.write();
    }

    Info<< endl;
}
