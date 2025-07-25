/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "nearWallFields.H"
#include "wordRes.H"
#include "findCellParticle.H"
#include "mappedPatchBase.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(nearWallFields, 0);
    addToRunTimeSelectionTable(functionObject, nearWallFields, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::nearWallFields::calcAddressing()
{
    // Count number of faces
    label nPatchFaces = 0;
    for (const label patchi : patchIDs_)
    {
        nPatchFaces += mesh_.boundary()[patchi].size();
    }

    // Global indexing
    const globalIndex globalWalls(nPatchFaces);

    DebugInFunction << "nPatchFaces: " << globalWalls.totalSize() << endl;

    // Start with empty cloud
    Cloud<findCellParticle> cloud(mesh_, Foam::zero{}, cloud::defaultName);

    // Add particles to track to sample locations
    nPatchFaces = 0;

    for (const label patchi : patchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const vectorField nf(patch.nf());
        const vectorField faceCellCentres(patch.patch().faceCellCentres());
        const labelUList& faceCells = patch.patch().faceCells();
        const vectorField::subField faceCentres = patch.patch().faceCentres();

        forAll(patch, patchFacei)
        {
            label meshFacei = patch.start()+patchFacei;

            // Find starting point on face (since faceCentre might not
            // be on face-diagonal decomposition)
            pointIndexHit startInfo
            (
                mappedPatchBase::facePoint
                (
                    mesh_,
                    meshFacei,
                    polyMesh::FACE_DIAG_TRIS
                )
            );


            // Starting point and tet
            point start;
            const label celli = faceCells[patchFacei];

            if (startInfo.hit())
            {
                // Move start point slightly in so it is inside the tet
                const face& f = mesh_.faces()[meshFacei];

                label tetFacei = meshFacei;
                label tetPti = (startInfo.index()+1) % f.size();

                start = startInfo.point();

                // Uncomment below to shift slightly in:
                tetIndices tet(celli, tetFacei, tetPti);
                start =
                    (1.0 - 1e-6)*startInfo.point()
                  + 1e-6*tet.tet(mesh_).centre();

                // Re-check that we have a valid location
                mesh_.findTetFacePt(celli, start, tetFacei, tetPti);
                if (tetFacei == -1)
                {
                    start = faceCellCentres[patchFacei];
                }
            }
            else
            {
                // Fallback: start tracking from neighbouring cell centre
                start = faceCellCentres[patchFacei];
            }


            // TBD: why use start? and not faceCentres[patchFacei]
            //const point end = start-distance_*nf[patchFacei];
            const point end = faceCentres[patchFacei]-distance_*nf[patchFacei];


            // Add a particle to the cloud with originating face as passive data
            cloud.addParticle
            (
                new findCellParticle
                (
                    mesh_,
                    start,
                    celli,
                    end,
                    globalWalls.toGlobal(nPatchFaces) // passive data
                )
            );

            nPatchFaces++;
        }
    }



    if (debug)
    {
        // Dump particles
        OBJstream str
        (
            mesh_.time().path()
           /"wantedTracks_" + mesh_.time().timeName() + ".obj"
        );
        InfoInFunction << "Dumping tracks to " << str.name() << endl;

        for (const findCellParticle& tp : cloud)
        {
            str.writeLine(tp.position(), tp.end());
        }
    }



    // Per cell: empty or global wall index and end location
    cellToWalls_.setSize(mesh_.nCells());
    cellToSamples_.setSize(mesh_.nCells());

    // Database to pass into findCellParticle::move
    findCellParticle::trackingData td(cloud, cellToWalls_, cellToSamples_);

    // Track all particles to their end position.
    scalar maxTrackLen = 2.0*mesh_.bounds().mag();


    //Debug: collect start points
    pointField start;
    if (debug)
    {
        start.setSize(nPatchFaces);
        nPatchFaces = 0;
        for (const findCellParticle& tp : cloud)
        {
            start[nPatchFaces++] = tp.position();
        }
    }


    cloud.move(cloud, td, maxTrackLen);


    // Rework cell-to-globalpatchface into a map
    List<Map<label>> compactMap;
    getPatchDataMapPtr_.reset
    (
        new mapDistribute
        (
            globalWalls,
            cellToWalls_,
            compactMap
        )
    );


    // Debug: dump resulting tracks
    if (debug)
    {
        getPatchDataMapPtr_().distribute(start);
        {
            OBJstream str
            (
                mesh_.time().path()
               /"obtainedTracks_" + mesh_.time().timeName() + ".obj"
            );
            InfoInFunction << "Dumping obtained to " << str.name() << endl;

            forAll(cellToWalls_, celli)
            {
                const List<point>& ends = cellToSamples_[celli];
                const labelList& cData = cellToWalls_[celli];
                forAll(cData, i)
                {
                    str.writeLine(ends[i], start[cData[i]]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::nearWallFields::nearWallFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::nearWallFields::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    fvMeshFunctionObject::read(dict);

    dict.readEntry("fields", fieldSet_);
    dict.readEntry("distance", distance_);

    // Can also use pbm.indices(), but no warnings...
    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();


    // Clear out any previously loaded fields
    vsf_.clear();
    vvf_.clear();
    vSpheretf_.clear();
    vSymmtf_.clear();
    vtf_.clear();
    fieldMap_.clear();
    reverseFieldMap_.clear();


    // Generate fields with mappedField boundary condition

    // Convert field to map
    fieldMap_.reserve(fieldSet_.size());
    reverseFieldMap_.reserve(fieldSet_.size());
    forAll(fieldSet_, seti)
    {
        const word& fldName = fieldSet_[seti].first();
        const word& sampleFldName = fieldSet_[seti].second();

        fieldMap_.insert(fldName, sampleFldName);
        reverseFieldMap_.insert(sampleFldName, fldName);
    }

    Info<< type() << " " << name()
        << ": Sampling " << fieldMap_.size() << " fields" << endl;

    // Do analysis
    calcAddressing();

    return true;
}


bool Foam::functionObjects::nearWallFields::execute()
{
    DebugInFunction << endl;

    if
    (
        fieldMap_.size()
     && vsf_.empty()
     && vvf_.empty()
     && vSpheretf_.empty()
     && vSymmtf_.empty()
     && vtf_.empty()
    )
    {
        Log << type() << " " << name()
            << ": Creating " << fieldMap_.size() << " fields" << endl;

        createFields(vsf_);
        createFields(vvf_);
        createFields(vSpheretf_);
        createFields(vSymmtf_);
        createFields(vtf_);

        Log << endl;
    }

    Log << type() << " " << name()
        << " write:" << nl
        << "    Sampling fields to " << time_.timeName()
        << endl;

    sampleFields(vsf_);
    sampleFields(vvf_);
    sampleFields(vSpheretf_);
    sampleFields(vSymmtf_);
    sampleFields(vtf_);

    return true;
}


bool Foam::functionObjects::nearWallFields::write()
{
    DebugInFunction << endl;

    Log << "    Writing sampled fields to " << time_.timeName()
        << endl;

    forAll(vsf_, i)
    {
        vsf_[i].write();
    }
    forAll(vvf_, i)
    {
        vvf_[i].write();
    }
    forAll(vSpheretf_, i)
    {
        vSpheretf_[i].write();
    }
    forAll(vSymmtf_, i)
    {
        vSymmtf_[i].write();
    }
    forAll(vtf_, i)
    {
        vtf_[i].write();
    }

    return true;
}


// ************************************************************************* //
