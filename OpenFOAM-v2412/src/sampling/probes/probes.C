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

#include "probes.H"
#include "dictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "IOmanip.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        probes,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::probes::createProbeFiles(const wordList& fieldNames)
{
    // Open new output streams

    bool needsNewFiles = false;
    for (const word& fieldName : fieldNames)
    {
        if (!probeFilePtrs_.found(fieldName))
        {
            needsNewFiles = true;
            break;
        }
    }

    if (needsNewFiles && Pstream::master())
    {
        DebugInfo
            << "Probing fields: " << fieldNames << nl
            << "Probing locations: " << *this << nl
            << endl;

        // Put in undecomposed case
        // (Note: gives problems for distributed data running)

        fileName probeDir
        (
            mesh_.time().globalPath()
          / functionObject::outputPrefix
          / name()/mesh_.regionName()
            // Use startTime as the instance for output files
          / mesh_.time().timeName(mesh_.time().startTime().value())
        );
        probeDir.clean();  // Remove unneeded ".."

        // Create directory if needed
        Foam::mkDir(probeDir);

        for (const word& fieldName : fieldNames)
        {
            if (probeFilePtrs_.found(fieldName))
            {
                // Safety
                continue;
            }

            auto osPtr = autoPtr<OFstream>::New(probeDir/fieldName);
            auto& os = *osPtr;

            probeFilePtrs_.insert(fieldName, osPtr);

            DebugInfo<< "open probe stream: " << os.name() << endl;

            const unsigned int width(IOstream::defaultPrecision() + 7);
            os.setf(std::ios_base::left);

            forAll(*this, probei)
            {
                os  << "# Probe " << probei << ' ' << operator[](probei);

                if (processor_[probei] == -1)
                {
                    os  << "  # Not Found";
                }
                // Only for patchProbes
                else if (probei < patchIDList_.size())
                {
                    const label patchi = patchIDList_[probei];
                    if (patchi != -1)
                    {
                        const polyBoundaryMesh& bm = mesh_.boundaryMesh();
                        if
                        (
                            patchi < bm.nNonProcessor()
                         || processor_[probei] == Pstream::myProcNo()
                        )
                        {
                            os  << " at patch " << bm[patchi].name();
                        }
                        os  << " with a distance of "
                            << mag(operator[](probei)-oldPoints_[probei])
                            << " m to the original point "
                            << oldPoints_[probei];
                    }
                }

                os  << nl;
            }

            os  << setw(width) << "# Time";

            forAll(*this, probei)
            {
                if (includeOutOfBounds_ || processor_[probei] != -1)
                {
                    os  << ' ' << setw(width) << probei;
                }
            }
            os  << endl;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh)
{
    DebugInfo<< "probes: resetting sample locations" << endl;

    elementList_.resize_nocopy(pointField::size());
    faceList_.resize_nocopy(pointField::size());
    processor_.resize_nocopy(pointField::size());
    processor_ = -1;

    forAll(*this, probei)
    {
        const point& location = (*this)[probei];

        const label celli = mesh.findCell(location);

        elementList_[probei] = celli;

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];
            const vector& cellCentre = mesh.cellCentres()[celli];
            scalar minDistance = GREAT;
            label minFaceID = -1;
            forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector dist = mesh.faceCentres()[facei] - cellCentre;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = facei;
                }
            }
            faceList_[probei] = minFaceID;
        }
        else
        {
            faceList_[probei] = -1;
        }

        if (debug && (elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "probes : found point " << location
                << " in cell " << elementList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(elementList_, probei)
    {
        const point& location = operator[](probei);
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        processor_[probei] = (celli != -1 ? Pstream::myProcNo() : -1);

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());
        reduce(processor_[probei], maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[probei] != -1 && elementList_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::probes::prepare(unsigned request)
{
    // Prefilter on selection
    HashTable<wordHashSet> selected =
    (
        loadFromFiles_
      ? IOobjectList(mesh_, mesh_.time().timeName()).classes(fieldSelection_)
      : mesh_.classes(fieldSelection_)
    );

    // Classify and count fields
    label nFields = 0;
    do
    {
        #undef  doLocalCode
        #define doLocalCode(InputType, Target)                                \
        {                                                                     \
            Target.clear();  /* Remove old values */                          \
            const auto iter = selected.cfind(InputType::typeName);            \
            if (iter.good())                                                  \
            {                                                                 \
                /* Add new (current) values */                                \
                Target.append(iter.val().sortedToc());                        \
                nFields += Target.size();                                     \
            }                                                                 \
        }

        doLocalCode(volScalarField, scalarFields_);
        doLocalCode(volVectorField, vectorFields_)
        doLocalCode(volSphericalTensorField, sphericalTensorFields_);
        doLocalCode(volSymmTensorField, symmTensorFields_);
        doLocalCode(volTensorField, tensorFields_);

        doLocalCode(surfaceScalarField, surfaceScalarFields_);
        doLocalCode(surfaceVectorField, surfaceVectorFields_);
        doLocalCode(surfaceSphericalTensorField, surfaceSphericalTensorFields_);
        doLocalCode(surfaceSymmTensorField, surfaceSymmTensorFields_);
        doLocalCode(surfaceTensorField, surfaceTensorFields_);
        #undef doLocalCode
    }
    while (false);


    // Adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields(2*nFields);
        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        DebugInfo
            << "Probing fields: " << currentFields << nl
            << "Probing locations: " << *this << nl
            << endl;

        // Close streams for fields that no longer exist
        forAllIters(probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                DebugInfo<< "close probe stream: " << iter()->name() << endl;

                probeFilePtrs_.remove(iter);
            }
        }

        if ((request & ACTION_WRITE) && !currentFields.empty())
        {
            createProbeFiles(currentFields.sortedToc());
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    pointField(),
    loadFromFiles_(loadFromFiles),
    fixedLocations_(true),
    includeOutOfBounds_(true),
    verbose_(false),
    onExecute_(false),
    fieldSelection_(),
    samplePointScheme_("cell")
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::probes::verbose(const bool on) noexcept
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


bool Foam::probes::read(const dictionary& dict)
{
    dict.readEntry("probeLocations", static_cast<pointField&>(*this));
    dict.readEntry("fields", fieldSelection_);

    dict.readIfPresent("fixedLocations", fixedLocations_);
    dict.readIfPresent("includeOutOfBounds", includeOutOfBounds_);

    verbose_ = dict.getOrDefault("verbose", false);
    onExecute_ = dict.getOrDefault("sampleOnExecute", false);

    if (dict.readIfPresent("interpolationScheme", samplePointScheme_))
    {
        if (!fixedLocations_ && samplePointScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations. InterpolationScheme "
                << "entry will be ignored"
                << endl;
        }
    }

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    // Close old (ununsed) streams
    prepare(ACTION_NONE);

    return true;
}


bool Foam::probes::performAction(unsigned request)
{
    if (!pointField::empty() && request && prepare(request))
    {
        performAction(scalarFields_, request);
        performAction(vectorFields_, request);
        performAction(sphericalTensorFields_, request);
        performAction(symmTensorFields_, request);
        performAction(tensorFields_, request);

        performAction(surfaceScalarFields_, request);
        performAction(surfaceVectorFields_, request);
        performAction(surfaceSphericalTensorFields_, request);
        performAction(surfaceSymmTensorFields_, request);
        performAction(surfaceTensorFields_, request);
    }
    return true;
}


bool Foam::probes::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


bool Foam::probes::write()
{
    return performAction(ACTION_ALL);
}


void Foam::probes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "probes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        findElements(mesh_);
    }
    else
    {
        DebugInfo<< "probes: remapping sample locations" << endl;

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                if (celli != -1)
                {
                    label newCelli = reverseMap[celli];
                    if (newCelli == -1)
                    {
                        // cell removed
                    }
                    else if (newCelli < -1)
                    {
                        // cell merged
                        elems.append(-newCelli - 2);
                    }
                    else
                    {
                        // valid new cell
                        elems.append(newCelli);
                    }
                }
                else
                {
                    // Keep -1 elements so the size stays the same
                    elems.append(-1);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            for (const label facei : faceList_)
            {
                if (facei != -1)
                {
                    label newFacei = reverseMap[facei];
                    if (newFacei == -1)
                    {
                        // face removed
                    }
                    else if (newFacei < -1)
                    {
                        // face merged
                        elems.append(-newFacei - 2);
                    }
                    else
                    {
                        // valid new face
                        elems.append(newFacei);
                    }
                }
                else
                {
                    // Keep -1 elements
                    elems.append(-1);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::probes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        findElements(mesh_);
    }
}


// ************************************************************************* //
