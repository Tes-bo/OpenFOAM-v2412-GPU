/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "polyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "SubField.H"
#include "entry.H"
#include "dictionary.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyPatch, 0);

    int polyPatch::disallowGenericPolyPatch
    (
        debug::debugSwitch("disallowGenericPolyPatch", 0)
    );

    defineRunTimeSelectionTable(polyPatch, word);
    defineRunTimeSelectionTable(polyPatch, dictionary);

    addToRunTimeSelectionTable(polyPatch, polyPatch, word);
    addToRunTimeSelectionTable(polyPatch, polyPatch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::polyPatch::movePoints(PstreamBuffers&, const pointField& p)
{
    primitivePatch::movePoints(p);
}


void Foam::polyPatch::updateMesh(PstreamBuffers&)
{
    primitivePatch::clearGeom();
    clearAddressing();
}


void Foam::polyPatch::clearGeom()
{
    primitivePatch::clearGeom();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, index),
    primitivePatch
    (
        faceSubList(bm.mesh().faces(), size, start),
        bm.mesh().points()
    ),
    start_(start),
    boundaryMesh_(bm)
{
    if (constraintType(patchType))
    {
        addGroup(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& physicalType,
    const wordList& inGroups
)
:
    patchIdentifier(name, index, physicalType, inGroups),
    primitivePatch
    (
        faceSubList(bm.mesh().faces(), size, start),
        bm.mesh().points()
    ),
    start_(start),
    boundaryMesh_(bm)
{}


Foam::polyPatch::polyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, dict, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            dict.get<label>("nFaces"),
            dict.get<label>("startFace")
        ),
        bm.mesh().points()
    ),
    start_(dict.get<label>("startFace")),
    boundaryMesh_(bm)
{
    if (constraintType(patchType))
    {
        addGroup(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(pp),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            pp.size(),
            pp.start()
        ),
        bm.mesh().points()
    ),
    start_(pp.start()),
    boundaryMesh_(bm)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    patchIdentifier(pp, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            newSize,
            newStart
        ),
        bm.mesh().points()
    ),
    start_(newStart),
    boundaryMesh_(bm)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    patchIdentifier(pp, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            mapAddressing.size(),
            newStart
        ),
        bm.mesh().points()
    ),
    start_(newStart),
    boundaryMesh_(bm)
{}


Foam::polyPatch::polyPatch(const polyPatch& p)
:
    patchIdentifier(p),
    primitivePatch(p),
    start_(p.start_),
    boundaryMesh_(p.boundaryMesh_)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& p,
    const labelList& faceCells
)
:
    polyPatch(p)
{
    faceCellsPtr_.reset(new labelList::subList(faceCells, faceCells.size()));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyPatch::~polyPatch()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyPatch::constraintType(const word& patchType)
{
    return
    (
        !patchType.empty()
     && pointPatchField<scalar>::patchConstructorTablePtr_
     && pointPatchField<scalar>::patchConstructorTablePtr_->found(patchType)
    );
}


Foam::wordList Foam::polyPatch::constraintTypes()
{
    const auto& cnstrTable = *dictionaryConstructorTablePtr_;

    wordList cTypes(cnstrTable.size());

    label i = 0;

    forAllConstIters(cnstrTable, iter)
    {
        if (constraintType(iter.key()))
        {
            cTypes[i++] = iter.key();
        }
    }

    cTypes.resize(i);

    return cTypes;
}


Foam::label Foam::polyPatch::offset() const noexcept
{
    // Same as start_ - polyMesh::nInternalFaces()
    return start_ - boundaryMesh_.start();
}


const Foam::polyBoundaryMesh& Foam::polyPatch::boundaryMesh() const noexcept
{
    return boundaryMesh_;
}


const Foam::vectorField::subField Foam::polyPatch::faceCentres() const
{
    return patchSlice(boundaryMesh().mesh().faceCentres());
}


const Foam::vectorField::subField Foam::polyPatch::faceAreas() const
{
    return patchSlice(boundaryMesh().mesh().faceAreas());
}


Foam::tmp<Foam::vectorField> Foam::polyPatch::faceCellCentres() const
{
    auto tcc = tmp<vectorField>::New(size());
    auto& cc = tcc.ref();

    // get reference to global cell centres
    const vectorField& gcc = boundaryMesh_.mesh().cellCentres();

    const labelUList& faceCells = this->faceCells();

    forAll(faceCells, facei)
    {
        cc[facei] = gcc[faceCells[facei]];
    }

    return tcc;
}


Foam::tmp<Foam::scalarField> Foam::polyPatch::areaFraction
(
    const pointField& points
) const
{
    auto tfraction = tmp<scalarField>::New(size());
    auto& fraction = tfraction.ref();

    const vectorField::subField faceAreas = this->faceAreas();

    forAll(*this, facei)
    {
        const face& f = this->operator[](facei);
        fraction[facei] = faceAreas[facei].mag()/(f.mag(points) + ROOTVSMALL);
    }

    return tfraction;
}


Foam::tmp<Foam::scalarField> Foam::polyPatch::areaFraction() const
{
    if (areaFractionPtr_)
    {
        return tmp<scalarField>(*areaFractionPtr_);
    }
    return nullptr;
}


void Foam::polyPatch::areaFraction(const scalar fraction)
{
    areaFractionPtr_ = std::make_unique<scalarField>(size(), fraction);
}


void Foam::polyPatch::areaFraction(const tmp<scalarField>& fraction)
{
    if (fraction)
    {
        // Steal or clone
        areaFractionPtr_.reset(fraction.ptr());
    }
    else
    {
        areaFractionPtr_.reset(nullptr);
    }
}


const Foam::labelUList& Foam::polyPatch::faceCells() const
{
    if (!faceCellsPtr_)
    {
        faceCellsPtr_.reset
        (
            new labelList::subList
            (
                patchSlice(boundaryMesh().mesh().faceOwner())
            )
        );
    }

    return *faceCellsPtr_;
}


const Foam::labelList& Foam::polyPatch::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_.reset
        (
            new labelList
            (
                primitivePatch::meshEdges
                (
                    boundaryMesh().mesh().edges(),
                    boundaryMesh().mesh().pointEdges()
                )
            )
        );
    }

    return *mePtr_;
}


void Foam::polyPatch::clearAddressing()
{
    primitivePatch::clearTopology();
    primitivePatch::clearPatchMeshAddr();
    faceCellsPtr_.reset(nullptr);
    mePtr_.reset(nullptr);
    areaFractionPtr_.reset(nullptr);
}


void Foam::polyPatch::write(Ostream& os) const
{
    os.writeEntry("type", type());
    patchIdentifier::write(os);
    os.writeEntry("nFaces", size());
    os.writeEntry("startFace", start());
}


void Foam::polyPatch::initOrder(PstreamBuffers&, const primitivePatch&) const
{}


bool Foam::polyPatch::order
(
    PstreamBuffers&,
    const primitivePatch&,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Nothing changed.
    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::polyPatch::operator=(const polyPatch& p)
{
    clearAddressing();

    patchIdentifier::operator=(p);
    primitivePatch::operator=(p);
    start_ = p.start_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyPatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
