/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "moleculeCloud.H"
#include "fvMesh.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<molecule>, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::moleculeCloud::buildConstProps()
{
    Info<< nl << "Reading moleculeProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size());

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    forAll(idList, i)
    {
        const word& id = idList[i];
        const dictionary& molDict = moleculePropertiesDict.subDict(id);

        List<word> siteIdNames
        (
            molDict.lookup("siteIds")
        );

        List<label> siteIds(siteIdNames.size());

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            siteIds[sI] = siteIdList.find(siteId);

            if (siteIds[sI] == -1)
            {
                FatalErrorInFunction
                    << siteId << " site not found."
                    << nl << abort(FatalError);
            }
        }

        molecule::constantProperties& constProp = constPropList_[i];

        constProp = molecule::constantProperties(molDict);

        constProp.siteIds() = siteIds;
    }
}


void Foam::moleculeCloud::setSiteSizesAndPositions()
{
    for (molecule& mol : *this)
    {
        const molecule::constantProperties& cP = constProps(mol.id());

        mol.setSiteSizes(cP.nSites());

        mol.setSitePositions(cP);
    }
}


void Foam::moleculeCloud::buildCellOccupancy()
{
    for (auto& list : cellOccupancy_)
    {
        list.clear();
    }

    for (molecule& mol : *this)
    {
        cellOccupancy_[mol.cell()].append(&mol);
    }

    for (auto& list : cellOccupancy_)
    {
        list.shrink();
    }
}


void Foam::moleculeCloud::calculatePairForce()
{
    PstreamBuffers pBufs;

    // Start sending referred data
    label startOfRequests = Pstream::nRequests();
    il_.sendReferredData(cellOccupancy(), pBufs);

    molecule* molI = nullptr;
    molecule* molJ = nullptr;

    {
        // Real-Real interactions

        const labelListList& dil = il_.dil();

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<molecule*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        evaluatePair(*molI, *molJ);
                    }
                }

                forAll(cellOccupancy_[d], cellIOtherMols)
                {
                    molJ = cellOccupancy_[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        evaluatePair(*molI, *molJ);
                    }
                }
            }
        }
    }

    // Receive referred data
    il_.receiveReferredData(pBufs, startOfRequests);

    {
        // Real-Referred interactions

        const labelListList& ril = il_.ril();

        List<IDLList<molecule>>& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            const List<label>& realCells = ril[r];

            IDLList<molecule>& refMols = referredMols[r];

            for (molecule& refMol : refMols)
            {
                forAll(realCells, rC)
                {
                    List<molecule*> celli = cellOccupancy_[realCells[rC]];

                    forAll(celli, cellIMols)
                    {
                        molI = celli[cellIMols];

                        evaluatePair(*molI, refMol);
                    }
                }
            }
        }
    }
}


void Foam::moleculeCloud::calculateTetherForce()
{
    const tetherPotentialList& tetherPot(pot_.tetherPotentials());

    for (molecule& mol : *this)
    {
        if (mol.tethered())
        {
            vector rIT = mol.position() - mol.specialPosition();

            label idI = mol.id();

            scalar massI = constProps(idI).mass();

            vector fIT = tetherPot.force(idI, rIT);

            mol.a() += fIT/massI;

            mol.potentialEnergy() += tetherPot.energy(idI, rIT);

            // What to do here?
            // mol.rf() += rIT*fIT;
        }
    }
}


void Foam::moleculeCloud::calculateExternalForce()
{
    for (molecule& mol : *this)
    {
        mol.a() += pot_.gravity();
    }
}


void Foam::moleculeCloud::removeHighEnergyOverlaps()
{
    Info<< nl << "Removing high energy overlaps, limit = "
        << pot_.potentialEnergyLimit()
        << nl << "Removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        Info<< ' ' << pot_.idList()[pot_.removalOrder()[rO]];
    }

    Info<< nl ;

    label initialSize = this->size();

    buildCellOccupancy();

    // Real-Real interaction

    molecule* molI = nullptr;
    molecule* molJ = nullptr;

    {
        DynamicList<molecule*> molsToDelete;

        const labelListList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<molecule*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        if (evaluatePotentialLimit(*molI, *molJ))
                        {
                            const label idI = molI->id();
                            const label idJ = molJ->id();

                            if
                            (
                                idI == idJ
                             || pot_.removalOrder().find(idJ)
                              < pot_.removalOrder().find(idI)
                            )
                            {
                                if (!molsToDelete.found(molJ))
                                {
                                    molsToDelete.append(molJ);
                                }
                            }
                            else if (!molsToDelete.found(molI))
                            {
                                molsToDelete.append(molI);
                            }
                        }
                    }
                }
            }

            forAll(cellOccupancy_[d], cellIOtherMols)
            {
                molJ = cellOccupancy_[d][cellIOtherMols];

                if (molJ > molI)
                {
                    if (evaluatePotentialLimit(*molI, *molJ))
                    {
                        const label idI = molI->id();
                        const label idJ = molJ->id();

                        if
                        (
                            idI == idJ
                         || pot_.removalOrder().find(idJ)
                          < pot_.removalOrder().find(idI)
                        )
                        {
                            if (!molsToDelete.found(molJ))
                            {
                                molsToDelete.append(molJ);
                            }
                        }
                        else if (!molsToDelete.found(molI))
                        {
                            molsToDelete.append(molI);
                        }
                    }
                }
            }
        }

        forAll(molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    PstreamBuffers pBufs;

    // Start sending referred data
    label startOfRequests = Pstream::nRequests();

    il_.sendReferredData(cellOccupancy(), pBufs);

        // Receive referred data
    il_.receiveReferredData(pBufs, startOfRequests);

    // Real-Referred interaction

    {
        DynamicList<molecule*> molsToDelete;

        const labelListList& ril(il_.ril());

        List<IDLList<molecule>>& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            IDLList<molecule>& refMols = referredMols[r];

            for (molecule& refMol : refMols)
            {
                molJ = &refMol;

                const List<label>& realCells = ril[r];

                forAll(realCells, rC)
                {
                    label celli = realCells[rC];

                    List<molecule*> cellIMols = cellOccupancy_[celli];

                    forAll(cellIMols, cIM)
                    {
                        molI = cellIMols[cIM];

                        if (evaluatePotentialLimit(*molI, *molJ))
                        {
                            const label idI = molI->id();
                            const label idJ = molJ->id();

                            if
                            (
                                pot_.removalOrder().find(idI)
                              < pot_.removalOrder().find(idJ)
                            )
                            {
                                if (!molsToDelete.found(molI))
                                {
                                    molsToDelete.append(molI);
                                }
                            }
                            else if
                            (
                                pot_.removalOrder().find(idI)
                             == pot_.removalOrder().find(idJ)
                            )
                            {
                                // Remove one of the molecules
                                // arbitrarily, assuring that a
                                // consistent decision is made for
                                // both real-referred pairs.

                                if (molI->origId() > molJ->origId())
                                {
                                    if (!molsToDelete.found(molI))
                                    {
                                        molsToDelete.append(molI);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        forAll(molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    // Start sending referred data
    startOfRequests = Pstream::nRequests();

    il_.sendReferredData(cellOccupancy(), pBufs);

    // Receive referred data
    il_.receiveReferredData(pBufs, startOfRequests);

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab << molsRemoved << " molecules removed" << endl;
}


void Foam::moleculeCloud::initialiseMolecules
(
    const IOdictionary& mdInitialiseDict
)
{
    Info<< nl
        << "Initialising molecules in each zone specified in mdInitialiseDict."
        << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (!cellZones.size())
    {
        FatalErrorInFunction
            << "No cellZones found in the mesh."
            << abort(FatalError);
    }

    for (const cellZone& zone : cellZones)
    {
        if (zone.size())
        {
            if (!mdInitialiseDict.found(zone.name()))
            {
                Info<< "No specification subDictionary for zone "
                    << zone.name() << " found, skipping." << endl;
            }
            else
            {
                const dictionary& zoneDict =
                    mdInitialiseDict.subDict(zone.name());

                const scalar temperature
                (
                    zoneDict.get<scalar>("temperature")
                );

                const vector bulkVelocity
                (
                    zoneDict.get<vector>("bulkVelocity")
                );

                List<word> latticeIds
                (
                    zoneDict.lookup("latticeIds")
                );

                List<vector> latticePositions
                (
                    zoneDict.lookup("latticePositions")
                );

                if (latticeIds.size() != latticePositions.size())
                {
                    FatalErrorInFunction
                        << "latticeIds and latticePositions must be the same "
                        << " size." << nl
                        << abort(FatalError);
                }

                diagTensor latticeCellShape
                (
                    zoneDict.lookup("latticeCellShape")
                );

                scalar latticeCellScale = 0.0;

                if (zoneDict.found("numberDensity"))
                {
                    const scalar numberDensity
                    (
                        zoneDict.get<scalar>("numberDensity")
                    );

                    if (numberDensity < VSMALL)
                    {
                        WarningInFunction
                            << "numberDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }

                    latticeCellScale = pow
                    (
                        latticeIds.size()/(det(latticeCellShape)*numberDensity),
                        (1.0/3.0)
                    );
                }
                else if (zoneDict.found("massDensity"))
                {
                    scalar unitCellMass = 0.0;

                    forAll(latticeIds, i)
                    {
                        label id = pot_.idList().find(latticeIds[i]);

                        const molecule::constantProperties& cP(constProps(id));

                        unitCellMass += cP.mass();
                    }

                    const scalar massDensity
                    (
                        zoneDict.get<scalar>("massDensity")
                    );

                    if (massDensity < VSMALL)
                    {
                        WarningInFunction
                            << "massDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }


                    latticeCellScale = pow
                    (
                        unitCellMass/(det(latticeCellShape)*massDensity),
                        (1.0/3.0)
                    );
                }
                else
                {
                    FatalErrorInFunction
                        << "massDensity or numberDensity not specified " << nl
                        << abort(FatalError);
                }

                latticeCellShape *= latticeCellScale;

                vector anchor(zoneDict.get<vector>("anchor"));

                bool tethered = false;

                if (zoneDict.found("tetherSiteIds"))
                {
                    tethered = bool
                    (
                        List<word>(zoneDict.lookup("tetherSiteIds")).size()
                    );
                }

                const vector orientationAngles
                (
                    zoneDict.get<vector>("orientationAngles")
                );

                const scalar phi(degToRad(orientationAngles.x()));
                const scalar theta(degToRad(orientationAngles.y()));
                const scalar psi(degToRad(orientationAngles.z()));

                const tensor R
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                );

                // Find the optimal anchor position.  Finding the approximate
                // mid-point of the zone of cells and snapping to the nearest
                // lattice location.

                vector zoneMin = VGREAT*vector::one;

                vector zoneMax = -VGREAT*vector::one;

                forAll(zone, cell)
                {
                    const point cellCentre = mesh_.cellCentres()[zone[cell]];

                    if (cellCentre.x() > zoneMax.x())
                    {
                        zoneMax.x() = cellCentre.x();
                    }
                    if (cellCentre.x() < zoneMin.x())
                    {
                        zoneMin.x() = cellCentre.x();
                    }
                    if (cellCentre.y() > zoneMax.y())
                    {
                        zoneMax.y() = cellCentre.y();
                    }
                    if (cellCentre.y() < zoneMin.y())
                    {
                        zoneMin.y() = cellCentre.y();
                    }
                    if (cellCentre.z() > zoneMax.z())
                    {
                        zoneMax.z() = cellCentre.z();
                    }
                    if (cellCentre.z() < zoneMin.z())
                    {
                        zoneMin.z() = cellCentre.z();
                    }
                }

                point zoneMid = 0.5*(zoneMin + zoneMax);

                point latticeMid = inv(latticeCellShape) & (R.T()
                & (zoneMid - anchor));

                point latticeAnchor
                (
                    label(latticeMid.x() + 0.5*sign(latticeMid.x())),
                    label(latticeMid.y() + 0.5*sign(latticeMid.y())),
                    label(latticeMid.z() + 0.5*sign(latticeMid.z()))
                );

                anchor += (R & (latticeCellShape & latticeAnchor));

                // Continue trying to place molecules as long as at
                // least one molecule is placed in each iteration.
                // The "|| totalZoneMols == 0" condition means that the
                // algorithm will continue if the origin is outside the
                // zone.

                label n = 0;

                label totalZoneMols = 0;

                label molsPlacedThisIteration = 0;

                while
                (
                    molsPlacedThisIteration != 0
                 || totalZoneMols == 0
                )
                {
                    label sizeBeforeIteration = this->size();

                    bool partOfLayerInBounds = false;

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // start of placement of molecules
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if (n == 0)
                    {
                        // Special treatment is required for the first position,
                        // i.e. iteration zero.

                        labelVector unitCellLatticePosition(0,0,0);

                        forAll(latticePositions, p)
                        {
                            label id = pot_.idList().find(latticeIds[p]);

                            const vector& latticePosition =
                                vector
                                (
                                    unitCellLatticePosition.x(),
                                    unitCellLatticePosition.y(),
                                    unitCellLatticePosition.z()
                                )
                              + latticePositions[p];

                            point globalPosition =
                                anchor
                              + (R & (latticeCellShape & latticePosition));

                            partOfLayerInBounds = mesh_.bounds().contains
                            (
                                globalPosition
                            );

                            const label cell =
                                mesh_.cellTree().findInside(globalPosition);

                            if (zone.found(cell))
                            {
                                createMolecule
                                (
                                    globalPosition,
                                    cell,
                                    id,
                                    tethered,
                                    temperature,
                                    bulkVelocity
                                );
                            }
                        }
                    }
                    else
                    {
                        // Place top and bottom caps.

                        labelVector unitCellLatticePosition(0,0,0);

                        for
                        (
                            unitCellLatticePosition.z() = -n;
                            unitCellLatticePosition.z() <= n;
                            unitCellLatticePosition.z() += 2*n
                        )
                        {
                            for
                            (
                                unitCellLatticePosition.y() = -n;
                                unitCellLatticePosition.y() <= n;
                                unitCellLatticePosition.y()++
                            )
                            {
                                for
                                (
                                    unitCellLatticePosition.x() = -n;
                                    unitCellLatticePosition.x() <= n;
                                    unitCellLatticePosition.x()++
                                )
                                {
                                    forAll(latticePositions, p)
                                    {
                                        const label id =
                                            pot_.idList().find(latticeIds[p]);

                                        const vector& latticePosition =
                                            vector
                                            (
                                                unitCellLatticePosition.x(),
                                                unitCellLatticePosition.y(),
                                                unitCellLatticePosition.z()
                                            )
                                          + latticePositions[p];

                                        point globalPosition =
                                            anchor
                                          + (
                                                R
                                              & (
                                                    latticeCellShape
                                                  & latticePosition
                                                )
                                            );

                                        partOfLayerInBounds =
                                            mesh_.bounds().contains
                                            (
                                                globalPosition
                                            );

                                        const label cell =
                                            mesh_.cellTree().findInside
                                            (
                                                globalPosition
                                            );

                                        if (zone.found(cell))
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }
                                }
                            }
                        }

                        for
                        (
                            unitCellLatticePosition.z() = -(n-1);
                            unitCellLatticePosition.z() <= (n-1);
                            unitCellLatticePosition.z()++
                        )
                        {
                            for (label iR = 0; iR <= 2*n -1; iR++)
                            {
                                unitCellLatticePosition.x() = n;

                                unitCellLatticePosition.y() = -n + (iR + 1);

                                for (label iK = 0; iK < 4; iK++)
                                {
                                    forAll(latticePositions, p)
                                    {
                                        const label id =
                                            pot_.idList().find(latticeIds[p]);

                                        const vector& latticePosition =
                                            vector
                                            (
                                                unitCellLatticePosition.x(),
                                                unitCellLatticePosition.y(),
                                                unitCellLatticePosition.z()
                                            )
                                          + latticePositions[p];

                                        point globalPosition =
                                            anchor
                                          + (
                                                R
                                              & (
                                                   latticeCellShape
                                                 & latticePosition
                                                )
                                            );

                                        partOfLayerInBounds =
                                            mesh_.bounds().contains
                                            (
                                                globalPosition
                                            );

                                        const label cell =
                                            mesh_.cellTree().findInside
                                            (
                                                globalPosition
                                            );

                                        if (zone.found(cell))
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }

                                    unitCellLatticePosition =
                                        labelVector
                                        (
                                          - unitCellLatticePosition.y(),
                                            unitCellLatticePosition.x(),
                                            unitCellLatticePosition.z()
                                        );
                                }
                            }
                        }
                    }

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // end of placement of molecules
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if
                    (
                        totalZoneMols == 0
                     && !partOfLayerInBounds
                    )
                    {
                        WarningInFunction
                            << "A whole layer of unit cells was placed "
                            << "outside the bounds of the mesh, but no "
                            << "molecules have been placed in zone '"
                            << zone.name()
                            << "'.  This is likely to be because the zone "
                            << "has few cells ("
                            << zone.size()
                            << " in this case) and no lattice position "
                            << "fell inside them.  "
                            << "Aborting filling this zone."
                            << endl;

                        break;
                    }

                    molsPlacedThisIteration =
                        this->size() - sizeBeforeIteration;

                    totalZoneMols += molsPlacedThisIteration;

                    n++;
                }
            }
        }
    }
}


void Foam::moleculeCloud::createMolecule
(
    const point& position,
    label cell,
    label id,
    bool tethered,
    scalar temperature,
    const vector& bulkVelocity
)
{
    point specialPosition(Zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = molecule::SPECIAL_TETHERED;
    }

    const molecule::constantProperties& cP(constProps(id));

    vector v = equipartitionLinearVelocity(temperature, cP.mass());

    v += bulkVelocity;

    vector pi = Zero;

    tensor Q = I;

    if (!cP.pointMolecule())
    {
        pi = equipartitionAngularMomentum(temperature, cP);

        const scalar phi(rndGen_.sample01<scalar>()*mathematical::twoPi);
        const scalar theta(rndGen_.sample01<scalar>()*mathematical::twoPi);
        const scalar psi(rndGen_.sample01<scalar>()*mathematical::twoPi);

        Q = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );
    }

    addParticle
    (
        new molecule
        (
            mesh_,
            position,
            cell,
            Q,
            v,
            Zero,
            pi,
            Zero,
            specialPosition,
            constProps(id),
            special,
            id
        )
    );
}


Foam::label Foam::moleculeCloud::nSites() const
{
    label n = 0;

    for (const molecule& mol : *this)
    {
        n += constProps(mol.id()).nSites();
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    const potential& pot,
    bool readFields
)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    cellOccupancy_(mesh_.nCells()),
    il_(mesh_, pot_.pairPotentials().rCutMax(), false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        molecule::readFields(*this);
    }

    buildConstProps();

    setSiteSizesAndPositions();

    removeHighEnergyOverlaps();

    calculateForce();
}


Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    const potential& pot,
    const IOdictionary& mdInitialiseDict,
    bool readFields
)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    il_(mesh_, 0.0, false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        molecule::readFields(*this);
    }

    clear();

    buildConstProps();

    initialiseMolecules(mdInitialiseDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeCloud::evolve()
{
    molecule::trackingData td0(*this, 0);
    Cloud<molecule>::move(*this, td0, mesh_.time().deltaTValue());

    molecule::trackingData td1(*this, 1);
    Cloud<molecule>::move(*this, td1, mesh_.time().deltaTValue());

    molecule::trackingData td2(*this, 2);
    Cloud<molecule>::move(*this, td2, mesh_.time().deltaTValue());

    calculateForce();

    molecule::trackingData td3(*this, 3);
    Cloud<molecule>::move(*this, td3, mesh_.time().deltaTValue());
}


void Foam::moleculeCloud::calculateForce()
{
    buildCellOccupancy();

    // Set accumulated quantities to zero
    for (molecule& mol : *this)
    {
        mol.siteForces() = Zero;

        mol.potentialEnergy() = 0.0;

        mol.rf() = Zero;
    }

    calculatePairForce();

    calculateTetherForce();

    calculateExternalForce();
}


void Foam::moleculeCloud::applyConstraintsAndThermostats
(
    const scalar targetTemperature,
    const scalar measuredTemperature
)
{
    scalar temperatureCorrectionFactor =
        sqrt(targetTemperature/max(VSMALL, measuredTemperature));

    Info<< "----------------------------------------" << nl
        << "Temperature equilibration" << nl
        << "Target temperature = "
        << targetTemperature << nl
        << "Measured temperature = "
        << measuredTemperature << nl
        << "Temperature correction factor = "
        << temperatureCorrectionFactor << nl
        << "----------------------------------------"
        << endl;

    for (molecule& mol : *this)
    {
        mol.v() *= temperatureCorrectionFactor;

        mol.pi() *= temperatureCorrectionFactor;
    }
}


void Foam::moleculeCloud::writeXYZ(const fileName& fName) const
{
    OFstream os(fName);

    os  << nSites() << nl << "moleculeCloud site positions in angstroms" << nl;

    for (const molecule& mol : *this)
    {
        const molecule::constantProperties& cP = constProps(mol.id());

        forAll(mol.sitePositions(), i)
        {
            const point& sP = mol.sitePositions()[i];

            os  << pot_.siteIdList()[cP.siteIds()[i]]
                << ' ' << sP.x()*1e10
                << ' ' << sP.y()*1e10
                << ' ' << sP.z()*1e10
                << nl;
        }
    }
}


// ************************************************************************* //
