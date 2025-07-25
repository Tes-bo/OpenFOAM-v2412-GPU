/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "singleProcessorFaceSetsConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(singleProcessorFaceSets);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        singleProcessorFaceSets,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionConstraints::singleProcessorFaceSets::printInfo() const
{
    for (const auto& nameAndProc : setNameAndProcs_)
    {
        Info<< "    all cells connected to faceSet "
            << nameAndProc.first()
            << " on processor " << nameAndProc.second() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::singleProcessorFaceSets::
singleProcessorFaceSets
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    setNameAndProcs_
    (
        coeffDict_.lookupCompat("sets", {{"singleProcessorFaceSets", 1806}})
    )
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        printInfo();
    }
}


Foam::decompositionConstraints::singleProcessorFaceSets::
singleProcessorFaceSets
(
    const List<Tuple2<word, label>>& setNameAndProcs
)
:
    decompositionConstraint(dictionary(), typeName),
    setNameAndProcs_(setNameAndProcs)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        printInfo();
    }
}


Foam::decompositionConstraints::singleProcessorFaceSets::
singleProcessorFaceSets
(
    Istream& is
)
:
    decompositionConstraint(dictionary(), typeName),
    setNameAndProcs_(is)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        printInfo();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::singleProcessorFaceSets::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.resize(mesh.nFaces(), true);

    // Mark faces already in set
    labelList faceToSet(mesh.nFaces(), -1);
    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& faceLabels = specifiedProcessorFaces[setI];
        for (const label facei : faceLabels)
        {
            if (faceToSet[facei] == -1)
            {
                faceToSet[facei] = setI;
            }
            else if (faceToSet[facei] != setI)
            {
                WarningInFunction << "Face " << facei
                    << " at " << mesh.faceCentres()[facei]
                    << " is already in existing constraint "
                    << faceToSet[facei]
                    << endl;
            }
        }
    }

    forAll(setNameAndProcs_, setI)
    {
        //Info<< "Keeping all cells connected to faceSet "
        //    << setNameAndProcs_[setI].first()
        //    << " on processor " << setNameAndProcs_[setI].second() << endl;

        const label destProcI = setNameAndProcs_[setI].second();

        // Read faceSet
        const faceSet fz(mesh, setNameAndProcs_[setI].first());

        // Check that it does not overlap with existing specifiedProcessorFaces
        labelList nMatch(specifiedProcessorFaces.size(), Zero);
        for (const label facei : fz)
        {
            const label seti = faceToSet[facei];
            if (seti != -1)
            {
                ++nMatch[seti];
            }
        }


        // Only store if all faces are not yet in specifiedProcessorFaces
        // (on all processors)
        bool store = true;

        forAll(nMatch, setI)
        {
            if (nMatch[setI] == fz.size())
            {
                // Full match
                store = false;
                break;
            }
            else if (nMatch[setI] > 0)
            {
                // Partial match
                store = false;
                break;
            }
        }

        if (returnReduceOr(store))
        {
            specifiedProcessorFaces.append(new labelList(fz.sortedToc()));
            specifiedProcessor.append(destProcI);
        }
    }


    // Unblock all point connected faces
    // 1. Mark all points on specifiedProcessorFaces
    boolList procFacePoint(mesh.nPoints(), false);
    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& faceLabels = specifiedProcessorFaces[setI];
        for (const label facei : faceLabels)
        {
            const face& f = mesh.faces()[facei];

            for (const label pointi : f)
            {
                procFacePoint[pointi] = true;
            }
        }
    }
    syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

    // 2. Unblock all faces on procFacePoint

    label nUnblocked = 0;

    forAll(procFacePoint, pointi)
    {
        if (procFacePoint[pointi])
        {
            const labelList& pFaces = mesh.pointFaces()[pointi];
            forAll(pFaces, i)
            {
                if (blockedFace[pFaces[i]])
                {
                    blockedFace[pFaces[i]] = false;
                    ++nUnblocked;
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : unblocked "
            << returnReduce(nUnblocked, sumOp<label>()) << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::singleProcessorFaceSets::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // For specifiedProcessorFaces rework the cellToProc to enforce
    // all on one processor since we can't guarantee that the input
    // to regionSplit was a single region.
    // E.g. faceSet 'a' with the cells split into two regions
    // by a notch formed by two walls
    //
    //          \   /
    //           \ /
    //    ---a----+-----a-----
    //
    //
    // Note that reworking the cellToProc might make the decomposition
    // unbalanced.
    label nChanged = 0;

    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& set = specifiedProcessorFaces[setI];

        // Get the processor to use for the set
        label procI = specifiedProcessor[setI];
        if (procI == -1)
        {
            // If no processor specified use the one from the
            // 0th element
            if (set.size())
            {
                procI = decomposition[mesh.faceOwner()[set[0]]];
            }
            reduce(procI, maxOp<label>());
        }

        // Get all points on the sets
        boolList procFacePoint(mesh.nPoints(), false);
        forAll(set, fI)
        {
            const face& f = mesh.faces()[set[fI]];
            forAll(f, fp)
            {
                procFacePoint[f[fp]] = true;
            }
        }
        syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

        // 2. Unblock all faces on procFacePoint
        forAll(procFacePoint, pointi)
        {
            if (procFacePoint[pointi])
            {
                const labelList& pFaces = mesh.pointFaces()[pointi];
                for (const label faceI : pFaces)
                {
                    const label own = mesh.faceOwner()[faceI];

                    if (decomposition[own] != procI)
                    {
                        decomposition[own] = procI;
                        ++nChanged;
                    }

                    if (mesh.isInternalFace(faceI))
                    {
                        const label nei = mesh.faceNeighbour()[faceI];
                        if (decomposition[nei] != procI)
                        {
                            decomposition[nei] = procI;
                            ++nChanged;
                        }
                    }
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : changed decomposition on "
            << returnReduce(nChanged, sumOp<label>()) << " cells" << endl;
    }
}


// ************************************************************************* //
