/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "cyclicPeriodicAMIPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// For debugging
#include "OBJstream.H"
#include "PatchTools.H"
#include "Time.H"
// Note: cannot use vtkSurfaceWriter here - circular linkage
// but foamVtkSurfaceWriter (vtk::surfaceWriter) would be okay.
//
// #include "foamVtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPeriodicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPeriodicAMIPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        cyclicPeriodicAMIPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cyclicPeriodicAMIPolyPatch::syncTransforms() const
{
    if (owner())
    {
        // At this point we guarantee that the transformations have been
        // updated. There is one particular case, where if the periodic patch
        // locally has zero faces then its transformation will not be set. We
        // need to synchronise the transforms over the zero-sized patches as
        // well.
        //
        // We can't put the logic into cyclicPolyPatch as
        // processorCyclicPolyPatch uses cyclicPolyPatch functionality.
        // Synchronisation will fail because processor-type patches do not exist
        // on all processors.
        //
        // The use in cyclicPeriodicAMI is special; we use the patch even
        // though we have no faces in it. Ideally, in future, the transformation
        // logic will be abstracted out, and we won't need a periodic patch
        // here. Until then, we can just fix the behaviour of the zero-sized
        // coupled patches here

        // Get the periodic patch
        const coupledPolyPatch& periodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        // If there are any zero-sized periodic patches
        if (returnReduceOr(size() && !periodicPatch.size()))
        {
            if (periodicPatch.separation().size() > 1)
            {
                FatalErrorInFunction
                    << "Periodic patch " << periodicPatchName_
                    << " has non-uniform separation vector "
                    << periodicPatch.separation()
                    << "This is not allowed inside " << type()
                    << " patch " << name()
                    << exit(FatalError);
            }

            if (periodicPatch.forwardT().size() > 1)
            {
                FatalErrorInFunction
                    << "Periodic patch " << periodicPatchName_
                    << " has non-uniform transformation tensor "
                    << periodicPatch.forwardT()
                    << "This is not allowed inside " << type()
                    << " patch " << name()
                    << exit(FatalError);
            }

            // Note that zero-sized patches will have zero-sized fields for the
            // separation vector, forward and reverse transforms. These need
            // replacing with the transformations from other processors.

            // Parallel in this context refers to a parallel transformation,
            // rather than a rotational transformation.

            // Note that a cyclic with zero faces is considered parallel so
            // explicitly check for that.

            if
            (
                returnReduceOr
                (
                    periodicPatch.size() && periodicPatch.parallel()
                )
            )
            {
                // Sync a list of separation vectors
                List<vectorField> sep(Pstream::nProcs());
                sep[Pstream::myProcNo()] = periodicPatch.separation();
                Pstream::allGatherList(sep);

                List<boolList> coll(Pstream::nProcs());
                coll[Pstream::myProcNo()] = periodicPatch.collocated();
                Pstream::allGatherList(coll);

                // If locally we have zero faces pick the first one that has a
                // separation vector
                if (!periodicPatch.size())
                {
                    forAll(sep, procI)
                    {
                        if (sep[procI].size())
                        {
                            const_cast<vectorField&>
                            (
                                periodicPatch.separation()
                            ) = sep[procI];
                            const_cast<boolList&>
                            (
                                periodicPatch.collocated()
                            ) = coll[procI];

                            break;
                        }
                    }
                }
            }
            else
            {
                // Sync a list of forward and reverse transforms
                List<tensorField> forwardT(Pstream::nProcs());
                forwardT[Pstream::myProcNo()] = periodicPatch.forwardT();
                Pstream::allGatherList(forwardT);

                List<tensorField> reverseT(Pstream::nProcs());
                reverseT[Pstream::myProcNo()] = periodicPatch.reverseT();
                Pstream::allGatherList(reverseT);

                // If locally we have zero faces pick the first one that has a
                // transformation vector
                if (!periodicPatch.size())
                {
                    forAll(forwardT, procI)
                    {
                        if (forwardT[procI].size())
                        {
                            const_cast<tensorField&>
                            (
                                periodicPatch.forwardT()
                            ) = forwardT[procI];
                            const_cast<tensorField&>
                            (
                                periodicPatch.reverseT()
                            ) = reverseT[procI];

                            break;
                        }
                    }
                }
            }
        }
    }
}


void Foam::cyclicPeriodicAMIPolyPatch::writeOBJ
(
    const primitivePatch& p,
    OBJstream& str
) const
{
    // Collect faces and points
    pointField allPoints;
    faceList allFaces;
    PatchTools::gatherAndMerge
    (
        -1.0,           // do not merge points
        p,
        allPoints,
        allFaces
    );

    if (Pstream::master())
    {
        // Write base geometry
        str.write(allFaces, allPoints);
    }
}


void Foam::cyclicPeriodicAMIPolyPatch::resetAMI() const
{
    if (owner())
    {
        // Get the periodic patch
        const coupledPolyPatch& periodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        // Synchronise the transforms
        syncTransforms();

        // Create copies of both patches' points, transformed to the owner
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(neighbPatch().localPoints());
        transformPosition(nbrPoints0);

        // Reset the stored number of periodic transformations to a lower
        // absolute value if possible
        if (nSectors_ > 0)
        {
            if (nTransforms_ > nSectors_/2)
            {
                nTransforms_ -= nSectors_;
            }
            else if (nTransforms_ < - nSectors_/2)
            {
                nTransforms_ += nSectors_;
            }
        }

        // Apply the stored number of periodic transforms
        for (label i = 0; i < nTransforms_; ++ i)
        {
            periodicPatch.transformPosition(thisPoints0);
        }
        for (label i = 0; i > nTransforms_; -- i)
        {
            periodicPatch.transformPosition(nbrPoints0);
        }

        autoPtr<OBJstream> ownStr;
        autoPtr<OBJstream> neiStr;
        if (debug)
        {
            const Time& runTime = boundaryMesh().mesh().time();

            const fileName dir(runTime.globalPath());
            const fileName postfix("_" + runTime.timeName()+"_expanded.obj");

            ownStr.reset(new OBJstream(dir/name() + postfix));
            neiStr.reset(new OBJstream(dir/neighbPatch().name() + postfix));

            InfoInFunction
                << "patch:" << name()
                << " writing accumulated AMI to "
                << ownStr().name().relative(dir)
                << " and " << neiStr().name().relative(dir) << endl;
        }

        // Create another copy
        pointField thisPoints(thisPoints0);
        pointField nbrPoints(nbrPoints0);

        // Create patches for all the points

        // Source patch at initial location
        const primitivePatch thisPatch0
        (
            SubList<face>(localFaces(), size()),
            thisPoints0
        );
        // Source patch that gets moved
        primitivePatch thisPatch
        (
            SubList<face>(localFaces(), size()),
            thisPoints
        );
        // Target patch at initial location
        const primitivePatch nbrPatch0
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints0
        );
        // Target patch that gets moved
        primitivePatch nbrPatch
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints
        );

        // Construct a new AMI interpolation between the initial patch locations
        AMIPtr_->setRequireMatch(false);
        AMIPtr_->calculate(thisPatch0, nbrPatch0, surfPtr());

        // Number of geometry replications
        label iter(0);
        label nTransformsOld(nTransforms_);

        if (ownStr)
        {
            writeOBJ(thisPatch0, *ownStr);
        }
        if (neiStr)
        {
            writeOBJ(nbrPatch0, *neiStr);
        }


        // Weight sum averages
        scalar srcSum(gAverage(AMIPtr_->srcWeightsSum()));
        scalar tgtSum(gAverage(AMIPtr_->tgtWeightsSum()));
        // Direction (or rather side of AMI : this or nbr patch) of
        // geometry replication
        bool direction = nTransforms_ >= 0;

        // Increase in the source weight sum for the last iteration in the
        // opposite direction. If the current increase is less than this, the
        // direction (= side of AMI to transform) is reversed.
        // We switch the side to replicate instead of reversing the transform
        // since at the coupledPolyPatch level there is no
        // 'reverseTransformPosition' functionality.
        scalar srcSumDiff = 0;

        DebugInFunction
            << "patch:" << name()
            << " srcSum:" << srcSum
            << " tgtSum:" << tgtSum
            << " direction:" << direction
            << endl;

        // Loop, replicating the geometry
        while
        (
            (iter < maxIter_)
         && (
                (1 - srcSum > matchTolerance())
             || (1 - tgtSum > matchTolerance())
            )
        )
        {
            if (direction)
            {
                periodicPatch.transformPosition(thisPoints);

                DebugInFunction
                    << "patch:" << name()
                    << " moving this side from:"
                    << gAverage(thisPatch.points())
                    << " to:" << gAverage(thisPoints) << endl;

                thisPatch.movePoints(thisPoints);

                DebugInFunction
                    << "patch:" << name()
                    << " appending weights with untransformed slave side"
                    << endl;

                AMIPtr_->append(thisPatch, nbrPatch0);

                if (ownStr)
                {
                    writeOBJ(thisPatch, *ownStr);
                }
            }
            else
            {
                periodicPatch.transformPosition(nbrPoints);

                DebugInFunction
                    << "patch:" << name()
                    << " moving neighbour side from:"
                    << gAverage(nbrPatch.points())
                    << " to:" << gAverage(nbrPoints) << endl;

                nbrPatch.movePoints(nbrPoints);

                AMIPtr_->append(thisPatch0, nbrPatch);

                if (neiStr)
                {
                    writeOBJ(nbrPatch, *neiStr);
                }
            }

            const scalar srcSumNew = gAverage(AMIPtr_->srcWeightsSum());
            const scalar srcSumDiffNew = srcSumNew - srcSum;

            if (srcSumDiffNew < srcSumDiff || srcSumDiffNew < SMALL)
            {
                direction = !direction;

                srcSumDiff = srcSumDiffNew;
            }

            srcSum = srcSumNew;
            tgtSum = gAverage(AMIPtr_->tgtWeightsSum());

            nTransforms_ += direction ? +1 : -1;

            ++iter;

            DebugInFunction
                << "patch:" << name()
                << " iteration:" << iter
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << " direction:" << direction
                << endl;
        }


        // Close debug streams
        ownStr.reset(nullptr);
        neiStr.reset(nullptr);


        // Average the number of transformations
        nTransforms_ = (nTransforms_ + nTransformsOld)/2;

        // Check that the match is complete
        if (iter == maxIter_)
        {
            // The matching algorithm has exited without getting the
            // srcSum and tgtSum above 1. This can happen because
            // - of an incorrect setup
            // - or because of non-exact meshes and truncation errors
            //   (transformation, accumulation of cutting errors)
            // so for now this situation is flagged as a SeriousError instead of
            // a FatalError since the default matchTolerance is quite strict
            // (0.001) and can get triggered far into the simulation.
            SeriousErrorInFunction
                << "Patches " << name() << " and " << neighbPatch().name()
                << " do not couple to within a tolerance of "
                << matchTolerance()
                << " when transformed according to the periodic patch "
                << periodicPatch.name() << "." << nl
                << "The current sum of weights are for owner " << name()
                << " : " << srcSum << " and for neighbour "
                << neighbPatch().name() << " : " << tgtSum << nl
                << "This is only acceptable during post-processing"
                << "; not during running. Improve your mesh or increase"
                << " the 'matchTolerance' setting in the patch specification."
                << endl;
        }

        // Check that both patches have replicated an integer number of times
        if
        (
            mag(srcSum - floor(srcSum + 0.5)) > srcSum*matchTolerance()
         || mag(tgtSum - floor(tgtSum + 0.5)) > tgtSum*matchTolerance()
        )
        {
            // This condition is currently enforced until there is more
            // experience with the matching algorithm and numerics.
            // This check means that e.g. different numbers of stator and
            // rotor partitions are not allowed.
            // Again see the section above about tolerances.
            SeriousErrorInFunction
                << "Patches " << name() << " and " << neighbPatch().name()
                << " do not overlap an integer number of times when transformed"
                << " according to the periodic patch "
                << periodicPatch.name() << "." << nl
                << "The current matchTolerance : " << matchTolerance()
                << ", sum of owner weights : " << srcSum
                << ", sum of neighbour weights : " << tgtSum
                << "." << nl
                << "This is only acceptable during post-processing"
                << "; not during running. Improve your mesh or increase"
                << " the 'matchTolerance' setting in the patch specification."
                << endl;
        }

        // Print some statistics
        if (returnReduceOr(size()))
        {
            scalarField srcWghtSum(size(), Zero);
            forAll(srcWghtSum, faceI)
            {
                srcWghtSum[faceI] = sum(AMIPtr_->srcWeights()[faceI]);
            }
            scalarField tgtWghtSum(neighbPatch().size(), Zero);
            forAll(tgtWghtSum, faceI)
            {
                tgtWghtSum[faceI] = sum(AMIPtr_->tgtWeights()[faceI]);
            }

            Info<< indent
                << "AMI: Patch " << name()
                << " sum(weights)"
                << " min:" << gMin(srcWghtSum)
                << " max:" << gMax(srcWghtSum)
                << " average:" << gAverage(srcWghtSum) << nl;
            Info<< indent
                << "AMI: Patch " << neighbPatch().name()
                << " sum(weights)"
                << " min:" << gMin(tgtWghtSum)
                << " max:" << gMax(tgtWghtSum)
                << " average:" << gAverage(tgtWghtSum) << nl;
        }
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    cyclicAMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        transform,
        faceAreaWeightAMI::typeName
    ),
    nTransforms_(0),
    nSectors_(0),
    maxIter_(36)
{
    AMIPtr_->setRequireMatch(false);
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch
    (
        name,
        dict,
        index,
        bm,
        patchType,
        faceAreaWeightAMI::typeName
    ),
    nTransforms_(dict.getOrDefault<label>("nTransforms", 0)),
    nSectors_(dict.getOrDefault<label>("nSectors", 0)),
    maxIter_(dict.getOrDefault<label>("maxIter", 36))
{
    AMIPtr_->setRequireMatch(false);
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{
    AMIPtr_->setRequireMatch(false);
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{
    AMIPtr_->setRequireMatch(false);
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{
    AMIPtr_->setRequireMatch(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::~cyclicPeriodicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicPeriodicAMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntryIfDifferent<label>("nTransforms", 0, nTransforms_);
    os.writeEntryIfDifferent<label>("nSectors", 0, nSectors_);
    os.writeEntryIfDifferent<label>("maxIter", 36, maxIter_);
}


// ************************************************************************* //
