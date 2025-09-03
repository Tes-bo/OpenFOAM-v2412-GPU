/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "ensightOutput.H"
#include "ensightPTraits.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class LabelListListType>
void Foam::ensightOutput::Detail::writeLabelListList
(
    ensightGeoFile& os,
    const LabelListListType& listOfLists,
    const label pointOffset
)
{
    const label off = (pointOffset + 1);  // 1-based for Ensight

    forAll(listOfLists, listi)
    {
        for (const label pointi : listOfLists[listi])
        {
            os.write(pointi + off);
        }
        os.newline();  // One list (cell/faces) per line (ASCII)
    }
}


template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::Detail::copyComponent
(
    const FieldContainer<Type>& input,
    const direction cmpt,
    UList<float>& cmptBuffer
)
{
    if (cmptBuffer.size() < input.size())
    {
        FatalErrorInFunction
            << "Component buffer too small: "
            << cmptBuffer.size() << " < " << input.size() << nl
            << exit(FatalError);
    }

    auto iter = cmptBuffer.begin();

    if (std::is_same<float, typename pTraits<Type>::cmptType>::value)
    {
        // Direct copy
        for (const Type& val : input)
        {
            *iter = component(val, cmpt);
            ++iter;
        }
    }
    else
    {
        // Copy with narrowing
        for (const Type& val : input)
        {
            *iter = narrowFloat(component(val, cmpt));
            ++iter;
        }
    }
}


template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::Detail::writeFieldComponents
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const char* key,
    const FieldContainer<Type>& fld,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    const label localSize = fld.size();

    // Gather sizes (offsets irrelevant)
    const globalIndex procAddr
    (
        parallel
      ? globalIndex(globalIndex::gatherOnly{}, localSize)
      : globalIndex(globalIndex::gatherNone{}, localSize)
    );

    if (Pstream::master() && key)
    {
        os.writeKeyword(key);
    }

    if (Pstream::master())
    {
        // Scratch buffer:
        // - allocate enough space to process an individual rank
        // - potentially enough to process multiple ranks before writing
        // - permit use of the full buffer capacity

        // Buffer size needed for an individual rank
        const label minSize(max(localSize, procAddr.maxSize()));

        // Maximum off-processor transfer size
        const label maxSize =
        (
            (ensightOutput::maxChunk_ > 0)
          ? min
            (
                static_cast<label>(ensightOutput::maxChunk_),
                (procAddr.totalSize() - localSize)
            )
          : scratch.capacity()
        );

        scratch.resize_nocopy
        (
            max(max(minSize, maxSize), scratch.capacity())
        );

        if (Pstream::master() && debug > 1)
        {
            Info<< "ensight";
            if (key)
            {
                Info<< " (" << key << ')';
            }

            Info<< " total-size:" << procAddr.totalSize()
                << " buf-size:" << scratch.size() << "/" << scratch.capacity()
                << " any-proc:" << minSize
                << " off-proc:" << (procAddr.totalSize() - localSize) << endl;

            Info<< "proc-sends: (";

            label nPending = localSize;

            Info<< (localSize ? '0' : '_');

            // Receive others, writing as needed
            for (const label proci : procAddr.subProcs())
            {
                const label procSize = procAddr.localSize(proci);

                if (procSize)
                {
                    if (nPending + procSize > scratch.size())
                    {
                        // Flush buffer
                        nPending = 0;
                        Info<< ") (";
                    }
                    else
                    {
                        Info<< ' ';
                    }

                    Info<< proci;
                    nPending += procSize;
                }
            }

            Info<< ')' << endl;
        }

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            // Master
            copyComponent(fld, cmpt, scratch);
            label nPending = localSize;

            // Receive others, writing as needed
            for (const label proci : procAddr.subProcs())
            {
                const label procSize = procAddr.localSize(proci);

                if (procSize)
                {
                    if (nPending + procSize > scratch.size())
                    {
                        // Flush buffer
                        os.writeList(SubList<float>(scratch, nPending));
                        nPending = 0;
                    }

                    SubList<float> slot(scratch, procSize, nPending);
                    nPending += procSize;

                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        proci,
                        slot.data_bytes(),
                        slot.size_bytes()
                    );
                }
            }

            if (nPending)
            {
                // Flush buffer
                os.writeList(SubList<float>(scratch, nPending));
            }
        }
    }
    else if (parallel && localSize)
    {
        scratch.resize_nocopy(localSize);

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            copyComponent(fld, cmpt, scratch);

            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                scratch.cdata_bytes(),
                scratch.size_bytes()
            );
        }
    }
}


template<template<typename> class FieldContainer>
bool Foam::ensightOutput::Detail::writeCoordinates
(
    ensightGeoFile& os,
    const label partId,
    const word& partName,
    const label nPoints,
    const FieldContainer<Foam::point>& fld,
    bool parallel
)
{
    if (Pstream::master())
    {
        os.beginPart(partId, partName);
        os.beginCoordinates(nPoints);
    }

    bool ok = (Pstream::master() && (nPoints > 0));
    if (parallel)
    {
        Pstream::broadcast(ok);
    }

    if (ok)
    {
        ensightOutput::floatBufferType scratch;
        ensightOutput::Detail::writeFieldComponents
        (
            scratch,
            os,
            nullptr,  // (no element type)
            fld,
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceSubField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Need geometry and field. part.total() is already reduced
    const bool good =
    (
        parallel
      ? (part.total() && returnReduceOr(fld.size()))
      : (part.size() && fld.size())
    );

    if (!good)
    {
        return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        // Write elements of this type
        if (parallel ? part.total(etype) : part.size(etype))
        {
            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightFaces::key(etype),
                SubField<Type>(fld, part.range(etype)),
                parallel
            );
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceLocalField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Need geometry and field. part.total() is already reduced
    const bool good =
    (
        parallel
      ? (part.total() && returnReduceOr(fld.size()))
      : (part.size() && fld.size())
    );

    if (!good)
    {
        return false;
    }


    bool validAddressing = (part.size() == part.faceOrder().size());

    if (parallel)
    {
        Pstream::reduceOr(validAddressing);
    }

    if (!validAddressing)
    {
        FatalErrorInFunction
            << "Called without faceOrder having been set" << nl
            << exit(FatalError);
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        // Write elements of this type
        if (parallel ? part.total(etype) : part.size(etype))
        {
            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightFaces::key(etype),
                UIndirectList<Type>(fld, part.faceOrder(etype)),
                parallel
            );
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const Field<Type>& fld,
    const ensightCells& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Need geometry and field. part.total() is already reduced
    const bool good =
    (
        parallel
      ? (part.total() && returnReduceOr(fld.size()))
      : (part.size() && fld.size())
    );

    if (!good)
    {
        return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        // Write elements of this type
        if (parallel ? part.total(etype) : part.size(etype))
        {
            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightCells::key(etype),
                UIndirectList<Type>(fld, part.cellIds(etype)),
                parallel
            );
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Need geometry and field. part.total() is already reduced
    const bool good =
    (
        parallel
      ? (part.total() && returnReduceOr(fld.size()))
      : (part.size() && fld.size())
    );

    if (!good)
    {
        return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        // Write elements of this type
        if (parallel ? part.total(etype) : part.size(etype))
        {
            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightFaces::key(etype),
                UIndirectList<Type>(fld, part.faceIds(etype)),
                parallel
            );
        }
    }

    return true;
}


// ************************************************************************* //
