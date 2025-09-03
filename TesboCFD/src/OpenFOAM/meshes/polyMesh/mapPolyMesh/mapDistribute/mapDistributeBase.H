/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

Class
    Foam::mapDistributeBase

Description
    Class containing processor-to-processor mapping information.

    We store mapping from the bits-to-send to the complete starting list
    (subXXXMap) and from the received bits to their location in the new
    list (constructXXXMap).

Note:
    Schedule is a list of processor pairs (one send, one receive. One of
    them will be myself) which forms a scheduled (i.e. non-buffered) exchange.
    See distribute on how to use it.
    Note2: number of items sent on one processor have to equal the number
    of items received on the other processor.

    To aid constructing these maps there are the constructors from global
    numbering, either with or without transforms.

    Constructors using compact numbering: layout is
    - all my own elements first (whether used or not)
    - followed by used-only remote elements sorted by remote processor.
    So e.g 4 procs and on proc 1 the compact
    table will first have all globalIndex.localSize() elements from proc1
    followed by used-only elements of proc0, proc2, proc3.
    The constructed mapDistributeBase sends the local elements from and
    receives the remote elements into their compact position.
    compactMap[proci] is the position of elements from proci in the compact
    map. compactMap[myProcNo()] is empty since trivial addressing.

    It rewrites the input global indices into indices into the constructed
    data.

    When constructing from components optionally a 'flip' on
    the maps can be specified. This will interpret the map
    values as index+flip, similar to e.g. faceProcAddressing. The flip
    will only be applied to fieldTypes (scalar, vector, .. triad)


SourceFiles
    mapDistributeBase.C
    mapDistributeBaseIO.C
    mapDistributeBaseTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mapDistributeBase_H
#define Foam_mapDistributeBase_H

#include "boolList.H"
#include "labelList.H"
#include "labelPair.H"
#include "Pstream.H"
#include "Map.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;
class dictionary;
class globalIndex;
class PstreamBuffers;
class mapPolyMesh;
class mapDistributeBase;

Istream& operator>>(Istream&, mapDistributeBase&);
Ostream& operator<<(Ostream&, const mapDistributeBase&);

template<>
Ostream& operator<<(Ostream&, const InfoProxy<mapDistributeBase>&);


/*---------------------------------------------------------------------------*\
                           Class mapDistributeBase Declaration
\*---------------------------------------------------------------------------*/

class mapDistributeBase
{
public:

    //- The map layout (eg, of the constructMap)
    enum class layoutTypes : char
    {
        linear,     //!< In processor-order
        localFirst  //!< Local processor first, others in linear order
    };


private:

    // Private Data

        //- Size of reconstructed data
        label constructSize_;

        //- Maps from subsetted data back to original data
        labelListList subMap_;

        //- Maps from subsetted data to new reconstructed data
        labelListList constructMap_;

        //- Whether subMap includes flip or not
        bool subHasFlip_;

        //- Whether constructMap includes flip or not
        bool constructHasFlip_;

        //- Communicator to use for parallel operations
        label comm_;

        //- Schedule
        mutable autoPtr<List<labelPair>> schedulePtr_;


protected:

    // Protected Member Functions

        //- Fatal if expected != received size
        static void checkReceivedSize
        (
            const label proci,
            const label expectedSize,
            const label receivedSize
        );

        //- Scan the maps for the max addressed index.
        //
        //  \param maps  The maps to scan
        //  \param hasFlip  True if maps has flip addressing
        //  \return max-size needed for addressing (eg, constructSize)
        static label getMappedSize
        (
            const labelListList& maps,
            const bool hasFlip
        );

        //- Construct per processor compact addressing of the global elements
        //  needed. The ones from the local processor are not included since
        //  these are always all needed.
        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelUList& elements,
            List<Map<label>>& compactMap
        ) const;

        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelListList& elements,
            List<Map<label>>& compactMap
        ) const;

        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelList& elements,
            List<Map<label>>& compactMap,
            labelList& compactStart
        );
        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelListList& elements,
            List<Map<label>>& compactMap,
            labelList& compactStart
        );

        //- Synchronize send/recv mask buffers as a 'copy' operation.
        //  Somewhat similar to Pstream::exchangeContainer
        //
        //  The masks must be properly sized by the caller, which avoids
        //  a needless all-to-all for the sizes and the sizing is already
        //  given by the maps.
        static void exchangeMasks
        (
            const UList<bitSet>& sendMasks,
            UList<bitSet>& recvMasks,
            const int tag,
            const label comm
        );

        //- Bi-direction sync of send/recv buffers using bitwise '&='
        //- combine operation.
        //
        //  The masks must be properly sized by the caller, which avoids
        //  a needless all-to-all for the sizes and the sizing is already
        //  given by the maps.
        static void unionCombineMasks
        (
            UList<bitSet>& sendMasks,
            UList<bitSet>& recvMasks,
            const int tag,
            const label comm
        );


        //- Combine field values (after any flip negation operation)
        //- into the specified mapped target locations
        template<class T, class CombineOp, class NegateOp>
        static void flipAndCombine
        (
            //! [in,out] The left of binary combine operation
            List<T>& lhs,
            //! The right of binary combine operation
            const UList<T>& rhs,
            //! The mapping indices
            const labelUList& map,
            //! Mapping indices include flip encoding
            const bool hasFlip,
            //! Binary combine operation
            const CombineOp& cop,
            //! Unary negation operation (for flipped indices)
            const NegateOp& negOp
        );

        //- Lookup field values at specified map indices and save
        //- after any flip negation operations
        template<class T, class NegateOp>
        static void accessAndFlip
        (
            //! [out] The result values
            List<T>& output,
            //! [out] The input values
            const UList<T>& values,
            //! The mapping indices
            const labelUList& map,
            //! Mapping indices include flip encoding
            const bool hasFlip,
            //! Unary negation operation (for flipped indices)
            const NegateOp& negOp
        );

        //- Lookup field values at specified indices and return
        //- after any flip negation operations
        template<class T, class NegateOp>
        static List<T> accessAndFlip
        (
            //! [out] The input values
            const UList<T>& values,
            //! The mapping indices
            const labelUList& map,
            //! Mapping indices include flip encoding
            const bool hasFlip,
            //! Unary negation operation (for flipped indices)
            const NegateOp& negOp
        );


private:

    // Private Member Functions

        //- Helper for compactData (private: filescope only!)
        //  Establishes the exact send/recv elements used after masking.
        //
        //  \param allowedLocalElems  Permissible local mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param allowedRemoteElems  Permissible remote mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param[out] sendMasks  Mask of local elements sent to procs.
        //  \param[out] recvMasks  Mask of remote elements received
        //      from procs
        //  \param tag  The message tag
        void calcCompactDataRequirements
        (
            const bitSet& allowedLocalElems,
            const bitSet& allowedRemoteElems,
            List<bitSet>& sendMasks,     // [out]
            List<bitSet>& recvMasks,     // [out]
            const int tag
        );

        //- Helper for compactLocalData (private: filescope only!)
        //  Establishes the exact send/recv elements used after masking.
        //
        //  \param allowedLocalElems  Permissible local mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param[out] sendMasks  Mask of local elements sent to procs.
        //  \param[out] recvMasks  Mask of remote elements received by proc.
        //      from procs
        //  \param tag  The message tag
        void calcCompactLocalDataRequirements
        (
            const bitSet& allowedLocalElems,
            List<bitSet>& sendMasks,    // [out]
            List<bitSet>& recvMasks,    // [out]
            const int tag
        );

        //- Helper for compactRemoteData (private: filescope only!)
        //  Establishes the exact send/recv elements used after masking.
        //
        //  \param allowedRemoteElems  Permissible remote mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param[out] sendMasks  Mask of local elements sent to procs.
        //  \param[out] recvMasks  Mask of remote elements received by proc.
        //  \param tag  The message tag
        void calcCompactRemoteDataRequirements
        (
            const bitSet& allowedRemoteElems,
            List<bitSet>& sendMasks,    // [out]
            List<bitSet>& recvMasks,    // [out]
            const int tag
        );

        //- Implementation for compact{Local,Remote}Data (private).
        //  Also renumbers the subMap/constructMap if oldToNew maps
        //  are notNull().
        //
        //  No communication
        void compactData
        (
            const UList<bitSet>& sendMasks,
            const UList<bitSet>& recvMasks,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1
        );

        //- Wrapper for compactData (private) that supplies oldToNew
        //- maps for renumbering if doRenumber is true.
        //  No communication
        void compactDataImpl
        (
            const UList<bitSet>& sendMasks,
            const UList<bitSet>& recvMasks,
            const bool doRenumber
        );


        //- Helper for renumbering compacted map elements and updating the
        //- supplied old-to-new mapping to account for the visit order of
        //- the original elements
        //
        //  \param origElements  The original elements visited (eg, meshPoints)
        //  \param[in,out] oldToNew The old-to-new mapping
        //  \param[in,out] mapElements  The map to be renumbered
        //  \param hasFlip  True if map has flip addressing
        static void renumberVisitOrder
        (
            const labelUList& origElements,
            labelList& oldToNew,
            labelListList& maps,
            const bool hasFlip
        );


public:

    // Declare name of the class and its debug switch
    ClassName("mapDistributeBase");


    // Constructors

        //- Default construct (uses worldComm)
        mapDistributeBase() noexcept;

        //- Default construct with specified communicator
        explicit mapDistributeBase(const label comm) noexcept;

        //- Copy construct
        mapDistributeBase(const mapDistributeBase& map);

        //- Move construct
        mapDistributeBase(mapDistributeBase&& map);

        //- Read construct from dictionary
        explicit mapDistributeBase
        (
            const dictionary& dict,
            const label comm = UPstream::worldComm
        );

        //- Move construct from components
        mapDistributeBase
        (
            const label constructSize,
            labelListList&& subMap,
            labelListList&& constructMap,
            const bool subHasFlip = false,
            const bool constructHasFlip = false,
            const label comm = UPstream::worldComm
        );

        //- Construct from reverse addressing: per data item the send
        //- processor and the receive processor.
        //
        //  \note data is not sorted per processor - cannot use printLayout!
        mapDistributeBase
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const label comm = UPstream::worldComm
        );

        //- Construct from list of (possibly) remote elements in globalIndex
        //- numbering (or -1).
        //
        //  Determines compact numbering (see above) and distribute map
        //  to get data into this ordering and renumbers the elements to
        //  be in compact numbering.
        mapDistributeBase
        (
            const globalIndex&,
            labelList& elements,
            List<Map<label>>& compactMap,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Special variant that works with the info sorted into bins
        //- according to local indices.
        //
        //  E.g. think cellCells where
        //  cellCells[localCellI] is a list of global cells
        mapDistributeBase
        (
            const globalIndex&,
            labelListList& cellCells,
            List<Map<label>>& compactMap,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Construct from my elements to send, targetting the specified
        //- constructMap layout
        explicit mapDistributeBase
        (
            const layoutTypes constructLayout,
            labelListList&& subMap,
            const bool subHasFlip = false,
            const bool constructHasFlip = false,
            const label comm = UPstream::worldComm
        );

        //- Construct from my elements to send.
        //  Assumes layout is my elements first followed by elements
        //  from all other processors in consecutive order.
        explicit mapDistributeBase
        (
            labelListList&& subMap,
            const bool subHasFlip = false,
            const bool constructHasFlip = false,
            const label comm = UPstream::worldComm
        );

        //- Construct from multiple maps and processor collation
        //  Assumes all local data first. Sorts contributions of maps
        //  in processor order i.e.  constructed map has all local data first.
        //  Returns
        //  - startOfLocal : per input map the start of the local data. Extends
        //                   one beyond number of maps so overall local size
        //                   is startOfLocal.last()
        //  - compactMaps  : per input map from slot position in the input map
        //                   to new slot position. (note there is no information
        //                   returned about which processor it is from)
        mapDistributeBase
        (
            const UPtrList<const mapDistributeBase>& maps,
            const labelList& localRanks,
            const label newComm,
            const labelListList& newToOldRanks, // from rank in newComm to
                                                // ranks in (old)comm
            labelList& startOfLocal,            // per map start of local data
            List<Map<label>>& compactMaps       // per map old slot to new slot
        );

        //- Construct from Istream
        explicit mapDistributeBase(Istream& is);


    // Static Functions

        //- Test for flip addressing, where flips are encoded as negative
        //- indices and non-flips are encoded as positive non-zero indices.
        //
        //  Exits early on the first detected zero or negative, which
        //  makes this more efficient than testing min(map) \< 0.
        //
        //  \note may return a false negative (ie, no flips detected)
        //     even when flip addressing is used, but the local map does not
        //     contain any flipped elements
        static bool hasFlipAddressing(const labelUList& map);

        //- Test for flip addressing, where flips are encoded as negative
        //- indices and non-flips are encoded as positive non-zero indices.
        //
        //  See notes above.
        static bool hasFlipAddressing(const labelListList& maps);

        //- Count the number of unmapped elements.
        //
        //  \param elements  The elements that are expected to be mapped
        //  \param maps  The maps to scan
        //  \param hasFlip  True if maps has flip addressing
        //  \return number of unmapped elements
        static label countUnmapped
        (
            const labelUList& elements,
            const labelListList& maps,
            const bool hasFlip
        );


    // Member Functions

    // Access

        //- Constructed data size
        label constructSize() const noexcept
        {
            return constructSize_;
        }

        //- Constructed data size
        label& constructSize() noexcept
        {
            return constructSize_;
        }

        //- From subsetted data back to original data
        const labelListList& subMap() const noexcept
        {
            return subMap_;
        }

        //- From subsetted data back to original data
        labelListList& subMap() noexcept
        {
            return subMap_;
        }

        //- From subsetted data to new reconstructed data
        const labelListList& constructMap() const noexcept
        {
            return constructMap_;
        }

        //- From subsetted data to new reconstructed data
        labelListList& constructMap() noexcept
        {
            return constructMap_;
        }

        //- Does subMap include a sign
        bool subHasFlip() const noexcept
        {
            return subHasFlip_;
        }

        //- Does subMap include a sign
        bool& subHasFlip() noexcept
        {
            return subHasFlip_;
        }

        //- Does constructMap include a sign
        bool constructHasFlip() const noexcept
        {
            return constructHasFlip_;
        }

        //- Does constructMap include a sign
        bool& constructHasFlip() noexcept
        {
            return constructHasFlip_;
        }

        //- The communicator used
        label comm() const noexcept
        {
            return comm_;
        }

        //- The number of sub-lists within the maps
        label nMaps() const noexcept
        {
            return constructMap_.size();
        }

        //- The sizes of the subMap lists
        labelList subMapSizes() const;

        //- The sizes of the constructMap lists
        labelList constructMapSizes() const;

        //- The sum of the subMap list sizes
        label subMapTotalSize() const noexcept;

        //- The sum of the constructMap list sizes
        label constructMapTotalSize() const noexcept;


    // Schedule

        //- Calculate a communication schedule. See above.
        static List<labelPair> schedule
        (
            const labelListList& subMap,
            const labelListList& constructMap,
            const int tag,  // Message tag: msgType()
            const label comm = UPstream::worldComm
        );

        //- Return a schedule. Demand driven. See above.
        const List<labelPair>& schedule() const;

        //- Return real or dummy schedule depending on the
        //- communication type
        const List<labelPair>& whichSchedule
        (
            const UPstream::commsTypes commsType
        ) const;


    // Other

        //- Reset to zero size, only retaining communicator
        void clear();

        //- Transfer the contents of the argument and annul the argument.
        void transfer(mapDistributeBase& rhs);

        //- Helper for construct from globalIndex. Renumbers element
        //  (in globalIndex numbering) into compact indices.
        static label renumber
        (
            const globalIndex&,
            const label comm,
            const List<Map<label>>& compactMap,
            const label globalElement
        );

        //- Helper for renumbering the (compacted) map elements
        //- using the supplied old-to-new mapping.
        //  Only compacts the maps, does not change the local layout.
        //
        //  \param[in,out] mapElements  The map to be renumbered
        //  \param oldToNew The old-to-new mapping
        //  \param hasFlip  True if map has flip addressing
        //
        //  \return max-size needed for new addressing (eg, constructSize)
        static label renumberMap
        (
            labelListList& mapElements,
            const labelUList& oldToNew,
            const bool hasFlip
        );

        //- Helper for renumbering the map elements. Assumes local elements
        //- are first, followed by any remote. Local elements get offset,
        //- remote elements are mapped.
        //
        //  \param[in,out] map  The map to be renumbered
        //  \param localSize    elements < localSize are offset
        //  \param offset       offset
        //  \param cMap         map for non-local elements
        //  \param hasFlip      True if map has flip addressing
        //
        //  \return max-size needed for new addressing (eg, constructSize)
        static label renumberMap
        (
            labelList& map,
            const label localSize,
            const label offset,
            const Map<label>& cMap,
            const bool hasFlip
        );

        //- Helper for a list of maps. Calls above renumberMap for all elements
        static label renumberMap
        (
            labelListList& mapElements,
            const label localSize,
            const label offset,
            const Map<label>& cMap,
            const bool hasFlip
        );


    // Compaction

        //- Compact send/receive maps based on selection of
        //- originating local (send) elements.
        //  Determines and removes the correspondingly unneeded elements
        //  in the send/receive maps.
        //  Only compacts the maps, does not change the local layout.
        //
        //  \param allowedLocalElems  Permissible local mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param tag  The message tag
        //  \param doRenumber  Apply oldToNew internally to renumber
        //     entries (uses renumberMap) and adjust the constructSize
        //
        //  \note generally preferable to compact based on remote data
        //      (ie, the actual receiver)
        void compactLocalData
        (
            const bitSet& allowedLocalElems,
            const int tag = UPstream::msgType(),
            const bool doRenumber = false
        );

        //- Compact send/receive maps based on selection of
        //- remote (receive) elements.
        //  Determines and removes the correspondingly unneeded elements
        //  in the send/receive maps.
        //  Only compacts the maps, does not change the local layout.
        //
        //  \param allowedRemoteElems  Permissible remote mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param tag  The message tag
        //  \param doRenumber  Apply oldToNew internally to renumber
        //     entries (uses renumberMap) and adjust the constructSize
        void compactRemoteData
        (
            const bitSet& allowedRemoteElems,
            const int tag = UPstream::msgType(),
            const bool doRenumber = false
        );


        //- Compact send/receive maps based on selection of
        //- originating local (send) elements.
        //- Returns compaction mappings for subMap and constructMap.
        //
        //  \param allowedLocalElems  Permissible local mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param[out] oldToNewSub  Old-to-new mapping: subMap
        //     Disabled if labelList::null(),
        //  \param[out] oldToNewConstruct  Old-to-new mapping: constructMap
        //     Disabled if labelList::null(),
        //  \param localSize  The max index for subMap (-1: auto-detect)
        //  \param tag  The message tag
        //
        //  \note Applies oldToNew to renumber entries
        //     (uses renumberMap) and adjust constructSize
        //
        //  \note generally preferable to compact based on remote data
        //      (ie, the actual receiver)
        void compactLocalData
        (
            const bitSet& allowedLocalElems,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1,
            const int tag = UPstream::msgType()
        );

        //- Compact send/receive maps based on selection of
        //- remote (receive) elements.
        //- Returns compaction mappings for subMap and constructMap.
        //
        //  \param allowedRemoteElems  Permissible remote mapped elements
        //     (true/false). Can be longer/shorter than actual number
        //     of mapped elements.
        //  \param[out] oldToNewSub  Old-to-new mapping: subMap
        //     Disabled if labelList::null(),
        //  \param[out] oldToNewConstruct  Old-to-new mapping: constructMap
        //     Disabled if labelList::null(),
        //  \param localSize  The max index for subMap (-1: auto-detect)
        //  \param tag  The message tag
        //
        //  \note Applies oldToNew to renumber entries
        //     (uses renumberMap) and adjust constructSize
        void compactRemoteData
        (
            const bitSet& allowedRemoteElems,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1,
            const int tag = UPstream::msgType()
        );


        //- Compact send/receive maps based on selection of
        //- originating local (send) and remote (receive) elements.
        //
        //  The resulting compact numbering:
        //  - \c subMap (and \c oldToNewSub) :
        //    will follow the original ordering of \c localElements.
        //  - \c constructMap (and \c oldToNewConstruct) :
        //    will follow the original ordering of \c remoteElements.
        //  .
        //  \warning ill-defined behaviour if \c localElements
        //  or \c remoteElements contains duplicates.
        void compactData
        (
            const labelUList& localElements,
            const labelUList& remoteElements,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1,
            const int tag = UPstream::msgType()
        );

        //- Compact send/receive maps based on selection of
        //- originating local (send) elements.
        //
        //  The resulting compact numbering:
        //  - \c subMap (and \c oldToNewSub) :
        //    will follow the original ordering of \c localElements.
        //  - \c constructMap (and \c oldToNewConstruct) :
        //    numbered in simple ascending order.
        //  .
        //  \warning ill-defined behaviour if \c localElements
        //  contains duplicates.
        //
        //  \note generally preferable to compact based on remote data
        //      (ie, the actual receiver)
        void compactLocalData
        (
            const labelUList& localElements,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1,
            const int tag = UPstream::msgType()
        );

        //- Compact send/receive maps based on selection of
        //- remote (receive) elements.
        //
        //  The resulting compact numbering:
        //  - \c subMap (and \c oldToNewSub) :
        //    numbered in simple ascending order.
        //  - \c constructMap (and \c oldToNewConstruct) :
        //    will follow the original ordering of \c remoteElements.
        //  .
        //  \warning ill-defined behaviour if \c remoteElements
        //  contains duplicates.
        void compactRemoteData
        (
            const labelUList& remoteElements,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const label localSize = -1,
            const int tag = UPstream::msgType()
        );


    // Distribute

        //- Distribute combine data with specified combine operation
        //- and negate operator (for flips).
        //
        //  If multiple processors write to same position,
        //  contributions are added using the combine cop.
        //
        //  \note schedule only used for UPstream::commsTypes::scheduled,
        //        others just use send-to-all, receive-from-all.
        template<class T, class CombineOp, class NegateOp>
        static void distribute
        (
            const UPstream::commsTypes commsType,
            const List<labelPair>& schedule,
            const label constructSize,
            const labelListList& subMap,
            const bool subHasFlip,
            const labelListList& constructMap,
            const bool constructHasFlip,
            List<T>& field,
            const T& nullValue,
            const CombineOp& cop,
            const NegateOp& negOp,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Distribute assign data with specified negate operator (for flips).
        //- Uses assignment for combine operation.
        //
        //  \note schedule only used for UPstream::commsTypes::scheduled,
        //        others just use send-to-all, receive-from-all.
        template<class T, class NegateOp>
        static void distribute
        (
            const UPstream::commsTypes commsType,
            const List<labelPair>& schedule,
            const label constructSize,
            const labelListList& subMap,
            const bool subHasFlip,
            const labelListList& constructMap,
            const bool constructHasFlip,
            List<T>& field,
            const NegateOp& negOp,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Distribute (simpler interface)

            //- Distribute List data using default commsType
            //- and the default flip/negate operator.
            template<class T>
            void distribute
            (
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute DynamicList data using default commsType
            //- and the default flip/negate operator.
            template<class T>
            void distribute
            (
                DynamicList<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute List data using specified commsType
            //- and the default flip/negate operator.
            template<class T>
            void distribute
            (
                const UPstream::commsTypes commsType,
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute DynamicList data using specified commsType
            //- and the default flip/negate operator.
            template<class T>
            void distribute
            (
                const UPstream::commsTypes commsType,
                DynamicList<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute data using default commsType
            //- and the specified negate operator (for flips).
            template<class T, class NegateOp>
            void distribute
            (
                List<T>& values,
                const NegateOp& negOp,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute data using specified commsType
            //- and the specified negate operator (for flips).
            //  Accepts a nullValue to initialize unmapped elements
            //  (ie, when the constructSize is larger than the number of
            //  mapped elements).
            template<class T, class NegateOp>
            void distribute
            (
                const UPstream::commsTypes commsType,
                List<T>& values,
                const NegateOp& negOp,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute data using specified commsType
            //- and the specified negate operator (for flips).
            //  Accepts a nullValue to initialize unmapped elements
            //  (ie, when the constructSize is larger than the number of
            //  mapped elements).
            template<class T, class NegateOp>
            void distribute
            (
                const UPstream::commsTypes commsType,
                const T& nullValue,
                List<T>& values,
                const NegateOp& negOp,
                const int tag = UPstream::msgType()
            ) const;


    // Reverse Distribute (simpler interface)

            //- Reverse distribute data using default commsType
            //- and the default flip/negate operator
            template<class T>
            void reverseDistribute
            (
                const label constructSize,
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using default commsType
            //- and the default flip/negate operator.
            //  Accepts a nullValue to initialize unmapped elements
            //  (ie, when the constructSize is larger than the subMap).
            template<class T>
            void reverseDistribute
            (
                const label constructSize,
                const T& nullValue,
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using specified commsType
            //- and the default flip/negate operator
            template<class T>
            void reverseDistribute
            (
                const UPstream::commsTypes commsType,
                const label constructSize,
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using specified commsType
            //- and the specified flip/negate operator.
            template<class T, class NegateOp>
            void reverseDistribute
            (
                const UPstream::commsTypes commsType,
                const label constructSize,
                List<T>& values,
                const NegateOp& negOp,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using specified commsType
            //- and the default flip/negate operator.
            //  Accepts a nullValue to initialize unmapped elements
            //  (ie, when the constructSize is larger than the subMap).
            template<class T>
            void reverseDistribute
            (
                const UPstream::commsTypes commsType,
                const label constructSize,
                const T& nullValue,
                List<T>& values,
                const int tag = UPstream::msgType()
            ) const;


    // Send/Receive

            //- Do all sends using PstreamBuffers
            template<class T>
            void send(PstreamBuffers& pBufs, const List<T>& field) const;

            //- Do all receives using PstreamBuffers
            template<class T>
            void receive(PstreamBuffers& pBufs, List<T>& field) const;


        // Low-level. Non-blocking. TBD: receive map is usually contiguous
        // - receive into slices

            //- Start all sends and receives (always non-blocking)
            template<class T, class negateOp>
            static void send
            (
                const labelListList& subMap,
                const bool subHasFlip,
                const labelListList& constructMap,
                const bool constructHasFlip,
                const UList<T>& field,
                labelRange& sendRequests,       //!< The send requests range
                PtrList<List<T>>& sendFields,   //!< The send fields storage
                labelRange& recvRequests,       //!< The recv requests range
                PtrList<List<T>>& recvFields,   //!< The recv fields storage
                const negateOp& negOp,
                const int tag,
                const label comm
            );

            //- Start all sends and receives (always non-blocking)
            template<class T>
            void send
            (
                const UList<T>& field,
                labelRange& sendRequests,       //!< The send requests range
                PtrList<List<T>>& sendFields,   //!< The send fields storage
                labelRange& recvRequests,       //!< The recv requests range
                PtrList<List<T>>& recvFields,   //!< The recv fields storage
                const int tag = UPstream::msgType()
            ) const;

            //- Wait for (receive) requests to finish and consume
            template<class T, class CombineOp, class negateOp>
            static void receive
            (
                const label constructSize,
                const labelListList& constructMap,
                const bool constructHasFlip,
                const labelRange& requests,
                const UPtrList<List<T>>& recvFields,
                List<T>& field,
                const CombineOp& cop,
                const negateOp& negOp,
                const int tag,
                const label comm
            );

            //- Wait for (receive) requests to finish and consume
            template<class T>
            void receive
            (
                const labelRange& requests,
                const UPtrList<List<T>>& recvFields,
                List<T>& field,
                const int tag = UPstream::msgType()
            ) const;


    // Other

        //- Debug: print layout. Can only be used on maps with sorted
        //  storage (local data first, then non-local data)
        void printLayout(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const mapDistributeBase& rhs);

        //- Move assignment
        void operator=(mapDistributeBase&& rhs);


    // IOstream Operators

        //- Return info proxy,
        //- used to print summary information to a stream
        InfoProxy<mapDistributeBase> info() const noexcept
        {
            return *this;
        }

        //- Read entries from dictionary format
        void readDict(const dictionary& dict);

        //- Write entries in dictionary format
        void writeEntries(Ostream& os) const;

        //- Read plain content (not dictionary) from Istream
        friend Istream& operator>>(Istream&, mapDistributeBase&);

        //- Write plain content (not dictionary) to Ostream
        friend Ostream& operator<<(Ostream&, const mapDistributeBase&);


    // Housekeeping

        //- No correction for topo change
        void updateMesh(const mapPolyMesh&)
        {
            NotImplemented;
        }

        //- OpenFOAM-v2112 and earlier naming for compactRemoteData()
        //- using boolList.
        void compact
        (
            const boolList& remoteElemUsed,
            const int tag = UPstream::msgType()
        );

        //- OpenFOAM-v2112 and earlier naming for compactRemoteData().
        //- using boolList.
        void compact
        (
            const boolList& remoteElemUsed,
            const label localSize,
            labelList& oldToNewSub,
            labelList& oldToNewConstruct,
            const int tag = UPstream::msgType()
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapDistributeBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
