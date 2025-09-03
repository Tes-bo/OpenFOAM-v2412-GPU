/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::ensightFaces

Description
    Sorting/classification of faces (2D) into corresponding ensight types.

    Some caution may be required when handling face addressing into a
    boundaryField. Since the face addressing is absolute, it will be
    necessary to work on a copy with local ids. For example,

    \code
        // Operate on a copy
        ensightFaces localPart(part);

        // Change from global faceIds to patch-local faceIds
        localPart.decrFaceIds(patchStart);

        // Can now address into boundaryField
    \endcode

    Additionally, for some uses (eg, finiteArea), the ordered
    list of faces can be used for a direct local lookup into the field
    instead of via the overall mesh face addressing.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightFaces_H
#define Foam_ensightFaces_H

#include "ensightPart.H"
#include "face.H"
#include "FixedList.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
template<class T> class InfoProxy;

/*---------------------------------------------------------------------------*\
                        Class ensightFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightFaces
:
    public ensightPart
{
public:

    // Public Data

        //- Supported ensight 'Face' element types.
        //  Must be zero-based since they are also for internal bookkeeping.
        enum elemType
        {
            TRIA3 = 0,  //!< "tria3"
            QUAD4,      //!< "quad4"
            NSIDED      //!< "nsided"
        };

        //- Number of 'Face' element types (3)
        static constexpr int nTypes = 3;

        //- The ensight 'Face' element type names
        static const char* elemNames[nTypes];


    // Static Functions

        //- The ensight element name for the specified 'Face' type
        static inline const char* key(const elemType etype) noexcept;


private:

    // Private Data

        //- The input face order, for indirect face lists
        labelList faceOrder_;

        //- List of face-flips (optional)
        boolList flipMap_;

        //- Begin/end offsets for address/flips of each element type
        FixedList<label, nTypes+1> offsets_;

        //- List of global sizes for each element type.
        //  Used temporarily for local sizes when building the element lists.
        FixedList<label, nTypes> sizes_;


    // Private Member Functions

        //- Low-level internal addition routine.
        //  \return insertion locaion
        inline label add(const elemType etype, label id, bool flip=false);

        //- Use temporarily stored sizes to redimension the element lists
        void resizeAll();


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("ensightFaces");


    // Constructors

        //- Default construct, with part index 0
        ensightFaces();

        //- Default construct, with description/partName
        explicit ensightFaces(const string& description);


    //- Destructor
    virtual ~ensightFaces() = default;


    // Member Functions

    // Access

        //- Processor-local size of all elements.
        using ensightPart::size;

        //- Processor-local size of the specified element type.
        inline label size(const elemType etype) const;

        //- Processor-local offset/size of element type.
        inline labelRange range(const elemType etype) const;

        //- The global size of all element types.
        //  This value is only meaningful after a reduce operation.
        label totalSize() const noexcept;

        //- Same as totalSize
        label total() const noexcept { return totalSize(); }

        //- The global size of the specified element type.
        //  This value is only meaningful after a reduce operation.
        inline label total(const elemType etype) const;

        //- The global sizes for each element type.
        //  This value is only meaningful after a reduce operation.
        inline const FixedList<label, nTypes>& totals() const;

        //- Processor-local sizes per element type.
        FixedList<label, nTypes> sizes() const;

        //- Processor-local face ids of all elements
        inline const labelList& faceIds() const noexcept;

        //- Processor-local face ids of the specified element type
        inline const labelUList faceIds(const elemType etype) const;

        //- Processor-local flip-map of all elements
        inline const boolList& flipMap() const;

        //- True for non-zero flip-map that spans the addresses
        inline bool usesFlipMap() const;

        //- Processor-local face order
        //- (where applicable)
        inline const labelList& faceOrder() const noexcept;

        //- Processor-local face order of specified element type
        //- (where applicable)
        inline const labelUList faceOrder(const elemType etype) const;


    // Edit

        //- Classify the face types and set the element lists.
        void classify(const UList<face>& faces);

        //- Classify face types (for a sublist) and set element lists.
        void classify(const UList<face>& faces, const labelRange& range);

        //- Classify the face types and set the element lists.
        //  The indirect addressing can be used when classifying groups of
        //  face (eg, from a faceZone etc) with an optional flipMap.
        //  The optional exclude marker can be used to skip faces on particular
        //  boundary types or regions.
        void classify
        (
            const UList<face>& faces,
            const labelUList& addr,
            const boolList& flipMap = boolList(),
            const bitSet& exclude = bitSet()
        );


        //- Clear any demand-driven data
        void clearOut();

        //- Set addressable sizes to zero, free up addressing memory.
        void clear();

        //- Sum element counts across all processes.
        void reduce();

        //- Inplace sort element lists numerically.
        void sort();


    // Advanced (use with caution)

        //- Increase face ids by specified offset value
        //  Eg, to change patch local Ids to global Ids
        void incrFaceIds(const label off);

        //- Decrease face ids by specified offset value
        //  Eg, to change global Ids to patch local Ids
        void decrFaceIds(const label off);


    // Output

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<ensightFaces> info() const noexcept { return *this; }


        //- Globally unique mesh points.
        //- Required when writing point fields.
        label uniqueMeshPoints
        (
            const polyMesh& mesh,
            labelList& uniqueMeshPointLabels,
            bool parallel
        ) const;

        //- Write information about the object as a dictionary,
        //- optionally write all element addresses
        virtual void writeDict(Ostream& os, const bool full=false) const;

        //- Write geometry, using a mesh reference
        //- No beginGeometry() marker.
        virtual void write
        (
            ensightGeoFile& os,
            const polyMesh& mesh,
            bool parallel
        ) const;
};


//- Print information via proxy
template<>
Ostream& operator<<(Ostream&, const InfoProxy<ensightFaces>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightFacesI.H"

#endif

// ************************************************************************* //
