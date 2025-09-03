/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::ensightCells

Description
    Sorting/classification of cells (3D) into corresponding ensight element
    types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightCells_H
#define Foam_ensightCells_H

#include "ensightPart.H"
#include "FixedList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;
class boundBox;
class polyMesh;
template<class T> class InfoProxy;

/*---------------------------------------------------------------------------*\
                        Class ensightCells Declaration
\*---------------------------------------------------------------------------*/

class ensightCells
:
    public ensightPart
{
public:

    // Public Data

        //- Supported ensight 'Cell' element types
        //  Must be zero-based since they are also for internal bookkeeping.
        enum elemType
        {
            TETRA4 = 0, //!< "tetra4"
            PYRAMID5,   //!< "pyramid5"
            PENTA6,     //!< "penta6"
            HEXA8,      //!< "hexa8"
            NFACED      //!< "nfaced"
        };

        //- Number of 'Cell' element types (5)
        static constexpr int nTypes = 5;

        //- The ensight 'Cell' element type names
        static const char* elemNames[nTypes];


    // Static Functions

        //- The ensight element name for the specified 'Cell' type
        inline static const char* key(const elemType etype) noexcept;


private:

    // Private Data

        //- Manifold cells detected
        bool manifold_;

        //- Begin/end offsets for address of each element type
        FixedList<label, nTypes+1> offsets_;

        //- List of global sizes for each element type.
        //  Used temporarily for local sizes when building the element lists.
        FixedList<label, nTypes> sizes_;


    // Private Member Functions

        //- Low-level internal addition routine
        //  \return insertion locaion
        inline label add(const elemType etype, label id);

        //- Use temporarily stored sizes to redimension the element lists
        void resizeAll();

        //- Classify cell types, set element lists for selection
        //  (implementation)
        template<class Addressing>
        void classifyImpl(const polyMesh& mesh, const Addressing& cellIds);


        label meshPointMapppings
        (
            const polyMesh& mesh,
            labelList& pointToGlobal,  // Can also be labelList::null()
            labelList& uniqueMeshPointLabels,
            bool parallel
        ) const;


        //- Write cell connectivity for polyhedral cells
        static void writePolysConnectivity
        (
            ensightGeoFile& os,
            const polyMesh& mesh,
            const ensightCells& part,
            const labelList& pointToGlobal,
            bool parallel
        );

        //- Write cell connectivity for specified (non-poly) type
        static void writeShapeConnectivity
        (
            ensightGeoFile& os,
            const polyMesh& mesh,
            const ensightCells::elemType etype,
            const ensightCells& part,
            const labelList& pointToGlobal,
            bool parallel
        );


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("ensightCells");


    // Constructors

        //- Default construct, with part index 0
        ensightCells();

        //- Default construct, with description/partName
        explicit ensightCells(const string& description);


    //- Destructor
    virtual ~ensightCells() = default;


    // Member Functions

    // Access

        //- Manifold mesh cells detected? Globally consistent quantity.
        inline bool manifold() const noexcept;

        //- Processor-local size of all elements.
        using ensightPart::size;

        //- Processor-local size of the specified element type.
        inline label size(const elemType etype) const;

        //- Processor-local offset/size of element type.
        inline labelRange range(const elemType etype) const;

        //- The global size of all element types.
        //  This value is only meaningful after a reduce operation.
        label totalSize() const noexcept;

        //- Same as totalSize()
        label total() const noexcept { return totalSize(); }

        //- The global size of the specified element type.
        //  This value is only meaningful after a reduce operation.
        inline label total(const elemType etype) const;

        //- The global sizes for each element type.
        //  This value is only meaningful after a reduce operation.
        inline const FixedList<label, nTypes>& totals() const;

        //- Processor-local sizes per element type.
        FixedList<label, nTypes> sizes() const;

        //- Processor-local cell ids of all elements
        inline const labelList& cellIds() const;

        //- Processor-local cell ids of the specified element type
        inline const labelUList cellIds(const elemType etype) const;


    // Addressing

        //- Mesh point map.
        //  Map mesh point index to local (compact) point index
        Map<label> meshPointMap(const polyMesh& mesh) const;


    // Edit

        //- Classify cell types and set the element lists.
        void classify(const polyMesh& mesh);

        //- Classify cell types and set element lists,
        //- using a subgroup of cells (eg, from a cellZone etc).
        void classify(const polyMesh& mesh, const labelUList& cellIds);

        //- Classify cell types and set element lists,
        //- using a subgroup of cells
        void classify(const polyMesh& mesh, const bitSet& selection);


        //- Clear any demand-driven data
        void clearOut();

        //- Set addressable sizes to zero, free up addressing memory.
        void clear();

        //- Sum element counts across all processes.
        void reduce();

        //- Sort element lists numerically.
        void sort();


    // Advanced (use with caution)

        //- Increase cell ids by specified offset value
        //  Eg, to change zone local Ids to global Ids
        inline void incrCellIds(const label off);

        //- Decrease face ids by specified offset value
        //  Eg, to change global Ids to zone local Ids
        inline void decrCellIds(const label off);



    // Output

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<ensightCells> info() const noexcept { return *this; }


        //- Globally unique mesh points. Required when writing point fields.
        label uniqueMeshPoints
        (
            const polyMesh& mesh,
            labelList& uniqueMeshPointLabels,
            bool parallel
        ) const;


        //- Write information about the object as a dictionary,
        //- optionally write all element addresses
        virtual void writeDict(Ostream& os, const bool full=false) const;

        //- Write geometry, using a mesh reference (serial only)
        //- No beginGeometry() marker.
        virtual void write
        (
            ensightGeoFile& os,
            const polyMesh& mesh,
            bool parallel
        ) const;

        //- Write bounding box geometry.
        //- All parameters are only relevant on master
        //- No beginGeometry() marker.
        static void writeBox
        (
            ensightGeoFile& os,
            const boundBox& bb,
            const label partIndex = 0,
            const word& partName = "geometry-box"
        );
};


//- Print information via proxy
template<>
Ostream& operator<<(Ostream&, const InfoProxy<ensightCells>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightCellsI.H"

#endif

// ************************************************************************* //
