/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::decompositionMethod

Description
    Abstract base class for domain decomposition

SourceFiles
    decompositionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_decompositionMethod_H
#define Foam_decompositionMethod_H

#include "polyMesh.H"
#include "CompactListList.H"
#include "decompositionConstraint.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class decompositionMethod Declaration
\*---------------------------------------------------------------------------*/

class decompositionMethod
{
    // Private Member Functions

        //- Check (and warn) about existence of old constraint syntax.
        //  The syntax changed MAY-2014.
        //
        //  \return true if this model name was found in the decompDict_
        //      but not previously added with the newer syntax.
        bool constraintCompat(const word& modelType) const;

        //- Set PtrList of constraints by reading decompDict_.
        void readConstraints();


protected:

        //- Selection type when handling the coefficients dictionary.
        //  To be used as a bit-mask for findCoeffsDict
        enum selectionType
        {
            DEFAULT = 0,   //!< Default request
            EXACT   = 1,   //!< No fallback to "coeffs" if main name not found
            MANDATORY = 2, //!< Fatal if dictionary could not be found
            NULL_DICT = 4, //!< On failure, return dictionary::null instead
                           //!< of the top-level enclosing dictionary.
        };


    // Protected Data

        //- Top-level decomposition dictionary (eg, decomposeParDict)
        const dictionary& decompDict_;

        //- Region-specific decomposition dictionary information
        const dictionary& decompRegionDict_;

        //- Number of domains for the decomposition
        label nDomains_;

        //- Optional constraints
        PtrList<decompositionConstraint> constraints_;


    // Protected Member Functions

        //- Locate coeffsName dictionary or the fallback "coeffs" dictionary
        //- within an enclosing dictionary.
        //
        //  \param select choose to include "coeffs" in the search, make
        //      failure a FatalError, return dictionary::null instead on
        //      failure.
        //
        //  \return the coefficients dictionary found. If nothing was found,
        //  return the enclosing dictionary or
        //  dictionary::null (depending on the select parameter).
        static const dictionary& findCoeffsDict
        (
            const dictionary& dict,
            const word& coeffsName,
            int select = selectionType::DEFAULT
        );


        //- Locate coeffsName dictionary or the fallback "coeffs" dictionary.
        //  Searches both the region-specific decomposition dictionary
        //  and the top-level decomposition dictionary.
        //
        //  \param select choose to include "coeffs" in the search, make
        //      failure a FatalError, return dictionary::null instead on
        //      failure.
        //
        //  \return the coefficients dictionary found. If nothing was found,
        //  return the top-level (non-region) dictionary or
        //  dictionary::null (depending on the select parameter).
        const dictionary& findCoeffsDict
        (
            const word& coeffsName,
            int select = selectionType::DEFAULT
        ) const;


    // Constructors

        //- Construct with specified number of domains,
        //- no coefficients or constraints
        explicit decompositionMethod(const label numDomains);


public:

    // Generated Methods

        //- No copy construct
        decompositionMethod(const decompositionMethod&) = delete;

        //- No copy assignment
        void operator=(const decompositionMethod&) = delete;


    //- Runtime type information
    TypeName("decompositionMethod");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionMethod,
            dictionary,
            (
                const dictionary& decompDict,
                const word& regionName
            ),
            (decompDict, regionName)
        );


    // Static Methods

        //- Return region-specific or top-level \c numberOfSubdomains entry.
        //  The region-specific version is found within the "regions"
        //  sub-dictionary.
        static label nDomains
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );

        //- Return an optional region-specific dictionary
        //- from "regions" sub-dictionary, or dictionary::null on failure
        static const dictionary& optionalRegionDict
        (
            const dictionary& decompDict,
            const word& regionName
        );


    // Selectors

        //- Return a reference to the selected decomposition method,
        //- optionally region-specific
        static autoPtr<decompositionMethod> New
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    // Constructors

        //- Construct given the decomposition dictionary,
        //- optionally region-specific
        explicit decompositionMethod
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~decompositionMethod() = default;


    // Member Functions

        //- Number of domains
        label nDomains() const noexcept
        {
            return nDomains_;
        }

        //- True if the method is purely geometric,
        //- often using cell centre points
        virtual bool geometric() const { return false; }

        //- Is method parallel aware?
        //  (i.e. does it synchronize domains across proc boundaries)
        virtual bool parallelAware() const = 0;

        // //- Is internally method parallel aware
        // virtual bool parallelNative() const  { return false; }


    // No topology (implemented by geometric decomposers)

        //- Return the wanted processor number for every coordinate,
        //- using uniform or specified point weights.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights = scalarField::null()
        ) const;


    // Topology provided by mesh

        //- Return for every coordinate the wanted processor number,
        //- using uniform or specified point weights.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights = scalarField::null()
        ) const = 0;

        //- Return for every coordinate the wanted processor number.
        //  Gets passed agglomeration map (from fine to coarse cells)
        //  and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively. Coarse cells are local to the processor
        //  (if in parallel). If you want to have coarse cells spanning
        //  processors use the globalCellCells instead.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& cellToRegion,
            const pointField& regionPoints,
            const scalarField& regionWeights = scalarField::null()
        ) const;


    // Topology provided explicitly

        //- Return for every coordinate the wanted processor number.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers
        //    (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const = 0;

        //- Return for every coordinate the wanted processor number.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers
        //    (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const = 0;


    // Other

        //- Helper: extract constraints:
        //  blockedface: existing faces where owner and neighbour on same
        //               proc
        //  explicitConnections: sets of boundary faces  ,,     ,,
        //  specifiedProcessorFaces: groups of faces with all cells on
        //  same processor.
        void setConstraints
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Helper: apply constraints to a decomposition.
        //  This gives constraints opportunity to modify decomposition in case
        //  the native decomposition method has not obeyed all constraints
        void applyConstraints
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& finalDecomp
        ) const;

        // Decompose a mesh with constraints:
        // - blockedFace : whether owner and neighbour should be on same
        //   processor
        // - specifiedProcessorFaces, specifiedProcessor : sets of faces
        //   that should go to same processor (as specified in
        //   specifiedProcessor, can be -1)
        // - explicitConnections : connections between baffle faces
        //   (blockedFace should be false on these). Owner and
        //   neighbour on same processor.
        // Set all to zero size to have unconstrained decomposition.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const scalarField& cellWeights,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections
        ) const;


        //- Decompose a mesh.
        //  Apply all constraints from decomposeParDict
        //  ('preserveFaceZones' etc). Calls either
        //  - no constraints, empty weights:
        //      decompose(mesh, cellCentres())
        //  - no constraints, set weights:
        //      decompose(mesh, cellCentres(), cellWeights)
        //  - valid constraints:
        //      decompose(mesh, cellToRegion, regionPoints, regionWeights)
        labelList decompose
        (
            const polyMesh& mesh,
            const scalarField& cWeights
        ) const;


    // Housekeeping

        //- Determine (local or global) cellCells from mesh agglomeration.
        //  Agglomeration is local to the processor.
        //  local  : connections are in local indices. Coupled across
        //           cyclics but not processor patches.
        //  global : connections are in global indices. Coupled across
        //            cyclics and processor patches.
        FOAM_DEPRECATED_STRICT(2023-11, "globalMeshData::calcCellCells()")
        static void calcCellCells
        (
            const polyMesh& mesh,
            const labelList& agglom,
            const label nLocalCoarse,
            const bool parallel,  //!< Use global cell ids in parallel
            CompactListList<label>& cellCells
        );

        //- Determine (local or global) cellCells and face weights
        //- from mesh agglomeration.
        //  Uses mag of faceArea as weights
        FOAM_DEPRECATED_STRICT(2023-11, "globalMeshData::calcCellCells()")
        static void calcCellCells
        (
            const polyMesh& mesh,
            const labelList& agglom,
            const label nLocalCoarse,
            const bool parallel,  //!< Use global cell ids in parallel
            CompactListList<label>& cellCells,
            CompactListList<scalar>& cellCellWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
