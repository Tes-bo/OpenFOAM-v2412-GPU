/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::faPatch

Description
    Finite area patch class.  Used for 2-D non-Euclidian finite area method.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faPatch.C
    faPatchNew.C
    faPatchTemplates.C
    faPatchFaMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faPatch_H
#define Foam_faPatch_H

#include "patchIdentifier.H"
#include "labelList.H"
#include "pointField.H"
#include "typeInfo.H"
#include "PtrList.H"
#include "faPatchFieldsFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faBoundaryMesh;
class faPatch;

//- Store lists of faPatch as a PtrList
typedef PtrList<faPatch> faPatchList;

Ostream& operator<<(Ostream&, const faPatch&);

/*---------------------------------------------------------------------------*\
                           Class faPatch Declaration
\*---------------------------------------------------------------------------*/

class faPatch
:
    public patchIdentifier,
    public labelList
{
    // Private Data

        //- Neighbour polyPatch index
        const label nbrPolyPatchId_;

        //- Reference to boundary mesh
        const faBoundaryMesh& boundaryMesh_;

        //- Demand-driven: edge-face addressing
        mutable std::unique_ptr<labelList::subList> edgeFacesPtr_;

        //- Demand-driven: local points labels
        mutable std::unique_ptr<labelList> pointLabelsPtr_;

        //- Demand-driven: point-edge addressing
        mutable std::unique_ptr<labelListList> pointEdgesPtr_;


    // Private Member Functions

        //- No copy construct
        faPatch(const faPatch&) = delete;

        //- No copy assignment
        void operator=(const faPatch&) = delete;

        //- Clear out topological patch data
        void clearOut();


protected:

    // Protected Member Functions

        //- The faPatch geometry initialisation is called by faBoundaryMesh
        friend class faBoundaryMesh;

        //- Calculate patch point labels
        void calcPointLabels() const;

        //- Calculate patch point-edge addressing
        void calcPointEdges() const;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patch after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&)
        {}


public:

    //- The boundary type associated with the patch
    typedef faBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faPatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const faBoundaryMesh& bm,
                const word& patchType
            ),
            (name, dict, index, bm, patchType)
        );


    // Constructors

        //- Construct from components
        faPatch
        (
            const word& name,
            const labelUList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchi,
            const word& patchType
        );

        //- Construct from dictionary
        faPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm,
            const word& patchType
        );

        //- Copy construct, resetting the boundary mesh
        faPatch(const faPatch& p, const faBoundaryMesh& bm);

        //- Copy construct, resetting boundary mesh and addressing
        faPatch
        (
            const faPatch& p,
            const faBoundaryMesh& bm,
            const label index,
            const labelUList& edgeLabels,
            const label nbrPolyPatchi
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<faPatch> clone(const faBoundaryMesh& bm) const
        {
            return autoPtr<faPatch>::New(*this, bm);
        }

        //- Construct and return a clone, resetting the edge list
        //- and boundary mesh
        virtual autoPtr<faPatch> clone
        (
            const faBoundaryMesh& bm,
            const labelUList& edgeLabels,
            const label index,
            const label nbrPolyPatchi
        ) const
        {
            return autoPtr<faPatch>::New
            (
                *this,
                bm,
                index,
                edgeLabels,
                nbrPolyPatchi
            );
        }


    // Selectors

        //- Return pointer to a new patch created on freestore from dictionary
        static autoPtr<faPatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );

        //- Return pointer to a new patch created on freestore from dictionary
        static autoPtr<faPatch> New
        (
            const word& patchType,
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );


    //- Destructor
    virtual ~faPatch();


    // Static Member Functions

        //- Return true if the given type is a constraint type
        static bool constraintType(const word& patchType);

        //- Return a list of all the constraint patch types
        static wordList constraintTypes();


    // Member Functions

        //- Return the list of edges
        const labelList& edgeLabels() const noexcept
        {
            return static_cast<const labelList&>(*this);
        }

        //- Number of patch points
        label nPoints() const
        {
            return pointLabels().size();
        }

        //- Number of edge labels (boundary edges) addressed by this patch
        label nEdges() const noexcept
        {
            return labelList::size();
        }

        //- The neighbour polyPatch index
        label ngbPolyPatchIndex() const noexcept
        {
            return nbrPolyPatchId_;
        }

        //- Return boundaryMesh reference
        const faBoundaryMesh& boundaryMesh() const noexcept;

        //- Return true if this patch is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- The offset where this patch starts in the boundary edge list.
        //  The value is the same as patch.start() - mesh.nInternalEdges()
        label offset() const;

        //- Patch start in edge list
        label start() const;

        //- Patch size is the number of edge labels, but can be overloaded
        virtual label size() const
        {
            return labelList::size();
        }

        //- Return label of edge in patch from global edge label
        label whichEdge(const label edgei) const
        {
            return edgei - start();
        }

        //- This patch slice from the complete list of values, which has
        //- size mesh::nEdges(), using the virtual patch size.
        template<class T>
        const typename List<T>::subList patchSlice(const List<T>& values) const
        {
            return typename List<T>::subList(values, size(), start());
        }

        //- This patch slice from the list of boundary values, which has
        //- size mesh::nBoundaryEdges(), using the virtual patch size.
        template<class T>
        const typename List<T>::subList
        boundarySlice(const List<T>& values) const
        {
            return typename List<T>::subList
            (
                values,
                size(),
                offset()
            );
        }


        //- Write
        virtual void write(Ostream&) const;


    // Topological information

        //- List of proc/face for the boundary edge neighbours
        //- in locally reordered edge numbering.
        List<labelPair> boundaryConnections() const;

        //- Boundary edge neighbour processors
        //- (does not include own proc)
        labelList boundaryProcs() const;

        //- List of proc/size for the boundary edge neighbour processors
        //- (does not include own proc)
        List<labelPair> boundaryProcSizes() const;


        // Access functions for geometrical data

            //- Return patch point labels
            const labelList& pointLabels() const;

            //- Return patch point-edge addressing
            const labelListList& pointEdges() const;

            //- Return normals of neighbour polyPatch faces
            //  Same as faMesh::haloFaceNormals()
            tmp<vectorField> ngbPolyPatchFaceNormals() const;

            //- Return normals of neighbour polyPatch joined points
            tmp<vectorField> ngbPolyPatchPointNormals() const;

            //- Return edge-face addressing
            const labelUList& edgeFaces() const;

            //- Return edge centres
            const vectorField& edgeCentres() const;

            //- Return edge length vectors, like the faMesh::Le() method
            const vectorField& edgeLengths() const;

            //- Return edge length magnitudes, like the faMesh::magLe() method
            const scalarField& magEdgeLengths() const;

            //- Return edge unit normals, like the faMesh::unitLe() method
            tmp<vectorField> edgeNormals() const;

            //- Return neighbour face centres
            tmp<vectorField> edgeFaceCentres() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorField> delta() const;


        // Access functions for demand driven data

            //- Make patch weighting factors
            virtual void makeWeights(scalarField&) const;

            //- Return patch weighting factors
            const scalarField& weights() const;

            //- Make patch geodesic distance between P and N
            virtual void makeLPN(scalarField&) const;

            //- Return patch geodesic distance between P and N
            const scalarField& lPN() const;

            //- Make patch edge - neighbour face distances
            virtual void makeDeltaCoeffs(scalarField&) const;

            //- Return patch edge - neighbour face distances
            const scalarField& deltaCoeffs() const;

            //- Make non-orthogonality correction vectors
            virtual void makeCorrectionVectors(vectorField&) const;

            //- Return patch skew-correction vectors
            const vectorField& skewCorrectionVectors() const;


        // Topological changes

            //- Reset the list of edges (use with caution)
            void resetEdges(const labelUList& newEdges);

            //- Reset the list of edges (use with caution)
            void resetEdges(labelList&& newEdges);


        // Evaluation

            //- Extract internal field next to patch using specified addressing
            // \param internalData The internal field to extract from
            // \param addressing Addressing from patch into internal field
            // \param [out] pfld The extracted patch field.
            //     It is always resized according to the patch size(),
            //     which can be smaller than the addressing size
            template<class Type>
            inline void patchInternalField
            (
                const UList<Type>& internalData,
                const labelUList& addressing,
                Field<Type>& pfld
            ) const;

            //- Extract internal field next to patch as patch field
            //- using edgeFaces() mapping.
            // \param internalData The internal field to extract from
            // \param [out] pfld The extracted patch field.
            //     It is always resized according to the patch size(),
            //     which can be smaller than the edgeFaces() size
            template<class Type>
            void patchInternalField
            (
                const UList<Type>& internalData,
                Field<Type>& pfld
            ) const;

            //- Return given internal field next to patch as patch field
            //- using edgeFaces() mapping.
            // \param internalData The internal field to extract from
            template<class Type>
            tmp<Field<Type>> patchInternalField
            (
                const UList<Type>& internalData
            ) const;

            //- Return the patch field of the GeometricField
            //- corresponding to this patch.
            template<class GeometricField, class AnyType = bool>
            const typename GeometricField::Patch& patchField
            (
                const GeometricField& gf
            ) const;

            //- Lookup the named field from the local registry and
            //- return the patch field corresponding to this patch.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            template<class GeometricField, class AnyType = bool>
            const typename GeometricField::Patch& lookupPatchField
            (
                const word& name,
                const GeometricField* = nullptr,
                const AnyType* = nullptr
            ) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const faPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
