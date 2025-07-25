/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022,2024 OpenCFD Ltd.
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
    Foam::oversetFvMeshBase

Description
    Support for overset functionality.

SourceFiles
    oversetFvMeshBase.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvMeshBase_H
#define oversetFvMeshBase_H

#include "fvMeshPrimitiveLduAddressing.H"
#include "lduInterfaceFieldPtrsList.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduPrimitiveProcessorInterface;
class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
                   Class oversetFvMeshBase Declaration
\*---------------------------------------------------------------------------*/

class oversetFvMeshBase
{
    // Private Member Functions

        //- No copy construct
        oversetFvMeshBase(const oversetFvMeshBase&) = delete;

        //- No copy assignment
        void operator=(const oversetFvMeshBase&) = delete;


protected:

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Select base addressing (false) or locally stored extended
        //  lduAddressing (true)
        mutable bool active_;

        //- Extended addressing (extended with local interpolation stencils)
        mutable autoPtr<fvMeshPrimitiveLduAddressing> lduPtr_;

        //- Added (processor)lduInterfaces for remote bits of stencil.
        mutable PtrList<const lduPrimitiveProcessorInterface>
            remoteStencilInterfaces_;

        //- Interfaces for above mesh. Contains both original and
        //- above added processorLduInterfaces
        mutable lduInterfacePtrsList allInterfaces_;

        //- Corresponding faces (in above lduPtr) to the stencil
        mutable labelListList stencilFaces_;

        //- Corresponding patches (in above lduPtr) to the stencil
        mutable labelListList stencilPatches_;

        //- From old to new face labels
        mutable labelList reverseFaceMap_;


   // Protected Member Functions

        //- Calculate the extended lduAddressing
        virtual bool updateAddressing() const;

        //- Debug: print matrix
        template<class Type>
        void write
        (
            Ostream&,
            const fvMatrix<Type>&,
            const lduAddressing&,
            const lduInterfacePtrsList&
        ) const;

        //- Freeze values at holes
        //template<class Type>
        //void freezeHoles(fvMatrix<Type>&) const;

         //- Scale coefficient depending on cell type
        template<class Type>
        void scaleConnection
        (
            Field<Type>& coeffs,
            const labelUList& types,
            const scalarList& factor,
            const bool setHoleCellValue,
            const label celli,
            const label facei
        ) const;

        //- Solve given dictionary with settings
        template<class Type>
        SolverPerformance<Type> solveOverset
        (
            fvMatrix<Type>& m,
            const dictionary&
        ) const;

        //- Debug: correct coupled bc
        template<class GeoField>
        static void correctCoupledBoundaryConditions(GeoField& fld);

        //- Average norm of valid neighbours
        scalar cellAverage
        (
            const labelList& types,
            const labelList& nbrTypes,
            const scalarField& norm,
            const scalarField& nbrNorm,
            const label celli,
            bitSet& isFront
        ) const;

        //- Debug: dump agglomeration
        void writeAgglomeration
        (
            const GAMGAgglomeration& agglom
        ) const;


public:

    //- Runtime type information
    TypeName("oversetFvMeshBase");


    // Constructors

        //- Construct from IOobject
        oversetFvMeshBase(const fvMesh& mesh, const bool doInit=true);


    //- Destructor
    virtual ~oversetFvMeshBase();


    // Member Functions

        // Extended addressing

            //- Return extended ldu addressing
            const fvMeshPrimitiveLduAddressing& primitiveLduAddr() const;

            //- Return ldu addressing. If active: is (extended)
            //  primitiveLduAddr
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set. If active:
            //  return additional remoteStencilInterfaces_
            virtual lduInterfacePtrsList interfaces() const;

            //- Return old to new face addressing
            const labelList& reverseFaceMap() const
            {
                return reverseFaceMap_;
            }

            //- Return true if using extended addressing
            bool active() const
            {
                return active_;
            }

            //- Enable/disable extended addressing
            void active(const bool f) const
            {
                active_ = f;

                if (active_)
                {
                    DebugInfo<< "Switching to extended addressing with nFaces:"
                        << primitiveLduAddr().lowerAddr().size()
                        << " nInterfaces:" << allInterfaces_.size()
                        << endl;
                }
                else
                {
                    DebugInfo<< "Switching to base addressing with nFaces:"
                        << mesh_.fvMesh::lduAddr().lowerAddr().size()
                        << " nInterfaces:" << mesh_.fvMesh::interfaces().size()
                        << endl;
                }
            }


        // Overset

            //- Manipulate the matrix to add the interpolation/set hole
            //  values
            template<class Type>
            void addInterpolation
            (
                fvMatrix<Type>& m,
                const scalarField& normalisation,
                const bool setHoleCellValue,
                const Type& holeCellValue
            ) const;


        //- Clear out local storage
        void clearOut();

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update fields when mesh is updated
        virtual bool interpolateFields();

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Debug: check halo swap is ok
        template<class GeoField>
        static void checkCoupledBC(const GeoField& fld);

        //- Correct boundary conditions of certain type (typeOnly = true)
        //- or explicitly not of the type (typeOnly = false)
        template<class GeoField, class PatchType>
        static void correctBoundaryConditions
        (
            typename GeoField::Boundary& bfld,
            const bool typeOnly
        );

        //- Determine normalisation for interpolation. This equals the
        //- original diagonal except stabilised for zero diagonals (possible
        //- in hole cells)
        template<class Type>
        tmp<scalarField> normalisation(const fvMatrix<Type>& m) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "oversetFvMeshBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
