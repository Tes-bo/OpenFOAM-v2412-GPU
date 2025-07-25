/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::radiation::viewFactor

Group
    grpRadiationModels

Description
    View factor radiation model. The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactor.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_viewFactor_H
#define radiation_viewFactor_H

#include "radiationModel.H"
#include "singleCellFvMesh.H"
#include "scalarMatrices.H"
#include "globalIndex.H"
#include "scalarListIOList.H"
#include "volFields.H"
#include "IOmapDistribute.H"
#include "solarLoad.H"

#include "lduPrimitiveMesh.H"
#include "lduPrimitiveProcessorInterface.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class viewFactor Declaration
\*---------------------------------------------------------------------------*/

class viewFactor
:
    public radiationModel
{
private:

    // Private Member Functions

        //- No copy construct
        viewFactor(const viewFactor&) = delete;

        //- No copy assignment
        void operator=(const viewFactor&) = delete;


protected:

    // Static data

        //- Static name for view factor walls
        static const word viewFactorWalls;


    // Private data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<IOmapDistribute> map_;

        //- Coarse mesh
        autoPtr<singleCellFvMesh> coarseMesh_;

        //- Net radiative heat flux [W/m2]
        volScalarField qr_;

        //- Coarse radiative heat flux
        List<scalarField> qrBandI_;

        //- View factor matrix
        autoPtr<scalarSquareMatrix> Fmatrix_;

        //- Inverse of C matrix
        autoPtr<scalarSquareMatrix> CLU_;

        //- Visible global faces
        autoPtr<labelListIOList> globalFaceFaces_;

        //- Selected patches
        labelList selectedPatches_;

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;

        //- Constant emissivity
        bool constEmissivity_;

        //- Iterations Counter
        label iterCounter_;

        //- Pivot Indices for LU decomposition
        labelList pivotIndices_;

        //- Use Solar Load model
        bool useSolarLoad_;

        //- Solar load radiation model
        autoPtr<solarLoad> solarLoad_;

        //-Number of bands
        label nBands_;

        //- Primitive addressing for lduMatrix
        autoPtr<lduPrimitiveMesh> lduPtr_;

        //- Matrix formed from view factors
        autoPtr<lduMatrix> matrixPtr_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //- for internal cells
        FieldField<Field, scalar> internalCoeffs_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //- for boundary cells
        FieldField<Field, scalar> boundaryCoeffs_;

        //- Rays on local proc
        edgeHashSet rays_;

        //- Map local-ray to j-column for F
        List<label> mapRayToFmy_;

        //- Local view factors
        autoPtr<scalarListIOList> FmyProc_;

        //- Map from proc to interafce
        labelList procToInterface_;

        //- Use direct or iterative solver
        bool useDirect_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Insert view factors into main matrix
        void insertMatrixElements
        (
            const globalIndex& index,
            const label fromProci,
            const labelListList& globalFaceFaces,
            const scalarListList& viewFactors,
            scalarSquareMatrix& matrix
        );


public:

    //- Runtime type information
    TypeName("viewFactor");


    // Constructors

        //- Construct from components
        viewFactor(const volScalarField& T);

        //- Construct from components
        viewFactor(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~viewFactor() = default;


    // Member functions

        // Edit

            //- Solve system of equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


    // Access

        //- Const access to total radiative heat flux field
        inline const volScalarField& qr() const;

        //- Use useSolarLoad
        inline bool useSolarLoad() const;

        //- Number of bands
        virtual label nBands() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "viewFactorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
