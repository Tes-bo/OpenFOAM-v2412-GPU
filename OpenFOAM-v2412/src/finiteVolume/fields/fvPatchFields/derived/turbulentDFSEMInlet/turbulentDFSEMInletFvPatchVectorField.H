/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::turbulentDFSEMInletFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    The \c turbulentDFSEMInlet is a synthesised-eddy based velocity inlet
    boundary condition to generate synthetic turbulence-alike time-series
    from a given set of turbulence statistics for LES and hybrid RANS-LES
    computations.

    Reference:
    \verbatim
        Standard model (tag:PCR):
            Poletto, R., Craft, T., & Revell, A. (2013).
            A new divergence free synthetic eddy method for
            the reproduction of inlet flow conditions for LES.
            Flow, turbulence and combustion, 91(3), 519-539.
            DOI:10.1007/s10494-013-9488-2

        Expression for the average length scale (tag:SST):
            Shur, M., Strelets, M., Travin, A.,
            Probst, A., Probst, S., Schwamborn, D., ... & Revell, A. (2018).
            Improved embedded approaches.
            In: Mockett C., Haase W., Schwamborn D. (eds)
            Go4Hybrid: Grey area mitigation for hybrid RANS-LES methods.
            Notes on Numerical Fluid Mechanics and Multidisciplinary Design.
            p. 51-87. Springer, Cham.
            DOI:10.1007/978-3-319-52995-0_3
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            turbulentDFSEMInlet;
        delta           <scalar>;
        R               <PatchFunction1>;
        U               <PatchFunction1>;
        L               <PatchFunction1>;

            // e.g.
            // R        uniform (<Rxx> <Rxy> <Rxz> <Ryy> <Ryz> <Rzz>);
            // U        uniform (<Ux> <Uy> <Uz>);
            // L        uniform <L>;

        // Optional entries
        d               <scalar>;
        nCellPerEddy    <label>;
        kappa           <scalar>;
        Uref            <scalar>;
        Lref            <scalar>;
        scale           <scalar>;
        m               <scalar>;
        writeEddies     <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type       | Reqd | Deflt
      type     | Type name: turbulentDFSEMInlet      | word       | yes  | -
      delta    | Characteristic length scale         | scalar     | yes  | -
      R        | Reynolds-stress tensor field        <!--
               -->                 | PatchFunction1\<symmTensor\> | yes  | -
      U        | Mean velocity field                 <!--
               -->                 | PatchFunction1<vector>       | yes  | -
      L        | Integral-length scale field         <!--
               -->                 | PatchFunction1<scalar>       | yes  | -
      d        | Ratio of sum of eddy volumes to eddy box volume  <!--
               --> i.e. eddy density (fill fraction) | scalar     | no   | 1.0
      nCellPerEddy | Minimum eddy length in units of number of cells     <!--
               -->                                   | label      | no   | 5
      kappa    | von Karman constant                 | scalar     | no   | 0.41
      Uref     | Normalisation factor for Reynolds-stress         <!--
               --> tensor and mean velocity          | scalar     | no   | 1.0
      Lref     | Normalisation factor for integral-length scale   <!--
               -->                                   | scalar     | no   | 1.0
      scale    | Heuristic scaling factor being applied           <!--
               --> on the normalisation factor C1    | scalar     | no   | 1.0
      m        | The power of V defined in C1        | scalar     | no   | 0.5
      writeEddies  | Flag to write eddies as OBJ file  | bool     | no   | false
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink
      - \link MappedFile.H \endlink

Note
  - The \c delta value typically represents the characteristic scale of flow
    or flow domain, e.g. a channel half height for plane channel flows.
  - \c nCellPerEddy and \c scale entries are heuristic entries
    which do not exist in the standard method, yet are necessary
    to reproduce the results published in the original journal paper.
  - In the original journal paper, \c C1 in Eq. 11 is not dimensionless.
    It is not clear whether this dimensionality issue was intentional.
    To alleviate this matter, users can alter the input entry \c m, which is
    the power of the eddy-box volume defined in the \c C1, to obtain a
    dimensionless \c C1 coefficient. The default value of \c m is 0.5,
    which is the value stated in the original journal paper,
    and \c m=1/3 leads to a dimensionless \c C1.

SourceFiles
    turbulentDFSEMInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFSEMInletFvPatchVectorField_H
#define turbulentDFSEMInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "eddy.H"
#include "pointIndexHit.H"
#include "PatchFunction1.H"
#include "labelledTri.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class turbulentDFSEMInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDFSEMInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Maximum number of attempts when seeding eddies
        static label seedIterMax_;

        //- Mean velocity field
        autoPtr<PatchFunction1<vector>> U_;

        //- Reynolds stress tensor field
        autoPtr<PatchFunction1<symmTensor>> R_;

        //- Integral-length scale field
        autoPtr<PatchFunction1<scalar>> L_;

        //- Characteristic length scale
        const scalar delta_;

        //- Ratio of sum of eddy volumes to eddy box volume, i.e. eddy density
        const scalar d_;

        //- von Karman constant
        const scalar kappa_;

        //- Normalisation factor for Reynolds-stress tensor and mean velocity
        const scalar Uref_;

        //- Normalisation factor for integral-length scale
        const scalar Lref_;

        //- Heuristic scaling factor being applied on the normalisation factor
        const scalar scale_;

        //- The power of V defined in C1
        const scalar m_;

        //- Minimum eddy length in units of number of cells
        const label nCellPerEddy_;


        // Patch information

            //- Patch area - total across all processors
            scalar patchArea_;

            //- The polyPatch faces as triangles, the index of each corresponds
            //- to the undecomposed patch face index.
            List<labelledTri> triFace_;

            //- Cumulative triangle area per triangle face (processor-local)
            scalarList triCumulativeMagSf_;

            //- Cumulative area fractions per processor
            scalarList sumTriMagSf_;

            //- Patch normal into the domain
            vector patchNormal_;

            //- Patch bounds (local processor)
            boundBox patchBounds_;


        // Eddy information

            //- List of eddies
            List<eddy> eddies_;

            //- Eddy box volume
            scalar v0_;

            //- Random number generator
            Random rndGen_;

            //- Integral-length scale per patch face
            scalarField sigmax_;

            //- Maximum integral-length scale (across all processors)
            scalar maxSigmaX_;

            //- Global number of eddies
            label nEddy_;

            //- Current time index (used for updating)
            label curTimeIndex_;

            //- Single processor contains all eddies (flag)
            bool singleProc_;

            //- Flag to write the eddies to file
            bool writeEddies_;


    // Private Member Functions

        //- Write Lumley coefficients to file
        void writeLumleyCoeffs() const;

        //- Write eddy info in OBJ format
        void writeEddyOBJ() const;

        //- Initialise info for patch point search
        void initialisePatch();

        //- Initialise the eddy box
        void initialiseEddyBox();

        //- Set a new eddy position
        pointIndexHit setNewPosition(const bool global);

        //- Initialise eddies
        void initialiseEddies();

        //- Convect the eddies with the bulk velocity
        void convectEddies(const vector& UBulk, const scalar deltaT);

        //- Return velocity fluctuation vector at a given point
        vector uPrimeEddy(const List<eddy>& eddies, const point& globalX) const;

        //- Return eddies from remote processors that interact with local
        //- processor
        void calcOverlappingProcEddies
        (
            List<List<eddy>>& overlappingEddies
        ) const;


public:

   //- Runtime type information
   TypeName("turbulentDFSEMInlet");


   // Constructors

        //- Construct from patch and internal field
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentDFSEMInletFvPatchVectorField
        //- onto a new patch
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~turbulentDFSEMInletFvPatchVectorField() = default;


    // Member Functions

        //- Check if input Reynolds stresses are valid
        //  Realizability conditions (tag:S):
        //      Schumann, U. (1977).
        //      Realizability of Reynolds‐stress turbulence models.
        //      The Physics of Fluids, 20(5), 721-725.
        //      DOI:10.1063/1.861942
        static void checkStresses(const symmTensorField& R);

        //- Check if input Reynolds stresses are valid
        static void checkStresses(const scalarField& R);


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField& ptf,
                const labelList& addr
            );


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // IO

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
