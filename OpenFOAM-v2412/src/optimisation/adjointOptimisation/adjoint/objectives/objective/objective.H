/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::objective

Description
    Abstract base class for objective functions. No point in making this
    runTime selectable since its children will have different constructors.

SourceFiles
    objective.C

\*---------------------------------------------------------------------------*/

#ifndef objective_H
#define objective_H

#include "localIOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "boundaryFieldsFwd.H"
#include "solverControl.H"
#include "objectiveFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class objective Declaration
\*---------------------------------------------------------------------------*/

class objective
:
    public localIOdictionary
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;
        const word adjointSolverName_;
        const word primalSolverName_;
        const word objectiveName_;
        bool computeMeanFields_;
        bool nullified_;
        bool normalize_;
        bool shouldWrite_;

        //- Objective function value and weight
        scalar J_;

        //- Average objective value
        scalar JMean_;

        //- Objective weight
        scalar weight_;

        //- Whether the objective is computed or not
        //  Some objective (e.g. geometric ones) might not change from one
        //  iteration of the primal solver to the next and might be expensive
        //  to evaluate. This can be used to compute them only once per
        //  optimisation cycle
        bool computed_;

        //- Normalization factor
        autoPtr<scalar> normFactor_;

        //- Target value, in case the objective is used as a constraint
        //  Should be used with caution and with updateMethods
        //  than get affected by the target value, without
        //  requiring a sqr (e.g. SQP, MMA)
        autoPtr<scalar> target_;

        //- Target on the left hand-side of a double inequality,
        //- for double sided constraints
        autoPtr<scalar> targetLeft_;

        //- Objective integration start and end times (for unsteady flows)
        autoPtr<scalar> integrationStartTimePtr_;
        autoPtr<scalar> integrationEndTimePtr_;

        //- List of adjoint fields for which this objective will contribute
        //- sources to their equations.
        //  Only for volume-based objectives for the moment
        wordList fieldNames_;

        //- Contribution to field sensitivity derivatives
        //  Topology optimisation or other variants with
        //  as many design variables as the mesh cells
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<volScalarField> dJdbPtr_;

        //- Contribution to sensitivity derivatives with a
        //- random number of design variables
        //- (neither surface, nor volume based)
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<scalarField> dJdbFieldPtr_;

        // Contribution to surface sensitivity derivatives
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Term from material derivative
        autoPtr<boundaryVectorField> bdJdbPtr_;

        //- Term multiplying delta(n dS)/delta b
        autoPtr<boundaryVectorField> bdSdbMultPtr_;

        //- Term multiplying delta(n)/delta b
        autoPtr<boundaryVectorField> bdndbMultPtr_;

        //- Term multiplying delta(x)/delta b at the boundary
        autoPtr<boundaryVectorField> bdxdbMultPtr_;

        //- Term multiplying delta(x)/delta b at the boundary
        //- for objectives that directly depend on x, e.g. moment
        //- Needed in both FI and SI computations
        autoPtr<boundaryVectorField> bdxdbDirectMultPtr_;

        //- Contribution located in specific parts of a patch.
        //- Mainly intended for patch boundary edges contributions, e.g.
        //- NURBS surfaces G1 continuity
        autoPtr<vectorField3> bEdgeContribution_;

        // Contribution to volume-based sensitivities from volume-based
        // objective functions
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Multiplier of d(Volume)/db
        autoPtr<volScalarField> divDxDbMultPtr_;

        //- Emerging from volume objectives that include spatial derivatives
        autoPtr<volTensorField> gradDxDbMultPtr_;

        //- Output file variables
        fileName objFunctionFolder_;

        //- File to keep the objective values after the end of the primal solver
        mutable autoPtr<OFstream> objFunctionFilePtr_;

        //- File to keep the objective values at each iteration of the primal
        //- solver
        mutable autoPtr<OFstream> instantValueFilePtr_;

        //- File to keep the average objective values after the end of the
        //- primal solver
        mutable autoPtr<OFstream> meanValueFilePtr_;

        //- Default width of entries when writing in the objective files
        unsigned int width_;


    // Protected Member Functions

        //- Set the output file ptr
        void setObjectiveFilePtr() const;

        //- Set the output file ptr for the instantaneous value
        void setInstantValueFilePtr() const;

        //- Set the output file ptr for the mean value
        void setMeanValueFilePtr() const;


private:

    // Private Member Functions

        //- No copy construct
        objective(const objective&) = delete;

        //- No copy assignment
        void operator=(const objective&) = delete;

        //- Make objective Function Folder
        void makeFolder();


public:

    //- Runtime type information
    TypeName("objective");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            objective,
            objective,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& adjointSolverName,
                const word& primalSolverName
            ),
            (mesh, dict, adjointSolverName, primalSolverName)
        );


    // Constructors

        //- Construct from components
        objective
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objective> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& objectiveType,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objective() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        //- Return the instantaneous objective function value
        virtual scalar J() = 0;

        //- Return the objective function of the optimisation cycle.
        //  This corresponds to the mean value, if it exists, or the
        //  instantaneous value otherwise
        virtual scalar JCycle(bool negate = false) const;

        //- Accumulate contribution for the mean objective value
        //  For steady-state runs
        virtual void accumulateJMean(solverControl& solverControl);
        

        //- Accumulate contribution for the mean objective value
        //  For unsteady runs
        virtual void accumulateJMean();

        //- Return the objective function weight
        virtual scalar weight() const;

        //- Is the objective normalized
        virtual bool normalize() const;

        //- Normalize all fields allocated by the objective
        virtual void doNormalization();

        //- Check whether this is an objective integration time
        virtual bool isWithinIntegrationTime() const;

        //- Increment integration times
        virtual void incrementIntegrationTimes(const scalar timeSpan);

        //- Contribution to field sensitivities
        inline const volScalarField& dJdb() const;

        //- Contribution to sensitivities with a random number of designVars
        inline const scalarField& dJdbField() const;

        //- Contribution to surface sensitivities for a specific patch
        inline const fvPatchVectorField& boundarydJdb(const label) const;

        //- Multiplier of delta(n dS)/delta b
        inline const fvPatchVectorField& dSdbMultiplier(const label) const;

        //- Multiplier of delta(n dS)/delta b
        inline const fvPatchVectorField& dndbMultiplier(const label) const;

        //- Multiplier of delta(x)/delta b
        inline const fvPatchVectorField& dxdbMultiplier(const label) const;

        //- Multiplier of delta(x)/delta b
        inline const fvPatchVectorField& dxdbDirectMultiplier
        (
            const label
        ) const;

        //- Multiplier located at patch boundary edges
        inline const vectorField& boundaryEdgeMultiplier
        (
            const label patchI,
            const label edgeI
        ) const;

        //- Contribution to surface sensitivities for all patches
        inline const boundaryVectorField& boundarydJdb() const;

        //- Multiplier of delta(n dS)/delta b for all patches
        inline const boundaryVectorField& dSdbMultiplier() const;

        //- Multiplier of delta(n dS)/delta b for all patches
        inline const boundaryVectorField& dndbMultiplier() const;

        //- Multiplier of delta(x)/delta b for all patches
        inline const boundaryVectorField& dxdbMultiplier() const;

        //- Multiplier of delta(x)/delta b for all patches
        inline const boundaryVectorField& dxdbDirectMultiplier() const;

        //- Multiplier located at patch boundary edges
        inline const vectorField3& boundaryEdgeMultiplier() const;

        //- Multiplier of grad( delta(x)/delta b) for volume-based sensitivities
        inline const volScalarField& divDxDbMultiplier() const;

        //- Multiplier of grad( delta(x)/delta b) for volume-based sensitivities
        inline const volTensorField& gradDxDbMultiplier() const;

        //- Update objective function derivatives
        virtual void update() = 0;

        //- Nullify adjoint contributions
        virtual void nullify();

        //- Update normalization factors, for objectives in
        //- which the factor is not known a priori
        virtual void updateNormalizationFactor();


        virtual void update_dJdb()
        {}

        virtual void update_dJdbField()
        {}

        //- Update objective function derivative term
        virtual void update_boundarydJdb()
        {}

        //- Update d (normal dS) / db multiplier. Surface-based sensitivity term
        virtual void update_dSdbMultiplier()
        {}

        //- Update d (normal) / db multiplier. Surface-based sensitivity term
        virtual void update_dndbMultiplier()
        {}

        //- Update d (x) / db multiplier. Surface-based sensitivity term
        virtual void update_dxdbMultiplier()
        {}

        //- Update d (x) / db multiplier. Surface and volume-based sensitivity
        //- term
        virtual void update_dxdbDirectMultiplier()
        {}

        //- Update boundary edge contributions
        virtual void update_boundaryEdgeContribution()
        {}

        //- Update div( dx/db multiplier). Volume-based sensitivity term
        virtual void update_divDxDbMultiplier()
        {}

        //- Update grad( dx/db multiplier). Volume-based sensitivity term
        virtual void update_gradDxDbMultiplier()
        {}


        //- Manipulate fvVectorMatrix through the objective
        virtual void addSource(fvVectorMatrix& matrix)
        {}

        //- Manipulate fvVectorMatrix through the objective
        virtual void addSource(fvScalarMatrix& matrix)
        {}

        //- Write objective function history
        virtual bool write(const bool valid = true) const;

        //- Write objective function history at each primal solver iteration
        virtual void writeInstantaneousValue() const;

        //- Append a blank line after the end of optimisation cycle to the
        //- file holding the instantaneous objective values to facilitate
        //- plotting
        virtual void writeInstantaneousSeparator() const;

        //- Write mean objective function history
        virtual void writeMeanValue() const;

        //- Write averaged objective for continuation
        virtual bool writeData(Ostream& os) const;

        // Helper write functions

            //- Write any information that needs to go the header of the file
            //  (e.g. targets, directions, etc)
            virtual void addHeaderInfo() const;

            //- Write headers for additional columns
            virtual void addHeaderColumns() const;

            //- Write information to additional columns
            virtual void addColumnValues() const;

        //- Return the objective name
        inline const word& objectiveName() const;

        //- Should the objective be written to file upon calling write()?
        inline bool shouldWrite() const;

        //- Should the objective be written to file upon calling write()?
        inline void setWrite(const bool shouldWrite);

        // Inline functions for checking whether pointers are set or not
        inline bool hasdJdb() const noexcept;
        inline bool hasdJdbField() const noexcept;
        inline bool hasBoundarydJdb() const noexcept;
        inline bool hasdSdbMult() const noexcept;
        inline bool hasdndbMult() const noexcept;
        inline bool hasdxdbMult() const noexcept;
        inline bool hasdxdbDirectMult() const noexcept;
        inline bool hasBoundaryEdgeContribution() const noexcept;
        inline bool hasDivDxDbMult() const noexcept;
        inline bool hasGradDxDbMult() const noexcept;

        // Inline functions for checking whether integration times are set
        inline bool hasIntegrationStartTime() const noexcept;
        inline bool hasIntegrationEndTime() const noexcept;

        //- Set the computed status of the objective
        inline void setComputed(const bool isComputed) noexcept;

        //- Return objective dictionary
        const dictionary& dict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectiveI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
