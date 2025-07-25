/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::solutionControl

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "fvMesh.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
:
    public regIOobject
{
public:

    //- Simple convenient storage of field residuals
    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


    // Static Member Functions

        //- Extract maximum residual for the solver performance entry,
        //- provided the corresponding volume field is available on the mesh.
        //
        //  \return initial residual as first member, the final residual
        //      as the second (or last) member
        static Pair<scalar> maxResidual
        (
            const fvMesh& fvmesh,
            const entry& dataDictEntry
        );


protected:

    // Protected data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- List of residual data per field
        List<fieldData> residualControl_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            bool frozenFlow_;


        // Evolution

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual bool read(const bool absTolOnly);

        //- Read controls from fvSolution dictionary
        virtual bool read();

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const word& fieldName,
            const bool useRegEx = true
        ) const;

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied() = 0;

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;

        //- Initial and final residual of the specified field-name,
        //- provided that the corresponding volume field is available
        //- on the fvMesh.
        //
        //  Populate residuals with initial residual as first member and
        //  the final residual as second (last) member.
        template<class Type>
        static bool maxTypeResidual
        (
            const fvMesh& fvmesh,
            const entry& solverPerfDictEntry,
            Pair<scalar>& residuals
        );

        //- Extract the maximum residual for the specified field
        //
        //  \return initial residual as first member, the final residual
        //      as second (last) member
        Pair<scalar> maxResidual(const entry& solverPerfDictEntry) const;

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        );

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;


private:

        //- No copy construct
        solutionControl(const solutionControl&) = delete;

        //- No copy assignment
        void operator=(const solutionControl&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");


    // Constructors

        //- Construct from mesh
        solutionControl(fvMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~solutionControl() = default;


    // Member Functions

        // Access

            //- Return the solution dictionary
            virtual const dictionary dict() const;

            //- Current corrector loop index
            inline label corr() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;

            //- Flag to indicate that the flow system of equations should not
            //  be evolved
            inline bool frozenFlow() const;


        // Evolution

            //- Main control loop
            virtual bool loop() = 0;

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
