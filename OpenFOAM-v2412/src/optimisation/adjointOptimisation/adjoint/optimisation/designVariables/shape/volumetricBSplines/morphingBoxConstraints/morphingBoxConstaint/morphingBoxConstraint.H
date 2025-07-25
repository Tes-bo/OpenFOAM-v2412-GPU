/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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
    Foam::morphingBoxConstraint

Description
    Abstract base class for defining constraints for the control points of
    volumetric B-Splines morphing boxes

SourceFiles
    morphingBoxConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef morphingBoxConstraint_H
#define morphingBoxConstraint_H

#include "volBSplinesBase.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class volumetricBSplinesDesignVariables;

/*---------------------------------------------------------------------------*\
                    Class morphingBoxConstraint Declaration
\*---------------------------------------------------------------------------*/

class morphingBoxConstraint
{
protected:

    // Protected data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Volumetric B-Splines variables dict
        const dictionary dict_;

        //- Reference to underlaying volumetric B-Splines morpher
        volumetricBSplinesDesignVariables& designVariables_;

        //- Easy access to the volBSplinesBase resting in the designVariables_
        volBSplinesBase& volBSplinesBase_;

        //- Initial CPs stored in scalarField
        scalarField initialCPs_;

        //- Initialise the design variables
        bool initialiseVars_;

        //- Folder holding the twist sensitivities
        fileName derivativesFolder_;


    // Protected Member Functions

        //- Compute sensitivities wrt the design variables (chain rule)
        void virtual computeDVsSensitivities
        (
            scalarField& dvSens,
            const scalarField& cpSens
        ) = 0;

        //- Write sensitivities w.r.t. the design variables
        void virtual writeDVSensitivities
        (
            const scalarField& sens,
            const word& name
        );


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        morphingBoxConstraint(const morphingBoxConstraint&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const morphingBoxConstraint&) = delete;


public:

    //- Runtime type information
    TypeName("morphingBoxConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            morphingBoxConstraint,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                volumetricBSplinesDesignVariables& designVariables
            ),
            (mesh, dict, designVariables)
        );


    // Constructors

        //- Construct from components
        morphingBoxConstraint
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volumetricBSplinesDesignVariables& designVariables
        );


    // Selectors

        //- Construct and return the selected morphingBoxConstraint
        static autoPtr<morphingBoxConstraint> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volumetricBSplinesDesignVariables& designVariables
        );


    //- Destructor
    virtual ~morphingBoxConstraint() = default;


    // Member Functions

        //- Compute the active design variables based on the IDs of the
        //- active control point coordinates
        virtual labelList computeActiveDesignVariables
        (
            const labelList& activeCPCoors
        ) = 0;

        //- Transform bounds from control points to design variables
        //  WIP
        virtual void computeBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds
        );

        //- Update the bounds of the design variables
        //  WIP
        virtual void updateBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds
        );

        //- Convert design variables to control points, stored in a scalarField
        virtual tmp<scalarField> designVariablesToControlPoints
        (
            const scalarField& designVariables
        ) = 0;

        //- Return the design variables corresponding to the given control
        //- points
        virtual tmp<scalarField> controlPointsToDesignVariables
        (
            const scalarField& cps
        ) = 0;

        //- Convert the correction of the design variables to the correction of
        //- the control points
        virtual tmp<scalarField> correctionCPs
        (
            const scalarField& correction
        ) = 0;

        //- Chain rule from control points to design variables
        virtual tmp<scalarField> postProcessSens
        (
            const scalarField& controlPointSens,
            const word& adjointSolverName
        );

        //- Compute eta if not set in the first step
        virtual scalar computeEta
        (
            scalarField& correction,
            const scalar maxInitChange
        );

        //- Initialise the design variables?
        inline bool initialiseVars() const
        {
            return initialiseVars_;
        }

        //- Append useful information to the design variables IOdictionary
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
