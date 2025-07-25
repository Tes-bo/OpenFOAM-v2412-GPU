/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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
    Foam::noConstraint

Description
    Applies no constraints to the control points.
    Enforces the non-overlapping bounds, if present.

SourceFiles
    noConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef noConstraint_H
#define noConstraint_H

#include "morphingBoxConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noConstraint Declaration
\*---------------------------------------------------------------------------*/

class noConstraint
:
    public morphingBoxConstraint
{

protected:

    // Protected Member Functions

        //- Compute sensitivities wrt the design variables (chain rule)
        void virtual computeDVsSensitivities
        (
            scalarField& dvSens,
            const scalarField& cpSens
        );

        //- Update the bounds of the internal control points
        void updateInternalBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds,
            const NURBS3DVolume& boxI,
            const label passed
        );

        //- Update the bounds of the boundary control points
        void updateBoundaryBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds,
            const NURBS3DVolume& boxI,
            const label passed
        );


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        noConstraint(const noConstraint&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const noConstraint&) = delete;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noConstraint
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volumetricBSplinesDesignVariables& designVariables
        );


    //- Destructor
    virtual ~noConstraint() = default;


    // Member Functions

        //- Compute the active design variables based on the IDs of the
        //- active control point coordinates
        virtual labelList computeActiveDesignVariables
        (
            const labelList& activeCPCoors
        );

        //- Transform bounds from control points to design variables
        //  Does nothing in this case
        virtual void computeBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds
        );

        //- Update the bounds of the design variables
        //  Will update the bound  values if nonOverlappingCPs is active
        virtual void updateBounds
        (
            autoPtr<scalarField>& lowerBounds,
            autoPtr<scalarField>& upperBounds
        );

        //- Convert design variables to control points, stored in a scalarField
        virtual tmp<scalarField> designVariablesToControlPoints
        (
            const scalarField& designVariables
        );

        //- Return the design variables corresponding to the given control
        //- points
        virtual tmp<scalarField> controlPointsToDesignVariables
        (
            const scalarField& cps
        );

        //- Convert the correction of the design variables to the correction of
        //- the control points
        virtual tmp<scalarField> correctionCPs
        (
            const scalarField& correctionDVs
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
