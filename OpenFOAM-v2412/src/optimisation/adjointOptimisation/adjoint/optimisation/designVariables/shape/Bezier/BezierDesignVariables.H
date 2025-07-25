/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::BezierDesignVariables

Description
    Bezier design variables for shape optimisation

SourceFiles
    BezierDesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef BezierDesignVariables_H
#define BezierDesignVariables_H

#include "shapeDesignVariables.H"
#include "Bezier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BezierDesignVariables Decleration
\*---------------------------------------------------------------------------*/

class BezierDesignVariables
:
    public shapeDesignVariables
{
protected:

    // Protected Data Members

        //- The Bezier control points and auxiliary functions
        Bezier bezier_;

        //- Boundary movement due to the change in Bezier control points
        pointVectorField dx_;


    // Protected Member Functions

        //- Read bounds for design variables, if present
        void readBounds
        (
            autoPtr<scalar> lowerBoundPtr = nullptr,
            autoPtr<scalar> upperBoundPtr = nullptr
        );

        //- Set uniform bounds for all control points
        void setBounds(autoPtr<scalarField>& bounds, const vector& cpBounds);

        //- Transform the correction of design variables to control points'
        //- movement
        tmp<vectorField> computeBoundaryDisplacement
        (
            const scalarField& correction
        );

        //- Decompose varID to cpID and direction
        void decomposeVarID(label& cpI, label& dir, const label varID) const;


private:

    // Private Member Functions

        //- No copy construct
        BezierDesignVariables(const BezierDesignVariables&) = delete;

        //- No copy assignment
        void operator=(const BezierDesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("Bezier");


    // Constructors

        //- Construct from components
        BezierDesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BezierDesignVariables() = default;


    // Member Functions

        //- Update design variables based on a given correction
        virtual void update(scalarField& correction);

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction);

        //- Whether to use global sum when computing matrix-vector products
        //  in update methods
        virtual bool globalSum() const;


        // Fields related to sensitivity computations

            //- Get dxdb for given design variable and patch
            virtual tmp<vectorField> dxdbFace
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dndb for given design variable and patch
            virtual tmp<vectorField> dndb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dSdb for given design variable and patch
            virtual tmp<vectorField> dSdb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dCdb for given design variable.
            //  Used for FI-based sensitivities
            virtual tmp<volVectorField> dCdb(const label varID) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
