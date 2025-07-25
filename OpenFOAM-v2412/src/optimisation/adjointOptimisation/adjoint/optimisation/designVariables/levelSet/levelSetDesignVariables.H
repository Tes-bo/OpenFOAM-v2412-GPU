/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 PCOpt/NTUA
    Copyright (C) 2022-2023 FOSS GP
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
    Foam::levelSetDesignVariables

Description
    Signed distance field design variables for level-set based topology
    optimization (topO).

    A great reference for various aspects of level-set-based topO can be found
    in

    Reference:
    \verbatim
        van Dijk, N.P., Maute, K., Langelaar, M., & & van Keulen, F. (2013).
        Level-set methods for structural topology optimization: a review.
        Structural and Multidisciplinary Optimization, 48, 437–472.
        https://doi.org/10.1007/s00158-013-0912-y
    \endverbatim

    The implemented approach borrows from the references therein, but does not
    follow any of them exactly. The underlying field of the design variables
    is not the level-set field but affects the latter through:
    a) a regularisation approach similar to that used in porosity-based topO
    (see the regularisation class and references therein) and
    b) a re-initialisation strategy. Since solving a modified Hamilton-Jacobi
    equation poses a number of numerical challenges, the signed distance field
    is obtained by identifying the zero level-set contour using the
    cutFaceIso, cutCellIso infrastructure (see references therein for the
    methodology) and using these "cut-faces" as seeds for meshWave.

SourceFiles
    levelSetDesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetDesignVariables_H
#define levelSetDesignVariables_H

#include "topOVariablesBase.H"
#include "designVariables.H"
#include "regularisationRadius.H"
#include "regularisationPDE.H"
#include "topOInterpolationFunction.H"
#include "adjointSensitivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class levelSetDesignVariables Declaration
\*---------------------------------------------------------------------------*/

class levelSetDesignVariables
:
    public topOVariablesBase,
    public designVariables
{
protected:

    // Protected Data Members

        //- The regularisation radius
        autoPtr<regularisationRadius> radius_;

        //- Regularisation mechanism
        autoPtr<regularisationPDE> regularisation_;

        //- The regularised field
        volScalarField aTilda_;

        //- The signed distances field
        volScalarField signedDistances_;

        //- Function to transorm signed distances to the indicator field beta_
        autoPtr<topOInterpolationFunction> interpolation_;

        //- The indicator field
        volScalarField beta_;

        //- Fix aTilda values in fixed{Zero}PorousZones and IOcells
        bool fixATildaValues_;

        //- Write all fields related to the distance calculation (debugging)
        bool writeAllDistanceFields_;

        //- Mesh faces acting as the source of MeshWave
        labelList changedFaces_;

        //- Seed distances to MeshWave and cell distances
        //  The data carried by each wallPoints corresponds to the origin
        //  mesh face ID
        List<wallPointData<label>> changedFacesInfo_;
        List<wallPointData<label>> allCellInfo_;


    // Protected Member Functions

        //- Read the design variables field
        void readField();

        //- Apply fixed values in certain zones
        void applyFixedPorosityValues();

        //- Determine which design variables are active
        void setActiveDesignVariables(bool activeIO = false);

        //- Make aTilda a signed distance field based on the zero iso-surface
        //- of the current aTilda field
        void updateSignedDistances();

        //- Update the beta field
        void updateBeta();


private:

    // Private Member Functions

        //- No copy construct
        levelSetDesignVariables(const levelSetDesignVariables&) = delete;

        //- No copy assignment
        void operator=(const levelSetDesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("levelSet");


    // Constructors

        //- Construct from dictionary
        levelSetDesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<levelSetDesignVariables> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~levelSetDesignVariables() = default;


    // Member Functions

        //- Const reference to the beta field
        virtual const volScalarField& beta() const;

        //- Update design variables based on a given correction
        virtual void update(scalarField& correction);

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction);

        //- Whether to use global sum when computing matrix-vector products
        //- in update methods
        virtual bool globalSum() const;

        //- Add contributions from the adjoint to the regularization PDE,
        //- the derivative of the interpolation function and the mesh volume
        //- to the final sensitivities
        virtual tmp<scalarField> assembleSensitivities
        (
            adjointSensitivity& adjointSens
        );

        //- Write useful quantities to files
        virtual void writeDesignVars();

        //- The writeData function required by the regIOobject write operation
        virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
