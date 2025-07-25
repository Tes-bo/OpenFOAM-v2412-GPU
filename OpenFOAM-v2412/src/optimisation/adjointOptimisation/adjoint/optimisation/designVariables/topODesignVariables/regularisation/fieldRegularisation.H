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
    fieldRegularisation

Description
    Regularises the given field using a PDE smoother and then sharpens
    (projects) the result using a smoothed Heaviside function. Used in topology
    optimisation to regularise the designed geometries and obtain 'grid
    independent' solutions.

    Reference:
    \verbatim
        Largely follows the PDE-based smoothing of

            Lazarov, B. S., & Sigmund O. (2010).
            Filters in topology optimization based on Helmholtz-type
            differential equations.
            International Journal for Numerical Methods in Engineering,
            86(6), 765-781.
            https://doi.org/10.1002/nme.3072
    \endverbatim

SourceFiles
    fieldRegularisation.C

\*---------------------------------------------------------------------------*/

#ifndef fieldRegularisation_H
#define fieldRegularisation_H

#include "fvMesh.H"
#include "volFields.H"
#include "topOZones.H"
#include "regularisationPDE.H"
#include "regularisationRadius.H"
#include "topOInterpolationFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fieldRegularisation Declaration
\*---------------------------------------------------------------------------*/

class fieldRegularisation
{

protected:

    // Protected data

        fvMesh& mesh_;
        dictionary dict_;

        //- Cell zones related to topology optimisation
        const topOZones& zones_;

        //- Perform regulaisation on alpha before inputing it on beta?
        bool regularise_;

        //- Perform the projection (sharpening) step?
        bool project_;

        //- Smoothing radius
        //  May be isotropic or differ per spatial direction
        autoPtr<regularisationRadius> radius_;

        //- Alpha field (design variables of topology optimisation)
        const scalarField& alpha_;

        //- The regularised alpha field, if regulatisation is performed
        autoPtr<volScalarField> alphaTilda_;

        //- Function used to sharpen the field after regularisation
        autoPtr<topOInterpolationFunction> sharpenFunction_;

        //- PDE used for the regularisation
        autoPtr<regularisationPDE> regularisationPDE_;

        //- Argument of the beta field.
        //  Is either alpha_ (no regularisation)
        //  or alphaTilda_ (with regularisation)
        const scalarField& betaArg_;

        //- Whether to apply a fixedValue BC or zeroGradient one to alphaTilda,
        //- when regularisation is performed
        bool growFromWalls_;

        //- Beta is the field used for all interpolations between fluid and
        //- solid in topology optimisation
        volScalarField beta_;


private:

    // Private Member Functions

        //- No copy construct
        fieldRegularisation(const fieldRegularisation&);

        //- No copy assignment
        void operator=(const fieldRegularisation&);


public:

    //- Runtime type information
    TypeName("fieldRegularisation");


    // Constructors

        //- Construct from components
        fieldRegularisation
        (
            fvMesh& mesh,
            const scalarField& alpha,
            const topOZones& zones,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldRegularisation() = default;


    // Member Functions

        //- Return beta field
        inline const volScalarField& beta() const
        {
            return beta_;
        }

        //- Grow beta from walls
        inline bool growFromWalls() const
        {
            return growFromWalls_;
        }

        //- Should regularisation be executed
        inline bool shouldRegularise() const
        {
            return regularise_;
        }

        //- Update the beta field.
        //  Performs regulaisation of alpha_, if necessary
        virtual void updateBeta();

        //- Regularise an externally provided radius
        void regularise
        (
            const scalarField& source,
            scalarField& result,
            const bool isTopoField,
            const regularisationRadius& radius
        );

        //- Regularise field (or sensitivities) using a Laplace PDE
        void regularise
        (
            const scalarField& source,
            scalarField& result,
            const bool isTopoField
        );

        //- Update part of fieldRegularisation to the sensitivitiy derivatives
        void postProcessSens(scalarField& sens);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
