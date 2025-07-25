/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::lineSearch

Description
    Abstract base class for line search methods

SourceFiles
    lineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef lineSearch_H
#define lineSearch_H

#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "stepUpdate.H"
#include "updateMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class lineSearch Declaration
\*---------------------------------------------------------------------------*/

class lineSearch
{
protected:

    // Protected data

        //- Subdict within updateMethod
        const dictionary dict_;

        //- IOdictionary under time/uniform for continuation
        IOdictionary lineSearchDict_;

        //- Directional derivative of the merit function
        scalar directionalDeriv_;

        //- Update direction
        scalarField direction_;

        //- Old merit value from this opt cycle
        scalar oldMeritValue_;

        //- New merit value from this opt cycle
        scalar newMeritValue_;

        //- Merit directional deriv from the previous opt cycle
        scalar prevMeritDeriv_;

        //- Correction multiplier at the first step of line search
        scalar initialStep_;

        //- Minimum allowed correction multiplier
        scalar minStep_;

        //- Correction multiplier
        scalar step_;

        //- Optimisation cycle
        label iter_;

        //- Inner line search iteration
        label innerIter_;

        //- Maximum line search iterations
        label maxIters_;

        //- Whether to extrapolate the correction multiplier for
        //- this optimisation cycle based on the previous ones.
        //  Useful for non-quasi Newton methods
        bool extrapolateInitialStep_;

        //- Mechanism to update method if line search conditions are not set
        autoPtr<stepUpdate> stepUpdate_;

        //- Reference to the update method related to the line search
        updateMethod& updateMethod_;


    // Protected Member Functions

        //- Optional coeffs dict
        const dictionary& coeffsDict();


private:

    // Private Member Functions

        //- No copy construct
        lineSearch(const lineSearch&) = delete;

        //- No copy assignment
        void operator=(const lineSearch&) = delete;


public:

    //- Runtime type information
    TypeName("lineSearch");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            lineSearch,
            dictionary,
            (
                const dictionary& dict,
                const Time& time,
                updateMethod& UpdateMethod
            ),
            (dict, time, UpdateMethod)
        );


    // Constructors

        //- Construct from components
        lineSearch
        (
            const dictionary& dict,
            const Time& time,
            updateMethod& UpdateMethod
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<lineSearch> New
        (
            const dictionary& dict,
            const Time& time,
            updateMethod& UpdateMethod
        );


    //- Destructor
    virtual ~lineSearch() = default;


    // Member Functions

        //- Set directional derivative
        virtual void setDeriv(const scalar deriv);

        //- Set new directional derivative.
        //  Does nothing in base. Only used by methods that require the
        //  gradient information to be computed at the new point
        virtual void setNewDeriv(const scalar deriv);

        //- Set direction
        inline void setDirection(const scalarField& direction)
        {
            direction_ = direction;
        }

        //- Set new objective value
        void setNewMeritValue(const scalar value);

        //- Set old objective value
        void setOldMeritValue(const scalar value);

        //- Reset step to initial value
        virtual void reset();

        //- Get inner line search iteration
        inline label innerIter() const
        {
            return innerIter_;
        }

        //- Get max number of iterations
        inline label maxIters() const
        {
            return maxIters_;
        }

        //- Get current step
        inline scalar step() const
        {
            return step_;
        }

        //- Return the correction of the design variables
        virtual bool converged() = 0;

        //- Update the line search step based on the specific line search
        //- strategy, e.g. bisection, quadratic fit, etc.
        virtual void updateStep() = 0;

        //- Update the step using the supplied value
        virtual void updateStep(const scalar newStep);

        //- Update the correction.
        //  Multiplies with step in base
        virtual void updateCorrection(scalarField& correction);

        //- Return true if lineSearch should continue and if so increment inner
        //  iteration
        virtual bool loop();

        //- Does line search need to update the gradient?
        virtual bool computeGradient() const;

        //- Execute steps at the end of the line search iterations
        virtual void postUpdate();


    // Member operators

        //- Increment iteration number and store merit value corresponding to
        //- the previous optimisation cycle
        virtual lineSearch& operator++();

        //- Postfix increment. Necessary for compilation
        virtual lineSearch& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
