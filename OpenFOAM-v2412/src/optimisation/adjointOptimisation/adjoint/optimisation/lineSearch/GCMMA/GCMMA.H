/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::GCMMA

Description
    Class implementing a workaround for using the Globally converging wrapper
    of the MMA update method. Implemented as a lineSearch, even though the
    search direction changes within each inner iteration of GCMMA, due the
    existing framework for re-solving the primal equations and checking the
    convergence of a merit function.

    Reference:
    \verbatim
        Svanberg, K. (2002)
        A class of globally convergent optimization methods based on
        conservative convex separable approximations.
        SIAM Journal of Optimization, 12), 555-573.
        https://doi.org/10.1137/S1052623499362822
    \endverbatim

SourceFiles
    GCMMA.C

\*---------------------------------------------------------------------------*/

#ifndef GCMMA_H
#define GCMMA_H

#include "lineSearch.H"
#include "MMA.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class GCMMA Declaration
\*---------------------------------------------------------------------------*/

class GCMMA
:
    public lineSearch
{
protected:

    // Protected data

        //- Cast of the update method to MMA
        MMA& mma_;

        //- File with rho values in each inner loop
        OFstream GCMMAFile_;

        //- File with CPU cost in each inner loop
        OFstream costFile_;

        //- Total cost, measured in Equivalent Flow Solutions, up to this cycle
        label cost_;

        //- Write the header of the output files
        bool writeHeader_;


    // Protected Member Functions

        //- Write rho and objective/constraint values & approx to files
        void writeToFiles(bool isConverged);


private:

    // Private Member Functions

        //- No copy construct
        GCMMA(const GCMMA&) = delete;

        //- No copy assignment
        void operator=(const GCMMA&) = delete;


public:

    //- Runtime type information
    TypeName("GCMMA");


    // Constructors

        //- Construct from components
        GCMMA
        (
            const dictionary& dict,
            const Time& time,
            updateMethod& UpdatheMethod
        );


    // Destructor
    virtual ~GCMMA() = default;


    // Member Functions

       //- Check whether linear search has converged
       virtual bool converged();

       //- Actually computes a new direction entirely,
       //- targeting the satisfaction of the GCMMA condition
       virtual void updateStep();

       //- Replace the correction with the one computed in updateStep
       virtual void updateCorrection(scalarField& correction);

       // Member operators

           //- Increment iteration number and store merit value corresponding to
           //- the previous optimisation cycle
           //virtual GCMMA& operator++();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
