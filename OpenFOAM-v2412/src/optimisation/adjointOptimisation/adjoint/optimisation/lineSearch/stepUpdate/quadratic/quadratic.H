/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::quadratic

Description
    Fits a quadratic polynomial of the merit function as a function of step
    and finds the "optimal" value

SourceFiles
    quadratic.C

\*---------------------------------------------------------------------------*/

#ifndef quadratic_H
#define quadratic_H

#include "stepUpdate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class quadratic Declaration
\*---------------------------------------------------------------------------*/

class quadratic
:
    public stepUpdate
{
protected:

    // Protected data

        scalar minRatio_;
        scalar firstMeritValue_;
        scalar secondMeritValue_;
        scalar meritDerivative_;


private:

    // Private Member Functions

        //- No copy construct
        quadratic(const quadratic&) = delete;

        //- No copy assignment
        void operator=(const quadratic&) = delete;


public:

    //- Runtime type information
    TypeName("quadratic");


    // Constructors

        //- Construct from components
        quadratic(const dictionary& dict);


    //- Destructor
    virtual ~quadratic() = default;


    // Member Functions

       //- Update step in given direction
       virtual void updateStep(scalar& step);

       //- Set objective derivative
       virtual void setDeriv(const scalar deriv);

       //- Set new merit value
       virtual void setNewMeritValue(const scalar value);

       //- Set old merit value
       virtual void setOldMeritValue(const scalar value);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
