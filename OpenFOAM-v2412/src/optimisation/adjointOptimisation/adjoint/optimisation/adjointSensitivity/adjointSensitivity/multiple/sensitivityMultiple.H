/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::sensitivityMultiple

Description
    Calculation of adjoint based sensitivities of multiple types

SourceFiles
    sensitivityMultiple.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityMultiple_H
#define sensitivityMultiple_H

#include "adjointSensitivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sensitivityMultiple Declaration
\*---------------------------------------------------------------------------*/

class sensitivityMultiple
:
    public adjointSensitivity
{
protected:

    // Protected data

        wordList sensTypes_;

        PtrList<adjointSensitivity> sens_;


private:

    // Private Member Functions

        //- No copy construct
        sensitivityMultiple(const sensitivityMultiple&) = delete;

        //- No copy assignment
        void operator=(const sensitivityMultiple&) = delete;


public:

    //- Runtime type information
    TypeName("multiple");


    // Constructors

        //- Construct from components
        sensitivityMultiple
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivityMultiple() = default;


    // Member Functions

       //- Read dict if changed
       virtual bool readDict(const dictionary& dict);

       //- Accumulate sensitivity integrands
       virtual void accumulateIntegrand(const scalar dt);

       //- Assemble sensitivities
       virtual void assembleSensitivities(autoPtr<designVariables>& designVars);

       //- Calculates sensitivities at wall surface points
       const scalarField& calculateSensitivities
       (
           autoPtr<designVariables>& designVars
       );

       //- Zero sensitivity fields and their constituents
       virtual void clearSensitivities();

       //- Write sensitivities to file
       virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
