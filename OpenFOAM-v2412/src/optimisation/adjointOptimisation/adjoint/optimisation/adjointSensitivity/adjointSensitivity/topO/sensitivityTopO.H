/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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


Description
    Calculation of adjoint based sensitivities for topology optimisation.
    This returns just the field part, with  contributions from regularisation
    and projection added by topODesignVariables.

Class
    Foam::sensitivityTopO

SourceFiles
    sensitivityTopO.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityTopO_H
#define sensitivityTopO_H

#include "adjointSensitivity.H"
#include "fvOptions.H"
#include "topOZones.H"
#include "adjointEikonalSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensitivityTopO Declaration
\*---------------------------------------------------------------------------*/

class sensitivityTopO
:
    public adjointSensitivity
{

protected:

    // Protected data

        //- Zones related to topology optimisation
        topOZones zones_;

        //- Name used as the argument for the post-processing of the
        //- sensitivities through fvOptions
        word designVariablesName_;


    // Protected Member Functions

        //- Zero sensitivities in fixed porous zones
        void zeroSensInFixedPorousZones(scalarField& sens);


private:

    // Private Member Functions

        //- No copy construct
        sensitivityTopO(const sensitivityTopO&);

        //- No copy assignment
        void operator=(const sensitivityTopO&);


public:

    //- Runtime type information
    TypeName("topO");


    // Constructors

        //- Construct from components
        sensitivityTopO
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivityTopO() = default;


    // Member Functions

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Accumulate sensitivity integrands
        virtual void accumulateIntegrand(const scalar dt);

        //- Assemble sensitivities
        virtual void assembleSensitivities
        (
            autoPtr<designVariables>& designVars
        );

        //- Add part of the sensitivities coming from fvOptions
        static void postProcessSens
        (
            scalarField& sens,
            scalarField& auxSens,
            fv::options& fvOptions,
            const word& fieldName,
            const word& designVariablesName
        );

        //- Add part of the sensitivities coming from fvOptions
        void postProcessSens
        (
            scalarField& sens,
            scalarField& auxSens,
            const word& fieldName
        );

        //- Potentially manipulate auxSens within the fvOption before adding
        //- the part related to the design variables
        template<class Type, template<class> class PatchField, class GeoMesh>
        static void postProcessAuxSens
        (
             const GeometricField<Type, PatchField, GeoMesh>& primal,
             const GeometricField<Type, PatchField, GeoMesh>& adjoint,
             fv::options& fvOptions,
             const word& fieldName,
             scalarField& auxSens
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sensitivityTopOTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
