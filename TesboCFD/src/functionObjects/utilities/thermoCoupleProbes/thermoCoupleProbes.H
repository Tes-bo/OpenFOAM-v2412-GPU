/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::thermoCoupleProbes

Group
    grpUtilitiesFunctionObjects

Description
    Sample probe for temperature using a thermocouple.

    Uses the correlation:

    \f[
        Nu = 2.0 + \left(0.4 Re^{0.5} + 0.06 Re^{2/3}\right)*Pr^{0.4}
    \f]

Usage
    Example of function object specification:
    \verbatim
    probes
    {
        type            thermoCoupleProbes;
        libs            (utilityFunctionObjects);
        writeControl    timeStep;
        writeInterval   1;

        solver          rodas23;
        absTol          1e-12;
        relTol          1e-8;

        interpolationScheme cellPoint;

        // thermocouple properties
        rho             8908;
        Cp              440;
        d               1e-3;
        epsilon         0.85;

        radiationField  G;

        probeLocations
        (
            (0.5 0.5 0.5)
        );
        fields
        (
            T
        );
    }
    \endverbatim


SourceFiles
    thermoCoupleProbes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_thermoCoupleProbes_H
#define Foam_functionObjects_thermoCoupleProbes_H

#include "probes.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "basicThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class thermoCoupleProbes Declaration
\*---------------------------------------------------------------------------*/

class thermoCoupleProbes
:
    public probes,
    public ODESystem
{
protected:

    // Protected Data

        //- Thermocouple density
        scalar rho_;

        //- Thermocouple heat capacity
        scalar Cp_;

        //- Thermocouple diameter
        scalar d_;

        //- Thermocouple emissivity
        scalar epsilon_;

        //- Name of the velocity field
        word UName_;

        //- Name of the incident radiation field
        word radiationFieldName_;

        //- Fluid thermo reference
        const fluidThermo& thermo_;

        //- ODESolver
        autoPtr<ODESolver> odeSolver_;

        //- Cached thermocouple temperature
        scalarField Ttc_;


private:

    // Private Member Functions

        //- Sample/write
        template<class Type>
        void writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalar timeValue
        );

        //- No copy construct
        thermoCoupleProbes(const thermoCoupleProbes&) = delete;

        //- No copy assignment
        void operator=(const thermoCoupleProbes&) = delete;


public:

    //- Runtime type information
    TypeName("thermoCoupleProbes");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        thermoCoupleProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~thermoCoupleProbes() = default;


    // ODE functions (overriding abstract functions in ODE.H)

        //- Number of ODE's to solve
        virtual label nEqns() const;

        virtual void derivatives
        (
            const scalar x,
            const scalarField& y,
            scalarField& dydx
        ) const;

        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const;


    // Public Member Functions

        //- Read
        virtual bool read(const dictionary&);

        //- Execute. Evaluates the ODESolver
        virtual bool execute();

        //- Sample and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thermoCoupleProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
