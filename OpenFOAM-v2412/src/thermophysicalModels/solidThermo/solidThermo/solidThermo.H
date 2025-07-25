/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::solidThermo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermo_H
#define solidThermo_H

#include "rhoThermo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermo
:
    public rhoThermo
{


public:

    //- Runtime type information
    TypeName("solidThermo");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict, const word& phaseName),
        (mesh, dict, phaseName)
    );

   //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        solidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        solidThermo
        (
            const fvMesh&,
            const dictionary& dict,
            const word& phaseName
        );

        //- Construct from mesh with dictionary name and without
        // constructing T phaseName
        solidThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );

        //- Return a pointer to a new solidThermo created from
        //  the solidThermophysicalProperties dictionary
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        //- Return a pointer to a new solidThermo created from
        //  local dictionary
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName=word::null
        );


        //- Return a pointer to a new solidThermo created from
        //  local dictionary and phase name
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~solidThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const = 0;

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const = 0;


        // Per patch calculation

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa
            (
                const label patchi
            ) const = 0;


    // I-O

        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
