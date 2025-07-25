/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::radiation::localDensityAbsorptionEmission

Description
    Constant radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    localDensityAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_localDensityAbsorptionEmission_H
#define radiation_localDensityAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
               Class localDensityAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class localDensityAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Name of phase fields
        List<word> alphaNames_;

        //- Absorption coefficients
        List<scalar> aCoeff_;

        //- Emission coefficients
        List<scalar> eCoeff_;

        //- Emission contributions
        List<scalar> ECoeff_;


    // Private Member Functions

        //- Retrieve the phase field from the database
        const volScalarField& alpha(word) const;


public:

    //- Runtime type information
    TypeName("localDensityAbsorptionEmission");


    // Constructors

        //- Construct from components
        localDensityAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~localDensityAbsorptionEmission() = default;


    // Member Functions

        //- Absorption coefficient for continuous phase
        tmp<volScalarField> aCont(const label bandI = 0) const;

        //- Emission coefficient for continuous phase
        tmp<volScalarField> eCont(const label bandI = 0) const;

        //- Emission contribution for continuous phase
        tmp<volScalarField> ECont(const label bandI = 0) const;


        //- Is grey
        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
