/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::radiation::greyMeanSolidAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    greyMeanSolidAbsorptionEmission radiation absorption and emission
    coefficients for solid mixture

SourceFiles
    greyMeanSolidAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef greyMeanSolidAbsorptionEmission_H
#define greyMeanSolidAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "solidThermo.H"
#include "basicSpecieMixture.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class greyMeanSolidAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class greyMeanSolidAbsorptionEmission
:
    public absorptionEmissionModel
{
private:

    // Private data

        //- Enumering of radiative properties
        enum radiativeProperties
        {
            absorptivity,
            emissivity
        };

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const solidThermo& thermo_;

        //- Hash table of species names
        HashTable<label> speciesNames_;

        //- Basic multicomponent mixture
        const basicSpecieMixture& mixture_;

        //- List of solid species data
        List<FixedList<scalar, 2>> solidData_;


    // Private member functions

       //- Calculate the volumetric fraction of Yj
       tmp<scalarField> X(const word specie) const;

       //- Calculate the property mixing
       tmp<volScalarField> calc(const label) const;


public:

    //- Runtime type information
    TypeName("greyMeanSolidAbsorptionEmission");


    // Constructors

        //- Construct from components
        greyMeanSolidAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~greyMeanSolidAbsorptionEmission() = default;


    // Member Functions

        //- Absorption coefficient for continuous phase
        tmp<volScalarField> aCont(const label bandI = 0) const;

        //- Emission coefficient for continuous phase
        tmp<volScalarField> eCont(const label bandI = 0) const;


    // Member Functions

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
