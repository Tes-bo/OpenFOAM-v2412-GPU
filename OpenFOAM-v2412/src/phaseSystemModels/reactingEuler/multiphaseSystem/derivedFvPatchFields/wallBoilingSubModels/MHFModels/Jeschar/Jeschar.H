/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::MHFModels::Jeschar

Description
    A model for minimum heat flux based on
    Jeschar et al. (1992) for boiling flows.

    Reference:
    \verbatim
        Jeschar, R., Specht, E., & Köhler, C. (1992).
        Heat transfer during cooling of heated
        metallic objects with evaporating liquids.
        In Theory and Technology of Quenching (pp. 73-92).
        Springer, Berlin, Heidelberg.
        DOI:10.1007/978-3-662-01596-4_4
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    MHFModel
    {
        // Mandatory entries
        type            Jeschar;

        // Optional entries
        Kmhf            <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: Jeschar                  | word | yes  | -
      Kmhf      | Burn out factor                   | scalar | no   | 1.0
    \endtable

SourceFiles
    Jeschar.C

\*---------------------------------------------------------------------------*/

#ifndef Jeschar_H
#define Jeschar_H

#include "MHFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class Jeschar Declaration
\*---------------------------------------------------------------------------*/

class Jeschar
:
    public MHFModel
{
    // Private Data

        //- Burn out factor
        scalar Kmhf_;


    // Private Member Functions

        //- No copy construct
        Jeschar(const Jeschar&) = delete;

        //- No copy assignment
        void operator=(const Jeschar&) = delete;


public:

    //- Runtime type information
    TypeName("Jeschar");


    // Constructors

        //- Construct from a dictionary
        Jeschar(const dictionary& dict);


    //- Destructor
    virtual ~Jeschar() = default;


    // Member Functions

        //- Calculate and return the minimum heat flux
        virtual tmp<scalarField> MHF
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CHFModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
