/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::CHFModels::Tatsumoto

Description
    A critical heat flux (CHF) sub-cooling correlation model.

Usage
    Example of the model specification:
    \verbatim
    CHFSubCoolModel
    {
        // Mandatory entries
        type            Tatsumoto;

        // Optional entries
        K               <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: Tatsumoto                | word | yes  | -
      K         | Model coefficient                   | scalar | no | 0.23
    \endtable

SourceFiles
    Tatsumoto.C

\*---------------------------------------------------------------------------*/

#ifndef Tatsumoto_H
#define Tatsumoto_H

#include "CHFSubCoolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class Tatsumoto Declaration
\*---------------------------------------------------------------------------*/

class Tatsumoto
:
    public CHFSubCoolModel
{
    // Private Data

        //- Model coefficient
        scalar K_;


    // Private Member Functions

        //- No copy construct
        Tatsumoto(const Tatsumoto&) = delete;

        //- No copy assignment
        void operator=(const Tatsumoto&) = delete;


public:

    //- Runtime type information
    TypeName("Tatsumoto");


    // Constructors

        //- Construct from a dictionary
        Tatsumoto(const dictionary& dict);


    //- Destructor
    virtual ~Tatsumoto() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> CHFSubCool
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

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
