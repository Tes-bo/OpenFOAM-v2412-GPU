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
    Foam::wallBoilingModels::CHFModels::Zuber

Description
    A critical heat flux (CHF) correlation model
    based on Zuber (1958) for boiling flows.

    Reference:
    \verbatim
        Zuber, N. (1958).
        On the stability of boiling heat transfer.
        Trans. Am. Soc. Mech. Engrs., 80.
        URL:https://www.osti.gov/biblio/4326542
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    CHFModel
    {
        // Mandatory entries
        type            Zuber;

        // Optional entries
        Cn              <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: Zuber                    | word | yes  | -
      Cn        | Model coefficient                 | scalar | no   | 0.131
    \endtable

SourceFiles
    Zuber.C

\*---------------------------------------------------------------------------*/

#ifndef Zuber_H
#define Zuber_H

#include "CHFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class Zuber Declaration
\*---------------------------------------------------------------------------*/

class Zuber
:
    public CHFModel
{
    // Private Data

        //- Model coefficient
        scalar Cn_;


    // Private Member Functions

        //- No copy construct
        Zuber(const Zuber&) = delete;

        //- No copy assignment
        void operator=(const Zuber&) = delete;


public:

    //- Runtime type information
    TypeName("Zuber");


    // Constructors

        //- Construct from a dictionary
        Zuber(const dictionary& dict);


    //- Destructor
    virtual ~Zuber() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> CHF
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        // I-O

            // Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CHFModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
