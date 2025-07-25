/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::wallBoilingModels::departureDiameterModels::KocamustafaogullariIshii

Description
    A correlation for bubble departure diameter modelling
    based on Kocamustafaogullari-Ishii (1983) for boiling flows.

    Reference:
    \verbatim
        Kocamustafaogullari, G., & Ishii, M. (1983).
        Interfacial area and nucleation site density in boiling systems.
        International Journal of Heat and Mass Transfer, 26(9), 1377-1387.
        DOI:10.1016/S0017-9310(83)80069-6
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    departureDiameterModel
    {
        // Mandatory entries
        type            KocamustafaogullariIshii;
        phi             <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: KocamustafaogullariIshii | word | yes  | -
      phi       | Contact angle [deg]                 | scalar | yes  | -
    \endtable

SourceFiles
    KocamustafaogullariIshii.C

\*---------------------------------------------------------------------------*/

#ifndef KocamustafaogullariIshii_H
#define KocamustafaogullariIshii_H

#include "departureDiameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{

/*---------------------------------------------------------------------------*\
                  Class KocamustafaogullariIshii Declaration
\*---------------------------------------------------------------------------*/

class KocamustafaogullariIshii
:
    public departureDiameterModel
{
    // Private Data

        //- Contact angle
        scalar phi_;


    // Private Member Functions

        //- No copy construct
        KocamustafaogullariIshii(const KocamustafaogullariIshii&) = delete;

        //- No copy assignment
        void operator=(const KocamustafaogullariIshii&) = delete;


public:

    //- Runtime type information
    TypeName("KocamustafaogullariIshii");

    // Constructors

        //- Construct from a dictionary
        KocamustafaogullariIshii(const dictionary& dict);


    //- Destructor
    virtual ~KocamustafaogullariIshii() = default;


    // Member Functions

        //- Calculate and return the departure diameter field
        virtual tmp<scalarField> dDeparture
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

} // End namespace departureDiameterModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
