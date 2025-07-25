/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::trimModel

Description
    Base class for trim models for handling
    blade characteristics and thrust-torque relations.

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        trimModel    <trimModel>;

        <trimModel>Coeffs
        {
            // Conditional mandatory entries (runtime modifiable)

                // when trimModel=fixedTrim
                theta0   5.0;
                theta1c  1.0;
                theta1s  4.0;
        }
    }

    where the entries mean:
    \table
      Property   | Description                          | Type  | Reqd | Dflt
      trimModel  | Trim model specification             | word  | yes  | -
      file  | Name of file containing profile characteristics | word   | no  | -
    \endtable

    Options for the \c trimModel entry:
    \verbatim
      fixed       | Fixed blade/angle characteristics
      targetForce | Target thrust/torque by using time-variant blade characs.
    \endverbatim

See also
  - Foam::fv::rotorDiskSource
  - Foam::fixedTrim
  - Foam::targetCoeffTrim

SourceFiles
    trimModel.C

\*---------------------------------------------------------------------------*/

#ifndef trimModel_H
#define trimModel_H

#include "rotorDiskSource.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class trimModel Declaration
\*---------------------------------------------------------------------------*/

class trimModel
{
protected:

    // Protected Data

        //- Reference to the rotor source model
        const fv::rotorDiskSource& rotor_;

        //- Name of model
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;


public:

    //- Run-time type information
    TypeName("trimModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            trimModel,
            dictionary,
            (
                const fv::rotorDiskSource& rotor,
                const dictionary& dict
            ),
            (rotor, dict)
        );


    // Constructors

        //- Construct from components
        trimModel
        (
            const fv::rotorDiskSource& rotor,
            const dictionary& dict,
            const word& name
        );


    // Selectors

        //- Return a reference to the selected trim model
        static autoPtr<trimModel> New
        (
            const fv::rotorDiskSource& rotor,
            const dictionary& dict
        );


    //- Destructor
    virtual ~trimModel() = default;


    // Member Functions

        //- Read
        virtual void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const = 0;

        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        ) = 0;

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
