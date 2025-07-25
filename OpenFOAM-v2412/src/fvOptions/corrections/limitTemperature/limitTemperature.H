/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::fv::limitTemperature

Group
    grpFvOptionsCorrections

Description
    Corrects temperature field (i.e. \c T) within a specified
    region by applying limits between a given minimum and maximum.

    Corrections applied to:
    \verbatim
      T    | Temperature                   [K]
    \endverbatim

    Required fields:
    \verbatim
      T    | Temperature                   [K]
      e/h  | Internal energy/Enthalphy     [m2/s2]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    limitTemperature1
    {
        // Mandatory entries (unmodifiable)
        type            limitTemperature;

        // Mandatory entries (runtime modifiable)
        min             <minValue>;
        max             <maxValue>;

        // Optional entries (runtime modifiable)
        phase           <phaseName>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Dflt
      type      | Type name: limitTemperature   | word   | yes   | -
      min       | Minimum temperature limit [K] | scalar | yes   | -
      max       | Maximum temperature limit [K] | scalar | yes   | -
      phase     | Name of operand phase field   | word   | no    | ""
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

See also
  - Foam::fv::fixedTemperatureConstraint

SourceFiles
    limitTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef fv_limitTemperature_H
#define fv_limitTemperature_H

#include "cellSetOption.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class limitTemperature Declaration
\*---------------------------------------------------------------------------*/

class limitTemperature
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;

        //- Optional phase name
        word phase_;


    // Protected Member Functions

        //- Write file header information
        void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("limitTemperature");


    // Constructors

        //- Construct from components
        limitTemperature
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        limitTemperature(const limitTemperature&) = delete;

        //- No copy assignment
        void operator=(const limitTemperature&) = delete;


    //- Destructor
    virtual ~limitTemperature() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the energy field
        virtual void correct(volScalarField& he);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
