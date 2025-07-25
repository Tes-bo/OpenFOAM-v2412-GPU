/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::fv::tabulatedAccelerationSource

Group
    grpFvOptionsSources

Description
    Solid-body 6-DoF acceleration source.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    tabulatedAccelerationSource1
    {
        // Mandatory entries (unmodifiable)
        type             tabulatedAccelerationSource;
        timeDataFileName "constant/acceleration-terms.dat";

        // Optional entries (unmodifiable)
        U                U;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type  | Reqd | Dflt
      type      | Type name: tabulatedAccelerationSource | word  | yes  | -
      timeDataFileName | Name of tabulated data file     | word  | yes  | -
      U         | Name of operand velocity field         | word  | no   | U
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

SourceFiles
    tabulatedAccelerationSource.C
    tabulatedAccelerationSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedAccelerationSource_H
#define tabulatedAccelerationSource_H

#include "fvOption.H"
#include "tabulated6DoFAcceleration.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class tabulatedAccelerationSource Declaration
\*---------------------------------------------------------------------------*/

class tabulatedAccelerationSource
:
    public fv::option
{
protected:

    // Protected Data

        //- Run-time selectable acceleration model
        tabulated6DoFAcceleration motion_;

        //- Name of operand velocity field
        word UName_;

        //- Gravitational accelaration
        dimensionedVector g0_;


private:

    // Private Member Functions

        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("tabulatedAccelerationSource");


    // Constructors

        //- Construct from components
        tabulatedAccelerationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        tabulatedAccelerationSource
        (
            const tabulatedAccelerationSource&
        ) = delete;

        //- No copy assignment
        void operator=(const tabulatedAccelerationSource&) = delete;


    //- Destructor
    virtual ~tabulatedAccelerationSource() = default;


    // Member Functions

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tabulatedAccelerationSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
