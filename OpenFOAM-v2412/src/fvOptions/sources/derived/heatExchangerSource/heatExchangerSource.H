/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::fv::heatExchangerSource

Group
    grpFvOptionsSources

Description
    Heat exchanger source model for compressible flows, where the heat
    exchanger is modelled as an energy source using a selection of cells.

    Sources applied to either of the below, if exist:
    \verbatim
      e         | Internal energy                            [m2/s2]
      h         | Enthalphy                                  [m2/s2]
    \endverbatim

    Required fields:
    \verbatim
      T         | Temperature                                [K]
      U         | Velocity                                   [m/s]
      phi       | Mass flux                                  [kg/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <name>
    {
        // Mandatory entries
        type                     heatExchangerSource;
        model                    <word>;
        faceZone                 <word>;

        // Conditional entries

            // Option-1: when model == effectivenessTable

            // Option-2: when model == referenceTemperature

        // Optional entries
        U                        <word>;
        T                        <word>;
        phi                      <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Deflt
      type      | Type name: heatExchangerSource    | word  | yes  | -
      model     | Name of the heat exchanger model  | word  | yes  | -
      faceZone  | Name of the faceZone at the heat exchanger inlet <!--
                -->                                 | word   | yes | -
      U         | Name of operand velocity field    | word   | no  | U
      T         | Name of operand temperature field | word   | no  | T
      phi       | Name of operand flux field        | word   | no  | phi
    \endtable

    Options for the \c model entry:
    \verbatim
      effectivenessTable    | Calculate heat exchange via an effectiveness table
      referenceTemperature  | Calculate heat exchange via a temperature table
    \endverbatim

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink
      - \link Function1.H \endlink

SourceFiles
    heatExchangerSource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_heatExchangerSource_H
#define Foam_fv_heatExchangerSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class heatExchangerModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class heatExchangerSource Declaration
\*---------------------------------------------------------------------------*/

class heatExchangerSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Heat exchanger model
        autoPtr<heatExchangerModel> heatExchangerModelPtr_;


public:

    //- Runtime type information
    TypeName("heatExchangerSource");


    // Constructors

        //- Construct from components
        heatExchangerSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        heatExchangerSource(const heatExchangerSource&) = delete;

        //- No copy assignment
        void operator=(const heatExchangerSource&) = delete;


    //- Destructor
    virtual ~heatExchangerSource() = default;


    // Member Functions

        //- Add explicit/implicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        )
        {
            NotImplemented;
        }

        //- Add explicit/implicit contribution
        //- to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
