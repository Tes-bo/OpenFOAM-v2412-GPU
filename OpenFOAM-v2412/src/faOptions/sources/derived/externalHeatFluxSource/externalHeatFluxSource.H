/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::fa::externalHeatFluxSource

Group
    grpFaOptionsSources

Description
    Applies a heat flux condition for a specified \c faMesh region
    to temperature on an external wall for compressible flows
    in one of three modes:

      - fixed power: supply \c Q
      - fixed heat flux: supply \c q
      - fixed heat transfer coefficient: supply \c h and \c Ta

    where
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    The ambient temperature \c Ta is specified
    as a \c Foam::Function1 of time but uniform in space.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    externalHeatFluxSource1
    {
        // Mandatory entries (unmodifiable)
        type            externalHeatFluxSource;

        // Mandatory entries (runtime modifiable)
        mode            <mode>;

        // Optional entries (runtime modifiable)
        T               <Tname>;
        emissivity      0;

        // Conditional mandatory entries (runtime modifiable)

            // when mode=power
            Q           1.0;

            // when mode=flux
            q           1.0;

            // when mode=coefficient
            h           1.0;
            Ta          <Function1>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                          | Type  | Reqd | Dflt
      type      | Type name: externalHeatFluxSource    | word  | yes  | -
      mode      | Mode of heat flux condition          | word  | yes  | -
      T         | Name of operand temperature field    | word  | no   | T
      emissivity | Surface emissivity for radiative flux to ambient <!--
                 -->                                   | scalar | no  | 0
      Q         | Fixed heat power [W]                 | Function1 | cndtnl | -
      q         | Fixed heat flux [W/m2]               | Function1 | cndtnl | -
      h         | Heat transfer coefficient [W/m^2/K]  | Function1 | cndtnl | -
      Ta        | Ambient temperature [K]           | Function1 | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
     - \link faOption.H \endlink
     - \link faceSetOption.H \endlink

    Options for the \c mode entry:
    \verbatim
      power       | Use fixed power (supply Q)
      flux        | Use fixed heat flux (supply q)
      coefficient | Use fixes heat transfer coefficient (supply h and T)
    \endverbatim

See also
  - Foam::Function1

SourceFiles
    externalHeatFluxSource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_externalHeatFluxSource_H
#define Foam_fa_externalHeatFluxSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                    Class externalHeatFluxSource Declaration
\*---------------------------------------------------------------------------*/

class externalHeatFluxSource
:
    public fa::faceSetOption
{
public:

    // Public Enumeration

        //- Options for the heat transfer condition mode
        enum operationMode
        {
            fixedPower,                 //!< Fixed heat power [W]
            fixedHeatFlux,              //!< Fixed heat flux [W/m2]
            fixedHeatTransferCoeff,     //!< Fixed heat transfer coefficient
        };

        //- Names for operationMode
        static const Enum<operationMode> operationModeNames;


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Name of temperature field
        word TName_;

        //- Heat power [W]
        autoPtr<Function1<scalar>> Q_;

        //- Heat flux [W/m2]
        autoPtr<Function1<scalar>> q_;

        //- Heat transfer coefficient [W/m2K]
        autoPtr<Function1<scalar>> h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("externalHeatFluxSource");


    // Constructors

        //- Construct from explicit source name and mesh
        externalHeatFluxSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        externalHeatFluxSource(const externalHeatFluxSource&) = delete;

        //- No copy assignment
        void operator=(const externalHeatFluxSource&) = delete;


    //- Destructor
    virtual ~externalHeatFluxSource() = default;


    // Member Functions

        // Evaluation

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
