/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fv::radialActuationDiskSource

Group
    grpFvOptionsSources

Description
    Applies sources on velocity field (i.e. \c U) within a specified region
    to enable actuator disk models with radial thrust distribution for
    aero/hydro thrust loading of horizontal axis turbines on surrounding
    flow field in terms of energy conversion processes.

    Sources applied to:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Thrust, which is computed by \link actuationDiskSource.H \endlink is
    distributed as a function of actuator disk radius through a given forth
    order polynomial function:

    \f[
        T(r) = T (C_0 + C_1 r^2 + C_2 r^4)
    \f]
    where
    \vartable
      T         | Thrust magnitude computed by actuationDiskSource setup
      T(r)      | Thrust magnitude as a function of "r"
      r         | Local actuator disk radius
      C_*       | Polynomial coefficients
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    radialActuationDiskSource1
    {
        // Mandatory entries (unmodifiable)
        type            radialActuationDiskSource;

        // Mandatory entries (runtime modifiable)
        coeffs          (0.1 0.5 0.01);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type | Reqd | Dflt
      type     | Type name: radialActuationDiskSource   | word | yes  | -
      coeffs   | Radial distribution function coefficients | vector | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink
     - \link writeFile.H \endlink
     - \link Function1.H \endlink
     - \link actuationDiskSource.H \endlink

Note
  - \c radialActuationDiskSource computes only
    \c Froude force computation method of \c actuationDiskSource.

SourceFiles
    radialActuationDiskSource.C
    radialActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef radialActuationDiskSource_H
#define radialActuationDiskSource_H

#include "actuationDiskSource.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class radialActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class radialActuationDiskSource
:
    public actuationDiskSource
{
    // Private Data

        //- Coefficients for the radial distribution
        FixedList<scalar, 3> radialCoeffs_;


    // Private Member Functions

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addRadialActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        );


public:

    //- Runtime type information
    TypeName("radialActuationDiskSource");


    // Constructors

        //- Construct from components
        radialActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        radialActuationDiskSource(const radialActuationDiskSource&) = delete;

        //- No copy assignment
        void operator=(const radialActuationDiskSource&) = delete;


    //- Destructor
    virtual ~radialActuationDiskSource() = default;


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
    #include "radialActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
