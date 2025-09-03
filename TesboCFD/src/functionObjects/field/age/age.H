/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::functionObjects::age

Description
    Calculates and writes out the time taken for a particle to travel from an
    inlet to the location. Solves the following equation when incompressible:

    \f[
        \div (\phi t) = 1
    \f]

    where
    \vartable
        t    | Age [s]
        \phi | Volumetric flux [m^3/s]
    \endvartable

    Boundary conditions are generated automatically as \c zeroGradient
    on all walls and \c inletOutlet everywhere else.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    age1
    {
        // Mandatory entries
        type            age;
        libs            (fieldFunctionObjects);

        // Optional entries
        phi             <word>;
        rho             <word>;
        schemesField    <word>;
        tolerance       <scalar>;
        nCorr           <int>;
        diffusion       <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                        | Type | Reqd | Deflt
      type      | Type name: age                     | word | yes  | -
      libs      | Library name: fieldFunctionObjects | word | yes  | -
      phi       | Name of flux field                 | word | no   | phi
      rho       | Name of density field              | word | no   | rho
      schemesField | Name of the field from which schemes are taken <!--
                -->                                  | word | no   | age
      tolerance | Solver residual control            | scalar | no | 1e-5
      nCorr     | Maximum number of correctors       | int  | no   | 5
      diffusion | Flag to turn on/off the diffusion term | bool  | no | false
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

See also
  - Foam::functionObjects::comfort

SourceFiles
    age.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_age_H
#define functionObjects_age_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class age Declaration
\*---------------------------------------------------------------------------*/

class age
:
    public fvMeshFunctionObject
{
    // Private Data

        //- The name of the flux field
        word phiName_;

        //- The name of the density field
        word rhoName_;

        //- Name of field from which schemes are taken
        word schemesField_;

        //- Convergence tolerance
        scalar tolerance_;

        //- Number of corrections
        int nCorr_;

        //- Flag to turn on/off the diffusion term
        bool diffusion_;


    // Private Member Functions

        //- The list of patch types for the age field
        wordList patchTypes() const;

        //- Return true if convergence is reached
        bool converged(const int nCorr, const scalar initialResidual) const;


public:

    //- Runtime type information
    TypeName("age");


    // Constructors

        //- Construct from Time and dictionary
        age
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~age() = default;


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
