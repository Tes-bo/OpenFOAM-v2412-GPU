/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::functionObjects::resolutionIndex

Group
    grpFieldFunctionObjects

Description
    Computes a single-mesh resolution index according to the specified model,
    which is used as a LES/DES quality/post-verification metric that does
    not require any experimental or DNS data.

    Operands:
    \table
      Operand        | Type                | Location
      input          | vol\<Type\>Field(s) | \<time\>/\<inpField\>
      output file    | -                   | -
      output field   | volScalarField      | \<time\>/\<outField\>
    \endtable

    References:
    \verbatim
        Governing equation (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531

        Governing equations (tag:CKJ):
            Celik, I., Klein, M., & Janicka, J. (2009).
            Assessment measures for engineering LES applications.
            Journal of fluids engineering, 131(3).
            DOI:10.1115/1.3059703
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    resolutionIndex1
    {
        // Mandatory entries
        type            resolutionIndex;
        libs            (fieldFunctionObjects);
        model           <word>;

        // Conditional entries

            // Option-1: when model == PopeIndex

            // Option-2: when model == CelikNuIndex

            // Option-3: when model == CelikEtaIndex

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: resolutionIndex         | word | yes  | -
      libs         | Library name: fieldFunctionObjects | word | yes  | -
      model        | Name of the resolutionIndex model  | word | yes  | -
    \endtable

    Options for the \c model entry:
    \verbatim
      PopeIndex    | Resolution method proposed by Pope (2000)
      CelikNuIndex | Resolution method proposed by Celik et al. (2009)
      CelikEtaIndex | Resolution method proposed by Celik et al. (2009)
    \endverbatim

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink

Note
  - The \c resolutionIndex throws runtime error
    when the turbulence model is RANS based.
  - Resolution-estimator models strictly assume that
    the flow is in a fully turbulent regime.

SourceFiles
    resolutionIndex.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_resolutionIndex_H
#define Foam_functionObjects_resolutionIndex_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class resolutionIndexModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class resolutionIndex Declaration
\*---------------------------------------------------------------------------*/

class resolutionIndex
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Resolution index model
        autoPtr<resolutionIndexModel> resolutionIndexModelPtr_;


public:

    //- Runtime type information
    TypeName("resolutionIndex");


    // Constructors

        //- Construct from Time and dictionary
        resolutionIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        resolutionIndex(const resolutionIndex&) = delete;

        //- No copy assignment
        void operator=(const resolutionIndex&) = delete;


    // Destructor
    virtual ~resolutionIndex();


    // Member Functions

        //- Read the resolutionIndex data
        virtual bool read(const dictionary& dict);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
