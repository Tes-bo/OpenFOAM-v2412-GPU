/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::functionObjects::DMD

Group
    grpFieldFunctionObjects

Description
    Computes a dynamic mode decomposition model on a specified field.

    Dynamic mode decomposition (i.e. DMD) is a data-driven
    dimensionality reduction method. DMD is being used as a mathematical
    processing tool to reveal dominant modes out of a given field (or dataset)
    each of which is associated with a constant frequency and decay rate,
    so that dynamic features of a given flow may become interpretable,
    tractable, and even reproducible without computing simulations.
    DMD only relies on input data, therefore it is an equation-free approach.

    References:
    \verbatim
        DMD characteristics:
            Brunton S. L. (2018).
            Dynamic mode decomposition overview.
            Seattle, Washington: University of Washington.
            youtu.be/sQvrK8AGCAo (Retrieved:24-04-20)
    \endverbatim

    Operands:
    \table
      Operand      | Type                 | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
                   -->                    | \<time\>/\<inputField\>(s)
      output file  | dat           | postProcessing/\<FO\>/\<time\>/\<file\>(s)
      output field | volVectorField(s)    | \<time\>/\<outputField\>(s)
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

    Output fields:
    \verbatim
      modeRe_<modeIndex>_<field>_<FO> | Real part of a mode field
      modeIm_<modeIndex>_<field>_<FO> | Imaginary part of a mode field
    \endverbatim

    Output files:
    \verbatim
      dynamics_<field>.dat            | Dynamics data for each mode
      filtered_dynamics_<field>.dat   | Filtered dynamics data for each mode
    \endverbatim

    wherein for each mode, the following quantities are output into files:
    \vartable
      freq        | Frequency
      mag         | Magnitude
      ampRe       | Amplitude coefficients (real part)
      ampIm       | Amplitude coefficients (imaginary part)
      evalRe      | Eigenvalue (real part)
      evalIm      | Eigenvalue (imaginary part)
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    DMD1
    {
        // Mandatory entries
        type                DMD;
        libs                (fieldFunctionObjects);
        DMDModel            <DMDModel>;
        field               <word>;

        // Optional entries

            // Option-1
            patch               <word>;

            // Option-2
            patches             (<wordRes>);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: DMD                     | word | yes  | -
      libs         | Library name: fieldFunctionObjects | word | yes  | -
      DMDModel     | Name of specified DMD model        | word | yes  | -
      field        | Name of operand field              | word | yes  | -
      patch        | Name of operand patch              | word | no   | null
      patches      | Names of operand patches        | wordRes | no   | null
    \endtable

    Options for the \c DMDModel entry:
    \verbatim
      STDMD    | Streaming total dynamic mode decomposition
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -fields '(U p)' -time '10:'
    \endverbatim

Note
  - Warning: DMD is an active research area at the time of writing;
    therefore, there could be cases whereat oddities can be seen.

See also
  - Foam::functionObjects::fvMeshFunctionObject
  - Foam::functionObjects::writeFile
  - Foam::DMDModel
  - Foam::DMDModels::STDMD

SourceFiles
    DMD.C
    DMDTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DMD_H
#define functionObjects_DMD_H

#include "fvMeshFunctionObject.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class DMDModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                              Class DMD Declaration
\*---------------------------------------------------------------------------*/

class DMD
:
    public fvMeshFunctionObject
{
    typedef RectangularMatrix<scalar> RMatrix;

    // Private Data

        //- Dynamic mode decomposition model
        autoPtr<DMDModel> DMDModelPtr_;

        //- Augmented snapshot matrix (effectively a column vector)
        //  Upper half = current-time snapshot slot
        //  Lower half = previous-time snapshot slot
        //  A snapshot is an input dataset to be processed per execution step
        //  A single snapshot is usually referred to as the snapshot of a single
        //  time step, an augmented snapshot is constructed of two snapshots
        RMatrix z_;

        //- Names of operand patches
        const wordRes patches_;

        //- Name of operand field
        const word fieldName_;

        //- Number of elements in a snapshot
        label nSnap_;

        //- Current execution-step index of DMD,
        //- not necessarily that of the simulation
        label step_;


    // Private Member Functions

        // Evaluation

            //- Initialise snapshot at the first-execution step
            //  Initialisation at the ctor or read level is not possible
            //  since the operand field is not available in the database
            void initialise();

            //- Create operand snapshot by using
            //- current-time and previous-time operand fields
            void snapshot();

            //- Get operand field based on its base type
            template<class Type>
            bool getSnapshot();

            //- Store operand field based on its geometric
            //- field type after few manipulations
            //  Move previous-time field into previous-time slot in snapshot
            //  copy new current-time field into current-time slot in snapshot
            template<class GeoFieldType>
            bool storeSnapshot();


        // Access

            //- Return number of components of the base type of a given field
            label nComponents(const word& fieldName) const;

            //- Get the number of components of the base type of a given field
            template<class Type>
            bool nComponents(const word& fieldName, label& nComps) const;


public:

    //- Runtime type information
    TypeName("DMD");


    // Constructors

        //- Construct from Time and dictionary
        DMD
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        DMD(const DMD&) = delete;

        //- No copy assignment
        void operator=(const DMD&) = delete;


    //- Destructor
    virtual ~DMD() = default;


    // Member Functions

        //- Read DMD settings
        virtual bool read(const dictionary& dict);

        //- Execute DMD
        virtual bool execute();

        //- Write DMD results
        virtual bool write();

        //- Write DMD results
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DMDTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
