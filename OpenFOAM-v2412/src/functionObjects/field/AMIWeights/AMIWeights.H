/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::AMIWeights

Group
    grpFieldFunctionObjects

Description
    Computes the min/max/average weights of arbitrary mesh interface (AMI)
    patches, and optionally reports to a text file or writes VTK surfaces of
    the sum of the weights and mask fields for arbitrarily coupled mesh
    interface (ACMI) patches.

    Operands:
    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | dat            | $POST/\<file\>
      output field      | vtp            | $POST/\<AMINames\>_{src,tgt}.vtp
    \endtable

    where \c $POST=$FOAM_CASE/postProcessing/\<FO\>/\<time\>.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    AMIWeights1
    {
        // Mandatory entries (unmodifiable)
        type        AMIWeights;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        writeFields false;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: AMIWeights              | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      writeFields  | Write weights as VTK fields        | bool |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func AMIWeights
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::AMIWeights

SourceFiles
    AMIWeights.C
    AMIWeightsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_AMIWeights_H
#define functionObjects_AMIWeights_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class AMIWeights Declaration
\*---------------------------------------------------------------------------*/

class AMIWeights
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Flag to write AMI fields (as VTK files)
        bool writeFields_;

        //- List of AMI patch IDs
        labelList patchIDs_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Helper function to report patch information
        virtual void reportPatch(const cyclicAMIPolyPatch& pp);

        //- Write weight field
        void writeWeightField
        (
            const cyclicAMIPolyPatch& cpp,
            const scalarField& weightSum,
            const word& side
        ) const;

        //- Write weight fields if writeFields=true
        void writeWeightFields(const cyclicAMIPolyPatch& cpp) const;


public:

    //- Runtime type information
    TypeName("AMIWeights");


    // Constructors

        //- Construct from Time and dictionary
        AMIWeights
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        AMIWeights(const AMIWeights&) = delete;

        //- No copy assignment
        void operator=(const AMIWeights&) = delete;


    //- Destructor
    virtual ~AMIWeights() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the AMIWeights
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
