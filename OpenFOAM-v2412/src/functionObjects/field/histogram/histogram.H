/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::functionObjects::histogram

Group
    grpFieldFunctionObjects

Description
    Computes the volume-weighted histogram of an input \c volScalarField.

    Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | \<time\>/\<inpField\>
      output file    | dat            | postProcessing/\<FO\>/\<time\>/histogram
      output field   | -              | -
    \endtable

    The data written contains four columns (from left to right):
    - time
    - mid-point of histogram bin
    - histogram counts - number of samples in each bin
    - volume-weighted histogram values

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    histogram1
    {
        // Mandatory entries
        type            histogram;
        libs            (fieldFunctionObjects);
        field           <word>;
        model           <word>;

        // Conditional entries

            // Option-1: when model == equalBinWidth

            // Option-2: when model == unequalBinWidth

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: histogram               | word | yes  | -
      libs         | Library name: fieldFunctionObjects | word | yes  | -
      field        | Name of operand field              | word | yes  | -
      model        | Name of the histogram model        | word | yes  | -
    \endtable

    Options for the \c model entry:
    \verbatim
      equalBinWidth    | Use equal-bin width
      unequalBinWidth  | Use unequal-bin widths
    \endverbatim

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink

SourceFiles
    histogram.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_histogram_H
#define Foam_functionObjects_histogram_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class histogramModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class histogram Declaration
\*---------------------------------------------------------------------------*/

class histogram
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Histogram model
        autoPtr<histogramModel> histogramModelPtr_;


public:

    //- Runtime type information
    TypeName("histogram");


    // Constructors

        //- Construct from Time and dictionary
        histogram
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        histogram(const histogram&) = delete;

        //- No copy assignment
        void operator=(const histogram&) = delete;


    // Destructor
    virtual ~histogram() = default;


    // Member Functions

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Execute (effectively no-op)
        virtual bool execute();

        //- Calculate the histogram and write
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
