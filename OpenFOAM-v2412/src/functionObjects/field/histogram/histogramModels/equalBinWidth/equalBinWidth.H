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
    Foam::histogramModels::equalBinWidth

Description
    Histogram model which groups data into bins of equal width.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    histogram1
    {
        // Inherited entries
        ...

        // Mandatory entries
        nBins       <label>;

        // Optional entries
        min         <scalar>;
        max         <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd | Deflt
      nBins     | Number of histogram bins          | label  | yes  | -
      min       | Minimum value of histogram data   | scalar | no   | -
      max       | Maximum value of histogram data   | scalar | no   | -
    \endtable

Note
  - If \c max is not provided, \c histogram will use operand
    field's min and max as the bin extremes.
  - If \c max is provided without \c min,
    \c histogram will use 0 for the \c min.

SourceFiles
    equalBinWidth.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_histogramModels_equalBinWidth_H
#define Foam_histogramModels_equalBinWidth_H

#include "histogramModel.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace histogramModels
{

/*---------------------------------------------------------------------------*\
                        Class equalBinWidth Declaration
\*---------------------------------------------------------------------------*/

class equalBinWidth
:
    public histogramModel
{
    // Private Data

        //- Number of bins
        label nBins_;

        //- Min/max for histogram data
        scalarMinMax range_;


public:

    //- Runtime type information
    TypeName("equalBinWidth");


    // Constructors

        //- Construct from components
        equalBinWidth
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        equalBinWidth(const equalBinWidth&) = delete;

        //- No copy assignment
        void operator=(const equalBinWidth&) = delete;


    // Destructor
    virtual ~equalBinWidth() = default;


    // Member Functions

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Write data to stream and files
        virtual bool write(const bool log);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace histogramModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
