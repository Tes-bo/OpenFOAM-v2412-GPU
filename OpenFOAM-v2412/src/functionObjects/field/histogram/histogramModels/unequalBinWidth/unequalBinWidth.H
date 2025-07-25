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
    Foam::histogramModels::unequalBinWidth

Description
    Histogram model which groups data into bins of unequal widths.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    histogram1
    {
        // Inherited entries
        ...

        // Mandatory entries
        ranges
        (
            // min       max
            (<scalar>    <scalar>)    // bin-0
            (<scalar>    <scalar>)    // bin-1
            ...
        );
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd | Deflt
      ranges    | Min-max values of histogram data <!--
                -->                 | List\<scalarMinMax\>  | no   | -
    \endtable

Note
  - All bins are half-open, that is [min, max).
  - Bins should be specified as consecutive, non-overlapping
    and adjacent intervals of the field variable. No warning
    or runtime error will be emitted, otherwise.

SourceFiles
    unequalBinWidth.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_histogramModels_unequalBinWidth_H
#define Foam_histogramModels_unequalBinWidth_H

#include "histogramModel.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace histogramModels
{

/*---------------------------------------------------------------------------*\
                        Class unequalBinWidth Declaration
\*---------------------------------------------------------------------------*/

class unequalBinWidth
:
    public histogramModel
{
    // Private Data

        //- Number of bins
        label nBins_;

        //- Lower and upper ranges of operand bins
        List<scalarMinMax> ranges_;


public:

    //- Runtime type information
    TypeName("unequalBinWidth");


    // Constructors

        //- Construct from components
        unequalBinWidth
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        unequalBinWidth(const unequalBinWidth&) = delete;

        //- No copy assignment
        void operator=(const unequalBinWidth&) = delete;


    // Destructor
    virtual ~unequalBinWidth() = default;


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
