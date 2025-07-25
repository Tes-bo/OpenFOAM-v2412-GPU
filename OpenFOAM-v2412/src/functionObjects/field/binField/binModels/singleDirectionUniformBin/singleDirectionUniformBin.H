/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::binModels::singleDirectionUniformBin

Description
    Calculates binned data in a specified direction.

    For example, a 10cm-long patch extending only in the x-direction
    can be binned into 5 bins in the same direction, so that
    local information can be output for each 2cm-long segment.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    binField1
    {
        // Other binField entries
        ...

        // Mandatory entries
        binModel          singleDirectionUniformBin;

        binData
        {
            // Mandatory entries
            nBin          <label>;
            direction     <vector>;

            // Optional entries
            cumulative    <bool>;
            min           <scalar>;
            max           <scalar>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd | Deflt
      binModel  | Type name: singleDirectionUniformBin  | word   | yes  | -
      binData   | Entries of the chosen bin model       | dict   | yes  | -
      nBin      | Number of bins in binning direction   | label  | yes  | -
      direction | Binning direction                     | vector | yes  | -
      cumulative | Flag to bin data accumulated with increasing distance <!--
                --> in binning direction                | bool   | no   | false
      min       | Min-bound in the binning direction with respect to   <!--
                --> the global coordinate system's origin | scalar | no | GREAT
      max       | Max-bound in the binning direction with respect to   <!--
                --> the global coordinate system's origin | scalar | no | GREAT
    \endtable

SourceFiles
    singleDirectionUniformBin.C
    singleDirectionUniformBinTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_binModels_singleDirectionUniformBin_H
#define Foam_binModels_singleDirectionUniformBin_H

#include "binModel.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{

/*---------------------------------------------------------------------------*\
                   Class singleDirectionUniformBin Declaration
\*---------------------------------------------------------------------------*/

class singleDirectionUniformBin
:
    public binModel
{
protected:

    // Protected Data

        //- Distance between bin divisions
        scalar binWidth_;

        //- The min/max bounds for the bins
        MinMax<scalar> binLimits_;

        //- Binning direction
        vector binDir_;


    // Protected Member Functions

        //- Write header for a binned-data file
        template<class Type>
        void writeFileHeader(OFstream& os) const;

        //- Initialise bin properties
        virtual void initialise();

        //- Apply the binning to field fieldi
        template<class Type>
        bool processField(const label fieldi);


public:

    //- Runtime type information
    TypeName("singleDirectionUniformBin");


    // Constructors

        //- Construct from components
        singleDirectionUniformBin
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& outputPrefix
        );

        //- No copy construct
        singleDirectionUniformBin(const singleDirectionUniformBin&) = delete;

        //- No copy assignment
        void operator=(const singleDirectionUniformBin&) = delete;


    //- Destructor
    virtual ~singleDirectionUniformBin() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Apply bins
        virtual void apply();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleDirectionUniformBinTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
