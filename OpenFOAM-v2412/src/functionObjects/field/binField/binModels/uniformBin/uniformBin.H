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
    Foam::binModels::uniformBin

Description
    Calculates binned data in multiple segments according to
    a specified Cartesian or cylindrical coordinate system.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    binField1
    {
        // Other binField entries
        ...

        // Mandatory entries
        binModel          uniformBin;

        binData
        {
            // Mandatory entries
            nBin          <Vector<label>>;

            // Optional entries
            cumulative    <bool>;
            minMax
            {
                e1        (<scalar> <scalar>); // (min max);
                e2        (<scalar> <scalar>);
                e3        (<scalar> <scalar>);
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd | Deflt
      binModel  | Type name: uniformBin                 | word   | yes  | -
      binData   | Entries of the chosen bin model       | dict   | yes  | -
      nBin      | Numbers of bins in specified directions | Vector\<label\> <!--
                -->                                              | yes  | -
      cumulative | Flag to bin data accumulated with increasing distance   <!--
                --> in binning direction                | bool   | no   | false
      minMax     | Min-max bounds in binning directions with respect to    <!--
                --> the coordinateSystem's origin       | dict   | no   | -
    \endtable

Note
  - The order of bin numbering is (e1, e2, e3), where the first
    varies the fastest. For example, for a cylindrical bin with
    \f$ nBin = (radial, azimuth, height) = (2, 4, 2) \f$, the bin indices
    may look like as follows - note the counterclockwise increments:
    \verbatim
                           |-------------------|
                           | 12 |         | 14 |
                                | 11 | 13 |
                                | 9  | 15 |
                           | 10 |         | 16 |
                           |-------------------|
                          /    /         /    /
                         /    /         /    /
                        |-------------------|
                        | 4  |         | 6  |
                             | 3  | 5  |
                             | 1  | 7  |
                        | 2  |         | 8  |
                        |-------------------|
    \endverbatim

SourceFiles
    uniformBin.C
    uniformBinTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_binModels_uniformBin_H
#define Foam_binModels_uniformBin_H

#include "binModel.H"
#include "writeFile.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{

/*---------------------------------------------------------------------------*\
                   Class uniformBin Declaration
\*---------------------------------------------------------------------------*/

class uniformBin
:
    public binModel
{
protected:

    // Protected Data

        //- Numbers of bins in binning directions
        Vector<label> nBins_;

        //- Equidistant bin widths in binning directions
        vector binWidth_;

        //- The geometric min/max bounds for the bins
        MinMax<vector> binLimits_;

        //- Face index to bin index addressing
        labelList faceToBin_;

        //- Cell index to bin index addressing
        labelList cellToBin_;


    // Protected Member Functions

        //- Write header for an binned-data file
        template<class Type>
        void writeFileHeader(OFstream& os) const;

        //- Initialise bin properties
        virtual void initialise();

        //- Return list of bin indices corresponding to positions given by d
        virtual labelList binAddr(const vectorField& d) const;

        //- Set/cache the bin addressing
        virtual void setBinsAddressing();

        //- Apply the binning to field fieldi
        template<class Type>
        bool processField(const label fieldi);


public:

    //- Runtime type information
    TypeName("uniformBin");


    // Constructors

        //- Construct from components
        uniformBin
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& outputPrefix
        );

        //- No copy construct
        uniformBin(const uniformBin&) = delete;

        //- No copy assignment
        void operator=(const uniformBin&) = delete;


    //- Destructor
    virtual ~uniformBin() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Apply bins
        virtual void apply();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformBinTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
