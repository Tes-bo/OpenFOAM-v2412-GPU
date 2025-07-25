/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::colourTable

Description
    Base class for generating a colour table from node points.

    Dictionary definition
    \table
        Property    | Description                           | Required | Default
        interpolate | rgb/hsv/diverging                     | no  | rgb
        table       | Node points for the colour table      | yes |
    \endtable

Predefined colour tables (in "etc/colourTables") include
"coolToWarm", "coldAndHot", "fire", "rainbow", "greyscale", "xray".

SourceFiles
    colourTable.C

\*---------------------------------------------------------------------------*/

#ifndef colourTable_H
#define colourTable_H

#include "Enum.H"
#include "List.H"
#include "Tuple2.H"
#include "vector.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class colourTable Declaration
\*---------------------------------------------------------------------------*/

class colourTable
{
public:

        //- Internal interpolation type
        enum interpolationType
        {
            RGB,
            HSV,
            DIVERGING
        };

        //- Enumeration of commonly used colour tables.
        //  The indices must match those in "etc/colourTables"
        enum predefinedType
        {
            COOL_WARM,        //!< "coolToWarm"
            COLD_HOT,         //!< "coldAndHot"
            FIRE,             //!< "fire" - ParaView "Black-Body Radiation"
            RAINBOW,          //!< "rainbow"
            GREYSCALE,        //!< greyscale - ParaView "Grayscale"
            XRAY              //!< "xray" - ParaView "X Ray"
        };


        //- Enumeration names for interpolationType
        static const Enum<interpolationType> interpolationTypeNames;

        //- Enumeration names for predefinedType
        static const Enum<predefinedType> predefinedNames;

        //- The data lookup type
        typedef Tuple2<scalar, vector> pair_type;


    // Lookup Colour Tables

        //- Look up pointer to colourTable by name, or nullptr on failure.
        static const colourTable* ptr(const word& tableName);

        //- Look up pointer to colourTable by type, or nullptr on failure.
        static const colourTable* ptr(const predefinedType tbl);

        //- Look up pointer to colourTable by name. Fatal on failure
        static const colourTable& ref(const word& tableName);

        //- Look up pointer to colourTable by type. Fatal on failure
        static const colourTable& ref(const predefinedType tbl);


private:

    // Private Static Data

        //- Predefined tables
        static HashPtrTable<colourTable> tables_;


    // Private Data

        //- The table control points
        List<pair_type> table_;

        //- Interpolator type
        interpolationType interp_;


    // Private Member Functions

        //- Construct from central "etc/colourTables" file.
        static void constructTables();


public:

    // Constructors

        //- Copy construct from table values
        explicit colourTable
        (
            const List<Tuple2<scalar, vector>>& values,
            const interpolationType interp = interpolationType::RGB
        );

        //- Copy construct from table values
        explicit colourTable
        (
            List<Tuple2<scalar, vector>>&& values,
            const interpolationType interp = interpolationType::RGB
        );

        //- Read construct from dictionary
        explicit colourTable
        (
            const dictionary& dict,
            const interpolationType interp = interpolationType::RGB
        );


    // Selectors

        //- Read as dictionary content
        static autoPtr<colourTable> New(Istream& is);


    //- Destructor
    virtual ~colourTable() = default;


    // Member Functions

    // Access

        //- Predefined tables
        static const HashPtrTable<colourTable>& tables();

        //- Return the colour at x. The input is clipped to 0-1 range.
        vector value(const scalar x) const;

        //- Return a discrete lookup table of colours
        List<Tuple2<scalar, vector>> table(const label nColours) const;


    // IO

        //- Write as dictionary format
        Ostream& writeDict(Ostream& os) const;
};


//- Write as dictionary format
Ostream& operator<<(Ostream& os, const colourTable& tbl);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
