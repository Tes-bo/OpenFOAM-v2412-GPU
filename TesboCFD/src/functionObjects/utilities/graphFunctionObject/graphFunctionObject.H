/*---------------------------------------------------------------------------*\
  =========                |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::functionObjects::graphFunctionObject

Description
    Accumulates function object result values and renders into a graph in
    SVG format.

    Operands:
    \table
    Operand      | Type                                 | Location
    input        | Function object results              | Memory; <!--
    --> $FOAM_CASE/\<time\>/uniform/functionObjects/functionObjectProperties
    output file  | SVG                                  | <!--
    --> $FOAM_CASE/postProcessing/<functionObject>/<time>/<functionObject>.svg
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions to plot the
    residuals from the \c solverInfo function object:

    \verbatim
    residualGraph
    {
        // Mandatory entries
        type            graphFunctionObject;
        libs            (utilityFunctionObjects);

        functions
        {
            entry
            {
                // Mandatory entries
                object      <word>;
                entry       <word>;

                // Optional entries
                title       <string>;
                colour      <labelVector>;
                dashes      <labelList>;
            }
            line1
            {
                object      solverInfo1;
                entry       Ux_initial;
            }
            line2
            {
                object      solverInfo1;
                entry       p_initial;
            }
        }

        // Optional entries
        xMin            <scalar>;
        xMax            <scalar>;
        yMin            <scalar>;
        yMax            <scalar>;
        xlabel          <string>;  // "Iteration";
        ylabel          <string>;  // "log10(Initial residual)";
        width           <label>;
        height          <label>;
        strokeWidth     <label>;
        logScaleX       <bool>;
        logScaleY       <bool>;
        drawGrid        <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
    Property     | Description                        | Type | Reqd | Deflt
    type         | Type name: graphFunctionObject     | word |  yes  | -
    libs         | Library name: utilityFunctionObjects | word | yes | -
    functions    | Dictionary of lines to draw        | dictionary | yes | -
    width        | Output SVG width in pixel          | label| no    | 800
    height       | Output SVG height in pixel         | label| no    | 600
    xMin         | User defined minimum x axis limit  | scalar | no | calculated
    xMax         | User defined maximum x axis limit  | scalar | no | calculated
    yMin         | User defined minimum y axis limit  | scalar | no | calculated
    yMax         | User defined maximum y axis limit  | scalar | no | calculated
    xLabel       | X axis label                 | string | no  | Iteration/Time
    yLabel       | Y axis label                 | string | no  | Property
    strokeWidth  | Line stroke width in pixel         | label  | no  | 2
    logScaleX    | Use log scale for x axis           | bool   | no  | false
    logScaleY    | Use log scale for y axis           | bool   | no  | false
    drawGrid     | Draw background grid               | bool   | no  | true
    \endtable

    The inherited entries are elaborated in:
      - \link stateFunctionObject.H \endlink
      - \link writeFile.H \endlink

    Each line corresponds to the history of function object result values, e.g.

    \verbatim
    line1
    {
        object      solverInfo1;
        entry       Ux_initial;
        colour      (255 0 0);
        dashes      (4 1);
        title       Ux;
    }
    \endverbatim

    where the entries mean:
    \table
    Property     | Description                        | Type | Reqd  | Deflt
    object       | Function object name               | word |  yes  | -
    entry        | Function object result entry name  | word |  yes  | -
    colour       | Line colour                        | label vector | no | auto
    dashes       | Line dash array                    | label vector | no | auto
    title        | Title                              | string | no | dict name
    \endtable

See also
    - Foam::functionObjects::solverInfo

SourceFiles
    graphFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_graphFunctionObject_H
#define Foam_functionObjects_graphFunctionObject_H

#include "stateFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class graphFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class graphFunctionObject
:
    public stateFunctionObject,
    public writeFile
{
    // Private Data

        //- List of default curve colours
        static wordList defaultColours;

        //- Names of function objects
        List<word> objects_;

        //- Function object entries
        List<word> entries_;

        //- Line titles
        List<string> titles_;

        //- Line colours - hex string
        List<word> colours_;

        //- Line dash array
        List<string> dashes_;

        //- Times
        DynamicList<scalar> times_;

        //- Time vs. flattened values
        DynamicList<DynamicList<scalar>> values_;

        //- Mapping from object to column index in values_
        List<DynamicList<label>> objectToCol_;

        //- User-supplied minimum x value
        scalar xMin_;

        //- User-supplied maximum x value
        scalar xMax_;

        //- User-supplied minimum y value
        scalar yMin_;

        //- User-supplied maximum y value
        scalar yMax_;

        //- X axis label
        const string xlabel_;

        //- Y axis label
        const string ylabel_;

        //- Width in px
        const label width_;

        //- Height in px
        const label height_;

        //- Line width in px
        const label strokeWidth_;

        //- Flag to use log scale on x-axis
        bool logScaleX_;

        //- Flag to use log scale on y-axis
        bool logScaleY_;

        //- Draw background grid
        const bool drawGrid_;


    // Private Functions

        //- Get the result value from the function object
        //  \returns true if the value was found
        template<class Type>
        bool getValue(const label objecti, label& valuei);

        //- Set axis min, max, tick
        //  \returns number of ticks
        label setAxisProps
        (
            const bool logScale,
            scalar& xmin,
            scalar& xmax,
            scalar& xtick
        ) const;

        //- No copy construct
        graphFunctionObject(const graphFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const graphFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("graphFunctionObject");

    //- Construct from Time and dictionary
    graphFunctionObject
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~graphFunctionObject() = default;


    // Member Functions

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