/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::text

Description
    Define text element for runTimePostProcessing

    Example of text object specification:
    \verbatim
    text1
    {
        string      "text to display";
        position    (0.1 0.05);
        size        18;
        // halign      left;  // (left | centre | right)
        bold        yes;
        // Optional entry

        shadow      false;

        // Optionally override default colour
        // colour    (0 1 1);

        // Optional entry
        timeStamp   yes; // Append solution time to string
    }
    \endverbatim

    Dictionary controls
    \table
        Property    | Description                          | Required | Default
        string      | Text to display                      | yes |
        position    | The (x y) viewport position          | yes |
        positions   | Multiple (x y) viewport positions    | no  |
        size        | The font size in points              | yes |
        halign      | Text justification (left/centre/ right) | no | left
        bold        | Use bold font                        | yes |
        italic      | Use italic font                      | no  | false
        shadow      | Add text shadow                      | no  | false
        colour      | Override default text colour         | no  |
        timeStamp   | Append solution timeName to string   | no  | false
    \endtable

    Inherited controls
    \table
        Property    | Description                          | Required | Default
        visible     | Display the object                   | no  | true
        opacity     | Object opacity                       | no  | 1.0
    \endtable

Note
    The string text is expanded on input using stringOps::inplaceExpand()
    to expand dictionary and environment variable entries.

SourceFiles
    text.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_text_H
#define functionObjects_runTimePostPro_text_H

#include "geometryBase.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkRenderer;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                            Class text Declaration
\*---------------------------------------------------------------------------*/

class text
:
    public geometryBase
{
public:

    // Public Enumerations

        //- Horizontal alignment type. Values to match VTK definitions
        enum halignType
        {
            LEFT   = 0,   //!< Left-justified text - default ("left")
            CENTER = 1,   //!< Centred text ("center", "centre")
            RIGHT  = 2    //!< Right-justified text  ("right")
        };

        //- Horizontal alignment names (includes "center" and "centre")
        static const Enum<halignType> halignTypeNames;


protected:

    // Protected Data

        //- Text
        string string_;

        //- Position(s)
        List<Tuple2<scalar, scalar>> positions_;

        //- Font size
        scalar size_;

        //- Text colour
        autoPtr<Function1<vector>> colour_;

        //- Horizontal alignment
        halignType halign_;

        //- Bold flag
        bool bold_;

        //- Italic flag
        bool italic_;

        //- Shadow flag
        bool shadow_;

        //- Time stamp flag
        bool timeStamp_;


    // Protected Member Functions

        //- No copy construct
        text(const text&) = delete;

        //- No copy assignment
        void operator=(const text&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        text
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~text();


    // Member Functions

        //- Add text to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- Clear files used to create the object(s) - no-op
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
