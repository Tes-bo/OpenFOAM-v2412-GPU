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
    Foam::functionObjects::runTimePostPro::scalarBar

Description
    Handling of scalar bar setup.

    Dictionary controls - scalar-bar entries (when \c visible is true)
    \table
        Property    | Description                           | Required | Default
        visible     | Display scalar bar                    | no  | true
        title       | The title for the scalar bar          | yes |
        position    | Viewport position (x y) of scalar bar | yes |
        vertical    | Vertical scalar bar                   | no  | true
        size        | Viewport size (x y) of scalar bar     | no  | auto
        fontSize    | Label size                            | no  | 0 == auto
        titleSize   | Title font size                       | no  | 0 == auto
        labelFormat | Label format string (eg, "%f")        | no  | "%f"
        numberOfLabels | Total number of labels             | no  | 5
        bold        | Title in bold                         | no  | yes
        italic      | Title in italic font                  | no  | no
        shadow      | Title with shadowont                  | no  | no
        titleHack   | Alternative placement strategy        | no  | no
    \endtable

Note
    The \c titleHack option is a leftover from much older VTK versions
    that had poor handling of the scalar bar text. When this is active,
    the normal scalar bar title is suppressed and placed manually
    (with a predefined ratio of the base font-size).

    When the titleHack is off, the label and title font sizes are
    automatically defined by the scalar bar size but can be overridden with
    provided values.

    The effect of fontSize and titleSize is unfortunately not consistent
    between titleHack on/off.

SourceFiles
    scalarBar.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_scalarBar_H
#define functionObjects_runTimePostPro_scalarBar_H

#include "dictionary.H"
#include "Tuple2.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkLookupTable;
class vtkRenderer;


namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class scalarBar Declaration
\*---------------------------------------------------------------------------*/

class scalarBar
{
    //- The default font size when titleHack is on
    static constexpr const label defaultFontSize = 18;

    //- The auto title/label ratio when titleHack is on
    static constexpr const float defaultTitleSizeFactor = 1.5;

    bool visible_;
    bool vertical_;
    bool bold_;
    bool italic_;
    bool shadow_;
    bool titleHack_;

    Tuple2<scalar, scalar> position_;

    Tuple2<scalar, scalar> size_;

    string title_;

    label fontSize_;

    label titleSize_;

    label nLabels_;

    string labelFormat_;


public:

    //- Default construct with sensible defaults
    scalarBar();

    //- Reset to sensible defaults
    void clear();

    //- Make non-visible
    void hide();

    //- Read dictionary settings
    void read(const dictionary& dict);

    //- Add colour bar, when visible.
    bool add
    (
        const vector& textColour,
        vtkRenderer* renderer,
        vtkLookupTable* lut
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
