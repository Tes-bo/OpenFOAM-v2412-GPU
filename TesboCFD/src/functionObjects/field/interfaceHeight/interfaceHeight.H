/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::interfaceHeight

Description
    This function object reports the height of the interface above a set of
    locations.

    For each location, it writes the vertical distance of the
    interface above both the location and the lowest boundary. It also writes
    the point on the interface from which these heights are computed. It uses
    an integral approach, so if there are multiple interfaces above or below a
    location then this method will generate average values.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file 1  | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/height
      output file 2  | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/position
      output field   | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    interfaceHeight1
    {
        // Mandatory entries (unmodifiable)
        type            interfaceHeight;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        locations       ((0 0 0) (10 0 0) (20 0 0));

        // Optional entries (runtime modifiable)
        alpha           alpha.water;
        liquid          true;
        direction       (1 0 0);
        interpolationScheme    cellPoint;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: interfaceHeight         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      locations    | Locations to report the height at  | vectorList | yes | -
      alpha        | Name of alpha field                | word |  no   | alpha
      liquid       | Flag if the alpha field that of the liquid | bool | no | true
      direction    | Direction of interface               | vector | no  | g
      interpolationScheme | Interpolation scheme        | word | no | cellPoint
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - Foam::functionObjects::logFile
    - ExtendedCodeGuide::functionObjects::field::interfaceHeight

SourceFiles
    interfaceHeight.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceHeight_H
#define interfaceHeight_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class interfaceHeight Declaration
\*---------------------------------------------------------------------------*/

class interfaceHeight
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Data

        //- Is the alpha field that of the liquid under the wave?
        bool liquid_;

        //- Name of the alpha field
        word alphaName_;

        //- Interpolation scheme
        word interpolationScheme_;

        //- Direction of interface motion
        vector direction_;

        //- List of locations to report the height at
        List<point> locations_;


    // Private Member Functions

        //- Output positions
        void writePositions();


    // Private Enumerations

        //- File enumeration
        enum class fileID
        {
            heightFile = 0,
            positionFile = 1
        };


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const fileID fid);


public:

    //- Runtime type information
    TypeName("interfaceHeight");


    // Constructors

        //- Construct from Time and dictionary
        interfaceHeight
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        interfaceHeight(const interfaceHeight&) = delete;

        //- No copy assignment
        void operator=(const interfaceHeight&) = delete;


    //- Destructor
    virtual ~interfaceHeight() = default;


    // Member Functions

        //- Inherit logFiles methods
        using logFiles::files;

        //- Return file corresponding to enumeration
        OFstream& files(const fileID fid)
        {
            return logFiles::files(label(fid));
        }


        //- Read
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
