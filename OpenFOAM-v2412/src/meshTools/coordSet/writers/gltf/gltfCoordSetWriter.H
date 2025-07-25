/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::coordSetWriters::gltfWriter

Description
    A coordSet(s) writer in glTF v2 format, which is particularly
    useful for writing track data.

    Two files are generated:
    - filename.bin  : a binary file containing all scene entities
    - filename.gltf : a JSON file that ties fields to the binary data

    The output can contain both geometry and fields, with additional support
    for colours using a user-supplied colour map, and animation of particle
    tracks.

    Controls are provided via the optional formatOptions dictionary.

    For non-particle track data:

    \verbatim
    formatOptions
    {
        gltf
        {
            // Apply colours flag (yes | no ) [optional]
            colour      yes;

            // List of options per field
            fieldInfo
            {
                p
                {
                    // Colour map [optional]
                    colourMap       <colourMap>;

                    // Colour map minimum and maximum limits [optional]
                    // Uses field min and max if not specified
                    min             0;
                    max             1;

                    // Alpha channel [optional] (<scalar>)
                    alpha           0.5;
                }
            }
        }
    }
    \verbatim

    For particle tracks:

    \verbatim
    formatOptions
    {
        gltf
        {
            // Apply colours flag (yes | no) [optional]
            colour      yes;

            // Animate tracks (yes | no) [optional]
            animate     yes;

            // Animation properties [optional]
            animationInfo
            {
                // Colour map [optional]
                colourMap       <colourMap>;

                // Colour [optional] (<vector> | uniform | field)
                colour          (1 0 0); // RGB in range [0-1]

                //colour          uniform;
                //colourValue     (1 0 0); // RGB in range [0-1]

                //colour          field;
                //colourField     d;

                // Colour map minimum and maximum limits [optional]
                // Note: for colour = field option
                // Uses field min and max if not specified
                min             0;
                max             1;

                // Alpha channel [optional] (<scalar>)
                alpha           0.5;
            }
        }
    }
    \endverbatim

Note
    When writing particle animations, the particle field and colour properties
    correspond to initial particle state (first data point) and cannot be
    animated (limitation of the file format).

    For more information on the specification see
    https://www.khronos.org/registry/glTF/

SourceFiles
    gltfCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_gltfWriter_H
#define Foam_coordSetWriters_gltfWriter_H

#include "coordSetWriter.H"
#include "colourTable.H"
#include "foamGltfFwd.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                         Class gltfWriter Declaration
\*---------------------------------------------------------------------------*/

class gltfWriter
:
    public coordSetWriter
{
public:

    // Enumerations

        //- Field option used for colours
        enum class fieldOption : char
        {
            NONE,       //!< Placeholder type (unnamed)
            UNIFORM,    //!< Uniform value
            FIELD       //!< Field value
        };


    //- Strings corresponding to the field options
    static const Enum<fieldOption> fieldOptionNames_;


private:

    // Private Data

        //- Backend output
        autoPtr<glTF::sceneWriter> writer_;

        //- Flag to animate - for particle tracks only
        bool animate_;

        //- Flag to add field colours
        bool colour_;

        //- Animation colour option
        fieldOption animateColourOption_;

        //- Animation colour field name
        word animateColourName_;

        //- Animation colour value
        vector animateColourValue_;

        //- Local field information
        const dictionary fieldInfoDict_;

        //- Animation information
        const dictionary animationDict_;

        //- The mesh indices in glTF scene
        labelList meshes_;


    // Private Member Functions

        //- Return the colour map name
        word getColourMap(const dictionary& dict) const;

        //- Return the colour table corresponding to the colour map
        const colourTable& getColourTable(const dictionary& dict) const;

        //- Return the named min/max field limits (from sub-dictionary)
        scalarMinMax getFieldLimits(const word& fieldName) const;

        //- Setup animation colour or field to search for
        void setupAnimationColour();

        //- Return the alpha field for mesh values
        tmp<scalarField> getAlphaField(const dictionary& dict) const;


        //- Templated write operation (static tracks)
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs
        );

        //- Write animated tracks
        template<class Type>
        fileName writeTemplate_animate
        (
            const word& fieldName,
            const UPtrList<const Field<Type>>& fieldPtrs
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,    //!< Name of field
            const Field<Type>& vals   //!< Local field values to write
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const List<Field<Type>>& fieldValues
        );


public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("gltf");


    // Constructors

        //- Default construct
        gltfWriter();

        //- Default construct with specified options
        explicit gltfWriter(const dictionary& options);

        //- Construct from components
        gltfWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        gltfWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~gltfWriter();


    // Member Functions

        //- Expected (characteristic) output file name - information only
        virtual fileName path() const;  // override

        //- Close and reset, clears backend.
        virtual void close(bool force = false);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const Time& t);  // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const instant& inst);  // override

        //- End time step. Clears existing backend.
        virtual void endTime();  // override


    // Write

        declareCoordSetWriterWriteMethod(label);
        declareCoordSetWriterWriteMethod(scalar);
        declareCoordSetWriterWriteMethod(vector);
        declareCoordSetWriterWriteMethod(sphericalTensor);
        declareCoordSetWriterWriteMethod(symmTensor);
        declareCoordSetWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSetWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
