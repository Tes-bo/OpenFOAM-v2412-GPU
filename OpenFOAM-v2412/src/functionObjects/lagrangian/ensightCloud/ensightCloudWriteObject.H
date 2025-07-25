/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::functionObjects::ensightCloudWriteObject

Group
    grpLagrangianFunctionObjects

Description
    This functionObject writes cloud(s) in ensight format

    Example of function object specification:
    \verbatim
    cloudWrite1
    {
        type            ensightCloud;
        libs            (lagrangianFunctionObjects);
        writeControl    writeTime;
        writeInterval   1;
        format          ascii;

        timeFormat      scientific;
        timePrecision   5;

        cloud           myCloud;
        fields          (T U rho);
        width           4;  // file-padding

        selection
        {
            stride
            {
                // every 10th parcelId
                action  add;
                source  stride;
                stride  10;
            }
            Umin
            {
                // Remove slow parcels
                action  subtract;
                source  field;
                field   U;
                accept  (less 1e-3);
            }
            diam
            {
                // Only particular diameter ranges
                action  subset;
                source  field;
                field   d;
                accept  (greater 1e-3) and (less 1e-3);
            }
        }
    }
    \endverbatim

    \heading Basic Usage
    \table
        Property     | Description                      | Required | Default
        type         | Type name: ensightCloud          | yes |
        clouds       | List of clouds (name or regex)   | no  |
        cloud        | Cloud name                       | no  |
        fields       | List of fields (name or regex)   | no  |
        selection    | Parcel selection control         | no  | empty-dict
    \endtable

    \heading Output Options
    \table
        Property     | Description                      | Required | Default
        format       | Format as ascii or binary        | no  | binary
        width        | Mask width for \c data/XXXX      | no  | 8
        directory    | The output directory name    | no | postProcessing/NAME
        overwrite    | Remove existing directory        | no  | false
        consecutive  | Consecutive output numbering     | no  | false
        width        | Padding width for file name      | no  | 8
        prune        | Suppress writing of empty clouds | no  | false
        timeFormat   | Time format (ensight case)       | no  | scientific
        timePrecision | Time precision (ensight case)   | no  | 5
        writeControl | Output control                   | recommended | timeStep
    \endtable

    The output filename and fields are added to the functionObjectProperties
    information. For the previous example specification:

    \verbatim
    cloudWrite1
    {
        myCloud
        {
            file    "<case>/simulation.case";
            fields  (T U rho);
        }
    }
    \endverbatim

Note
    The selection dictionary can be used for finer control of the parcel
    output. It contains a set of (add,subtract,subset,clear,invert)
    selection actions and sources.
    Omitting the selection dictionary is the same as specifying the
    conversion of all parcels (in the selected clouds).
    More syntax details are to be found in the corresponding
    Foam::Detail::parcelSelection class.

See also
    Foam::Detail::parcelSelection
    Foam::functionObjects::vtkCloud
    Foam::functionObjects::ensightWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    ensightCloudWriteObject.cxx
    ensightCloudWriteObjectImpl.cxx

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ensightCloudWriteObject_H
#define functionObjects_ensightCloudWriteObject_H

#include "fvMeshFunctionObject.H"
#include "ensightCase.H"
#include "globalIndex.H"
#include "parcelSelectionDetail.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class ensightCloudWriteObject Declaration
\*---------------------------------------------------------------------------*/

class ensightCloudWriteObject
:
    public fvMeshFunctionObject,
    public Foam::Detail::parcelSelection
{
    // Private Data

        //- Ensight output options
        ensightCase::options caseOpts_;

        //- Output directory
        fileName outputDir_;

        //- Consecutive output numbering
        bool consecutive_;

        //- Suppress writing of empty clouds
        bool pruneEmpty_;

        //- Apply output filter (for the current cloud)
        bool applyFilter_;

        //- Sizing of selected parcels (including any filtering)
        globalIndex procAddr_;

        //- Requested names of clouds to process
        wordRes selectClouds_;

        //- Subset of cloud fields to process
        wordRes selectFields_;

        //- Ensight case handler
        autoPtr<ensightCase> ensCase_;


    // Private Member Functions

        //- Ensight case handler
        ensightCase& ensCase() { return *ensCase_; }

        //- Write a cloud to disk (creates parent directory),
        //- and record on the cloud OutputProperties.
        //  \param file is the output file name, with extension.
        bool writeCloud(const word& cloudName);

        //- Write fields of IOField<Type>
        template<class Type>
        wordList writeFields
        (
            const word& cloudName,
            const objectRegistry& obrTmp
        );


        //- No copy construct
        ensightCloudWriteObject(const ensightCloudWriteObject&) = delete;

        //- No copy assignment
        void operator=(const ensightCloudWriteObject&) = delete;


public:

    //- Runtime type information
    TypeName("ensightCloud");


    // Constructors

        //- Construct from Time and dictionary
        ensightCloudWriteObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ensightCloudWriteObject() = default;


    // Member Functions

        //- Read the ensightCloud specification
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
