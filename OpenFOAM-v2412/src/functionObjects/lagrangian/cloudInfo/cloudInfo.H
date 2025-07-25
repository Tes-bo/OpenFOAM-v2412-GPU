/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::functionObjects::cloudInfo

Group
    grpLagrangianFunctionObjects

Description
    Outputs Lagrangian cloud information to a file.

    The current outputs include:
    - total current number of parcels
    - total current mass of parcels

Usage
    Example of function object specification:
    \verbatim
    cloudInfo1
    {
        type        cloudInfo;
        libs        (lagrangianFunctionObjects);
        ...
        clouds
        (
            kinematicCloud1
            thermoCloud1
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                          | Required | Default
        type         | type name: cloudInfo                 | yes |
        clouds       | list of clouds names to process      | yes |
        selection    | Parcel selection control         | no  | empty-dict
        sampleOnExecute| Sample/report (on execute) without writing | no | false
    \endtable

    The output data of each cloud is written to a file named \<cloudName\>.dat

See also
    Foam::functionObject
    Foam::functionObjects::regionFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    cloudInfo.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_cloudInfo_H
#define Foam_functionObjects_cloudInfo_H

#include "regionFunctionObject.H"
#include "logFiles.H"
#include "parcelSelectionDetail.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class cloudInfo Declaration
\*---------------------------------------------------------------------------*/

class cloudInfo
:
    public functionObjects::regionFunctionObject,
    public functionObjects::logFiles,
    public Foam::Detail::parcelSelection
{
protected:

    // Data Types

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_ALL = 0xF
        };


    // Protected Data

        //- Additional verbosity
        bool verbose_;

        //- Perform sample actions on execute as well
        bool onExecute_;

        //- List of cloud names
        wordList cloudNames_;

        //- Output file per cloud
        PtrList<OFstream> filePtrs_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Perform operation report/write
        bool performAction(unsigned request);

        //- No copy construct
        cloudInfo(const cloudInfo&) = delete;

        //- No copy assignment
        void operator=(const cloudInfo&) = delete;


public:

    //- Runtime type information
    TypeName("cloudInfo");


    // Constructors

        //- Construct from Time and dictionary
        cloudInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~cloudInfo() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
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
