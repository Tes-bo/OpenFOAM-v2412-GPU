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
    Foam::functionObjects::runTimePostPro::functionObjectBase

Description
    Base class for function object visualisation

    Dictionary controls
    \table
        Property     | Description                          | Required | Default
        functionObject | The data source                    | yes |
        clearObjects | Remove file after use                | no  | no
        liveObject   | Prefer simulation data source        | no  | true
    \endtable

    The "live" keyword indiates that from within the simulation (in memory)
    is preferred over data from disk (for example).

SourceFiles
    functionObjectBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_functionObjectBase_H
#define functionObjects_runTimePostPro_functionObjectBase_H

#include "fieldVisualisationBase.H"
#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectBase Declaration
\*---------------------------------------------------------------------------*/

class functionObjectBase
:
    public fieldVisualisationBase
{
protected:

    // Protected Data

        //- Reference to the state
        const stateFunctionObject& state_;

        //- The function object name that provides the source data
        word functionObjectName_;

        //- Flag to indicate "live" (simulation) data source should be used
        bool liveObject_;

        //- Flag to indicate that source data should be cleared after use
        bool clearObjects_;


    // Protected Member Functions

        //- Retrieve the filename used to create the scene object
        //- using information stored via the stateFunctionObject.
        //
        //  \param keyword is normally "file"
        //  \param subDictName is the sub-dictionary name, which is often
        //     the fieldName when one file is used per field.
        //
        //  Eg,
        //  \verbatim
        //      T
        //      {
        //          file    "<case>/relpath/T_object.vtk";
        //      }
        //      defaultCloud
        //      {
        //          file    "<case>/relpath/cloud_0001.vtp";
        //          fields  (T U);
        //      }
        //  \endverbatim
        fileName getFileName
        (
            const word& keyword,
            const word& subDictName
        ) const;

        //- Remove file used to create the scene object.
        //  Same naming semantics as getFileName.
        //
        //  \note does not change the stateFunctionObject
        bool removeFile(const word& keyword, const word& subDictName);


        //- No copy construct
        functionObjectBase(const functionObjectBase&) = delete;

        //- No copy assignment
        void operator=(const functionObjectBase&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        functionObjectBase
        (
            const stateFunctionObject& state,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~functionObjectBase() = default;


    // Member Functions

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
