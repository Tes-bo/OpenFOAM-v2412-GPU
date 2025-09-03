/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::DESModelRegions

Group
    grpFieldFunctionObjects

Description
    Computes an indicator field for detached eddy simulation (DES) turbulence
    calculations, where the values of the indicator mean:

    \verbatim
      0 = Reynolds-averaged Navier-Stokes (RAS) regions
      1 = Large eddy simulation (LES) regions
    \endverbatim

    Operands:
    \table
      Operand         | Type              | Location
      input           | -                 | -
      output file     | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field    | volScalarField    | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    DESModelRegions1
    {
        // Mandatory entries (unmodifiable)
        type        DESModelRegions;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        result      DESField;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: DESModelRegions         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      result       | Name of DES indicator field        | word |  no   | \<FO\>
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    \c DESModelRegions function object can only be executed for DES simulations.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::DESModelRegions

SourceFiles
    DESModelRegions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DESModelRegions_H
#define functionObjects_DESModelRegions_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class DESModelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESModelRegions
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Name of DES indicator field
        word resultName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("DESModelRegions");


    // Constructors

        //- Construct from Time and dictionary
        DESModelRegions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        DESModelRegions(const DESModelRegions&) = delete;

        //- No copy assignment
        void operator=(const DESModelRegions&) = delete;


    //- Destructor
    virtual ~DESModelRegions() = default;


    // Member Functions

        //- Read the DESModelRegions data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the DESModelRegions and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
