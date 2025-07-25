/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::derivedFields

Group
    grpFieldFunctionObjects

Description
    Computes two predefined derived fields, i.e. \c rhoU, and \c pTotal, where
    the defined fields are hard-coded as follows:

    \vartable
      rhoU              | \f$ \rho \vec U \f$
      pTotal            | \f$ p + 1/2 \rho \, mag(\vec U)^2 \f$
    \endvartable

    Operands:
    \table
      Operand       | Type                    | Location
      input         | vol{Scalar,Vector}Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file   | -                       | -
      output field  | vol{Scalar,Vector}Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    derivedFields1
    {
        // Mandatory entries (unmodifiable)
        type            derivedFields;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        derived         (rhoU pTotal);

        // Optional entries (runtime modifiable)
        rhoRef          1.0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: derivedFields              | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      derived      | Names of operand fields (rhoU/pTotal) | word |  yes  | -
      rhoRef  | Reference density for incompressible flows | scalar | no  | 1.0
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::derivedFields

SourceFiles
    derivedFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_derivedFields_H
#define functionObjects_derivedFields_H

#include "fvMeshFunctionObject.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class derivedFields Declaration
\*---------------------------------------------------------------------------*/

class derivedFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        //- Options for the derived/calculated field type
        enum derivedType
        {
            NONE = 0,           //!< "none"
            MASS_FLUX,          //!< "rhoU"
            TOTAL_PRESSURE,     //!< "pTotal"
            UNKNOWN
        };

        //- Names for derivedType
        static const Enum<derivedType> knownNames;


protected:

    // Read from dictionary

        //- List of derived field (types) to create
        List<derivedType> derivedTypes_;

        //- Reference density (to convert from kinematic to static pressure)
        scalar rhoRef_;


    // Protected Member Functions

        //- Hard-coded derived field (rho * U)
        //  \return true if field did not previously exist
        bool add_rhoU(const word& derivedName);

        //- Hard-coded derived field (p + 1/2 * rho * U)
        //  \return true if field did not previously exist
        bool add_pTotal(const word& derivedName);


public:

    //- Run-time type information
    TypeName("derivedFields");


    // Constructors

        //- Construct from Time and dictionary
        derivedFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        derivedFields(const derivedFields&) = delete;

        //- No copy assignment
        void operator=(const derivedFields&) = delete;


    //- Destructor
    virtual ~derivedFields() = default;


    // Member Functions

        //- Remove (checkOut) derived fields from the object registry
        void removeDerivedFields();

        //- Read the data
        virtual bool read(const dictionary& dict);

        //- Calculate the derived fields
        virtual bool execute();

        //- Write derived fields
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& m);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
