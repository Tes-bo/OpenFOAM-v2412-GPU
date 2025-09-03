/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::functionObjects::forceCoeffs

Group
    grpForcesFunctionObjects

Description
    Computes force and moment coefficients over a given
    list of patches, and optionally over given porous zones.
    The following coefficients can be selected and output:
    \verbatim
      Cd      | Drag coefficient
      Cs      | Side-force coefficient
      Cl      | Lift coefficient
      CmRoll  | Roll-moment coefficient
      CmPitch | Pitch-moment coefficient
      CmYaw   | Yaw-moment coefficient
    \endverbatim

    The force coefficients can also be optionally output
    in terms of their front and rear axle constituents:
    \verbatim
        Cd{f,r} = 0.5*Cd {+,-} CmRoll
        Cs{f,r} = 0.5*Cs {+,-} CmYaw
        Cl{f,r} = 0.5*Cl {+,-} CmPitch
    \endverbatim
    where \c f and \c r represent front and rear axles, respectively.

    Force and moment coefficients are output
    in their total and constituent components:
    - total force and moment coefficients
    - pressure contributions
    - viscous contributions
    - porous resistance contributions (optional)

    Force and moment coefficients can be computed and output in:
    - the global Cartesian coordinate system (default)
    - a user-defined Cartesian coordinate system

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat            | postProcessing/\<FO\>/\<time\>/\<file\>s
      output field  | volVectorField | \<time\>/\<outField\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      coefficient.dat  | Integrated coefficients over all patches
    \endverbatim

    where \c \<outField\>s:
    \verbatim
      <namePrefix>:forceCoeff   | Force coefficient field
      <namePrefix>:momentCoeff  | Moment coefficient field
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    forceCoeffs1
    {
        // Mandatory entries
        type                forceCoeffs;
        libs                (forces);
        patches             (<wordRes>); // (wall1 "(wall2|wall3)");
        magUInf             <scalar>;
        lRef                <scalar>;
        Aref                <scalar>;

        // Optional entries
        coefficients        (<wordHashSet>);
        directForceDensity  <bool>;
        porosity            <bool>;
        writeFields         <bool>;
        useNamePrefix       <bool>;

        // Conditional mandatory entries

            // Cartesian coordinate system specification when evaluating
            // force and moment coefficients, either of the below

            // Define the centre of rotation
            // with implicit directions e1=(1 0 0) and e3=(0 0 1)
            CofR            (0 0 0); // Centre of rotation

            // Define local coordinate system by origin + axes
            origin          (0 0 0);
            e1              (1 0 0);
            e3              (0 0 1); // (e1, e2) or (e2, e3) or (e3, e1)

            // General coordinate system specification (always cartesian)
            coordinateSystem
            {
                origin      (0 0 0);
                rotation
                {
                    type    ...;
                    ...
                }
            }

        // Conditional optional entries

            // if directForceDensity == true
            fD              <word>;

            // if directForceDensity == false
            p               <word>;
            U               <word>;
            rho             <word>;
            rhoInf          <scalar>; // redundant for incompressible-flow cases
            pRef            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                     | Type     | Reqd | Deflt
      type       | Type name: forceCoeffs          | word     | yes  | -
      libs       | Library name: forces            | word     | yes  | -
      patches    | Names of operand patches        | wordRes  | yes  | -
      coefficients | Names of operand coefficients | wordHashSet | no | -
      magUInf    | Reference velocity magnitude       | scalar | yes | -
      lRef       | Reference length scale for moment calculations   <!--
                 -->                                  | scalar | yes | -
      Aref       | Reference area                     | scalar | yes | -
      directForceDensity | Flag to directly supply force density    <!--
                 -->                                  | bool   | no  | false
      porosity   | Flag to include porosity contributions | bool | no | false
      writeFields | Flag to write force and moment fields | bool | no | false
      useNamePrefix | Flag to include prefix for field names | bool | no | false
      CofR    | Centre of rotation                    | vector | cndtnl | -
      origin  | Origin of coordinate system           | vector | cndtnl | -
      e3      | e3 coordinate axis                    | vector | cndtnl | -
      e1      | e1 coordinate axis                    | vector | cndtnl | -
      coordinateSystem | Coordinate system specifier  | dictionary | cndtnl | -
      fD      | Name of force density field  | word   | cndtnl-no | fD
      p       | Name of pressure field       | word   | cndtnl-no | p
      U       | Name of velocity field       | word   | cndtnl-no | U
      rho     | Name of density field        | word   | cndtnl-no | rho
      rhoInf  | Value of reference density   | scalar | cndtnl-yes | -
      pRef    | Value of reference pressure  | scalar | cndtnl-no | 0
    \endtable

    Options for the \c coefficients entry:
    \verbatim
      Cd      | Drag coefficient
      Cs      | Side-force coefficient
      Cl      | Lift coefficient
      CmRoll  | Roll-moment coefficient
      CmPitch | Pitch-moment coefficient
      CmYaw   | Yaw-moment coefficient
      Cd(f)   | Front-axle drag coefficient
      Cs(f)   | Front-axle side-force coefficient
      Cl(f)   | Front-axle lift coefficient
      Cd(r)   | Rear-axle drag coefficient
      Cs(r)   | Rear-axle side-force coefficient
      Cl(r)   | Rear-axle lift coefficient
    \endverbatim

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink
      - \link forces.H \endlink

Note
  - \c rhoInf is always redundant for incompressible computations.
    That is, \c rhoInf is always equal to 1 in incompressible
    computations no matter which input value is assigned to \c rhoInf.
    The value of \c rhoInf is only used for compressible computations.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force and moment files and fields.
  - Input for force/moment coefficient directions
    require an origin and two orthogonal directions where
    the remaining orthogonal direction is automatically computed.
  - The default direction relations are shown below:

    \table
      Property     | Description           | Alias | Direction
      dragDir      | Drag direction        | e1    | (1 0 0)
      sideDir      | Side force direction  | e2    | (0 1 0)
      liftDir      | Lift direction        | e3    | (0 0 1)
      rollAxis     | Roll axis             | e1    | (1 0 0)
      pitchAxis    | Pitch axis            | e2    | (0 1 0)
      yawAxis      | Yaw axis              | e3    | (0 0 1)
    \endtable

  - If a \c coordinateSystem entry exists, it is taken in favour of \c CofR.

SourceFiles
    forceCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_forceCoeffs_H
#define Foam_functionObjects_forceCoeffs_H

#include "forces.H"
#include "HashSet.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffs Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffs
:
    public forces
{
public:

    // Container for force and moment components
    class forceComponents
    {
        // Private Data

        //- Pressure force/moment component
        vector pressure_;

        //- Viscous force/moment component
        vector viscous_;

        //- Internal force/moment component
        vector internal_;


    public:

        // Constructors

        //- Default construct
        forceComponents()
        :
            pressure_(Zero),
            viscous_(Zero),
            internal_(Zero)
        {}


        // Member Functions

        //- Return the sum of the components
        const vector total() const noexcept
        {
            return pressure_ + viscous_ + internal_;
        }

        //- Reset the components to zero
        void reset() noexcept
        {
            pressure_ = Zero;
            viscous_ = Zero;
            internal_ = Zero;
        }

        //- Reset the components to given values
        void reset
        (
            const vector& pressure,
            const vector& viscous,
            const vector& internal
        ) noexcept
        {
            pressure_ = pressure;
            viscous_ = viscous;
            internal_ = internal;
        }


        // Member Operators

        //- Return components in a given direction
        vector operator[](const label cmpt) const
        {
            return vector(pressure_[cmpt], viscous_[cmpt], internal_[cmpt]);
        }
    };


    //- Coefficients description
    struct coeffDesc
    {
        enum splitType
        {
            stFront,
            stRear,
            stNone
        };

        string desc_;
        word name_;
        label c0_;
        label c1_;
        bool active_;
        splitType splitType_;

        // Constructors

        //- Construct from components
        coeffDesc
        (
            const string& description,
            const word& name,
            const label c0,
            const label c1 = -1
        )
        :
            desc_(description),
            name_(name),
            c0_(c0),
            c1_(c1),
            active_(false),
            splitType_(stNone)
        {}

        //- Return name with front-name appendix
        const word frontName() const noexcept
        {
            return name_ + "(f)";
        }

        //- Return name with rear-name appendix
        const word rearName() const noexcept
        {
            return name_ + "(r)";
        }

        //- Print description
        void print(Ostream& os) const
        {
            os  << desc_.c_str() << ": " << name_;
        }

        //- Return force/moment components based on the specified split type
        vector value(const forceComponents& f, const forceComponents& m) const
        {
            if (c1_ == -1)
            {
                return m[c0_];
            }
            else
            {
                switch (splitType_)
                {
                    case stFront:
                    {
                        return 0.5*f[c0_] + m[c1_];
                    }
                    case stRear:
                    {
                        return 0.5*f[c0_] - m[c1_];
                    }
                    case stNone:
                    {
                        return f[c0_];
                    }
                }
            }

            FatalErrorInFunction
                << "Cannot determine value"
                << abort(FatalError);

            return vector::zero;
        }

        //- Return front-axle coefficient description
        coeffDesc front() const
        {
            coeffDesc coeff(*this);
            coeff.name_ = coeff.frontName();
            coeff.desc_ = coeff.desc_ + " front";
            coeff.splitType_ = stFront;
            return coeff;
        }

        //- Return rear-axle coefficient description
        coeffDesc rear() const
        {
            coeffDesc coeff(*this);
            coeff.name_ = coeff.rearName();
            coeff.desc_ = coeff.desc_ + " rear";
            coeff.splitType_ = stRear;
            return coeff;
        }
    };


private:

    // Private Data

        //- Force components
        forceComponents Cf_;

        //- Moment components
        forceComponents Cm_;

        //- Operand force and moment coefficients
        HashTable<coeffDesc> coeffs_;


        // File streams

            //- File stream for integrated operand coefficients
            autoPtr<OFstream> coeffFilePtr_;


        // Reference scales

            //- Reference velocity magnitude [m/s]
            scalar magUInf_;

            //- Reference length scale [m]
            scalar lRef_;

            //- Reference area [m^2]
            scalar Aref_;


        //- Flag of initialisation (internal)
        bool initialised_;


protected:

    // Protected Member Functions

        //- Initialise containers and fields
        void initialise();

        //- Return access to the force coefficients field
        volVectorField& forceCoeff();

        //- Return access to the moment coefficients field
        volVectorField& momentCoeff();


        //- Reset containers and fields
        void reset();


    // Evaluation

        //- Return the operand coefficients
        HashTable<coeffDesc> selectCoeffs() const;

        //- Calculate force coefficients
        void calcForceCoeffs();

        //- Calculate moment coefficients
        void calcMomentCoeffs();

        //- Return integrated {total, pressure, viscous, porous} components
        List<scalar> integrateData(const List<Field<scalar>>& coeff) const;


    // I-O

        //- Create the integrated-coefficient file
        void createIntegratedDataFile();

        //- Write header to the integrated-coefficient file
        void writeIntegratedDataFileHeader
        (
            const word& header,
            OFstream& os
        ) const;

        //- Write integrated coefficients to the integrated-coefficient file
        void writeIntegratedDataFile();


public:

    //- Runtime type information
    TypeName("forceCoeffs");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffs
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        forceCoeffs(const forceCoeffs&) = delete;

        //- No copy assignment
        void operator=(const forceCoeffs&) = delete;


    //- Destructor
    virtual ~forceCoeffs() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();
};

} // End namespace functionObjects


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Output operator
inline Ostream& operator<<
(
    Ostream& os,
    const functionObjects::forceCoeffs::coeffDesc& coeff
)
{
    coeff.print(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
