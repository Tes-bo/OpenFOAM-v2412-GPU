/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "simpleRegIOobject.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

// Since dimensionSystems() can be reread we actually store a copy of
// the controlDict subDict (v.s. a reference to the subDict for e.g.
// dimensionedConstants)
dictionary* dimensionSystemsPtr_(nullptr);
HashTable<dimensionedScalar>* unitSetPtr_(nullptr);
dimensionSets* writeUnitSetPtr_(nullptr);

//! \cond file_scope

// Helper class to
//   register re-reader
//   deallocate demand-driven data
class addDimensionSetsToDebug
:
    public ::Foam::simpleRegIOobject
{
public:

    addDimensionSetsToDebug(const addDimensionSetsToDebug&) = delete;
    void operator=(const addDimensionSetsToDebug&) = delete;

    explicit addDimensionSetsToDebug(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addDimensionSetObject, name)
    {}

    virtual ~addDimensionSetsToDebug()
    {
        deleteDemandDrivenData(dimensionSystemsPtr_);
        deleteDemandDrivenData(unitSetPtr_);
        deleteDemandDrivenData(writeUnitSetPtr_);
    }

    virtual void readData(Foam::Istream& is)
    {
        deleteDemandDrivenData(dimensionSystemsPtr_);
        deleteDemandDrivenData(unitSetPtr_);
        deleteDemandDrivenData(writeUnitSetPtr_);
        dimensionSystemsPtr_ = new dictionary(is);
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << dimensionSystems();
    }
};

addDimensionSetsToDebug addDimensionSetsToDebug_("DimensionSets");

//! \endcond


dictionary& dimensionSystems()
{
    if (!dimensionSystemsPtr_)
    {
        dictionary* cachedPtr(nullptr);
        dimensionSystemsPtr_ = new dictionary
        (
            debug::switchSet
            (
                "DimensionSets",
                cachedPtr
            )
        );
    }
    return *dimensionSystemsPtr_;
}


const HashTable<dimensionedScalar>& unitSet()
{
    if (!unitSetPtr_)
    {
        const dictionary& dict = dimensionSystems();

        if (!dict.found("unitSet"))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find unitSet in dictionary " << dict.name()
                << exit(FatalIOError);
        }

        const word unitSetCoeffs(dict.get<word>("unitSet") + "Coeffs");

        const dictionary* unitDictPtr = dict.findDict(unitSetCoeffs);

        if (!unitDictPtr)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find " << unitSetCoeffs << " in dictionary "
                << dict.name() << nl
                << exit(FatalIOError);
        }

        const dictionary& unitDict = *unitDictPtr;

        unitSetPtr_ = new HashTable<dimensionedScalar>(2*unitDict.size());

        wordList writeUnitNames;

        for (const entry& dEntry : unitDict)
        {
            if ("writeUnits" == dEntry.keyword())
            {
                dEntry.readEntry(writeUnitNames);
            }
            else
            {
                dimensionedScalar dt;
                dt.read(dEntry.stream(), unitDict);

                bool ok = unitSetPtr_->insert(dEntry.keyword(), dt);
                if (!ok)
                {
                    FatalIOErrorInFunction(dict)
                        << "Duplicate unit " << dEntry.keyword()
                        << " in DimensionSets dictionary"
                        << exit(FatalIOError);
                }
            }
        }

        if (writeUnitNames.size() != 0 && writeUnitNames.size() != 7)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find entry \"writeUnits\" in " << unitDict.name()
                << " or it is not a wordList of size 7"
                << exit(FatalIOError);
        }

        writeUnitSetPtr_ = new dimensionSets(*unitSetPtr_, writeUnitNames);
    }

    return *unitSetPtr_;
}


const dimensionSets& writeUnitSet()
{
    if (!writeUnitSetPtr_)
    {
        (void)unitSet();
    }
    return *writeUnitSetPtr_;
}


const dimensionSet dimless;

const dimensionSet dimMass(1, 0, 0, 0, 0, 0, 0);
const dimensionSet dimLength(0, 1, 0, 0, 0, 0, 0);
const dimensionSet dimTime(0, 0, 1, 0, 0, 0, 0);
const dimensionSet dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet dimMoles(0, 0, 0, 0, 1, 0, 0);
const dimensionSet dimCurrent(0, 0, 0, 0, 0, 1, 0);
const dimensionSet dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const dimensionSet dimArea(sqr(dimLength));
const dimensionSet dimVolume(pow3(dimLength));
const dimensionSet dimVol(dimVolume);

const dimensionSet dimVelocity(dimLength/dimTime);
const dimensionSet dimAcceleration(dimVelocity/dimTime);

const dimensionSet dimDensity(dimMass/dimVolume);
const dimensionSet dimForce(dimMass*dimAcceleration);
const dimensionSet dimEnergy(dimForce*dimLength);
const dimensionSet dimPower(dimEnergy/dimTime);

const dimensionSet dimPressure(dimForce/dimArea);
const dimensionSet dimCompressibility(dimDensity/dimPressure);
const dimensionSet dimGasConstant(dimEnergy/dimMass/dimTemperature);
const dimensionSet dimSpecificHeatCapacity(dimGasConstant);
const dimensionSet dimViscosity(dimArea/dimTime);
const dimensionSet dimDynamicViscosity(dimDensity*dimViscosity);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSets::dimensionSets
(
    const HashTable<dimensionedScalar>& units,
    const wordList& unitNames
)
:
    units_(unitNames.size()),
    conversion_(unitNames.size()),
    conversionPivots_(unitNames.size()),
    valid_(false)
{
    forAll(unitNames, i)
    {
        units_.set
        (
            i,
            new dimensionedScalar
            (
                units[unitNames[i]]
            )
        );
    }

    if (unitNames.size() == 7)
    {
        valid_ = true;

        // Determine conversion from basic units to write units
        for (label rowI = 0; rowI < conversion_.m(); rowI++)
        {
            scalar* row = conversion_[rowI];

            for (label columnI = 0; columnI < conversion_.n(); columnI++)
            {
                const dimensionedScalar& dSet = units_[columnI];
                row[columnI] = dSet.dimensions()[rowI];
            }
        }
        conversionPivots_.setSize(conversion_.m());
        LUDecompose(conversion_, conversionPivots_);
    }
}


void Foam::dimensionSets::coefficients(scalarField& exponents) const
{
    LUBacksubstitute(conversion_, conversionPivots_, exponents);
}


// ************************************************************************* //
