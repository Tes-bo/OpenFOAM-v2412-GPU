/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "adiosReader.H"
#include "stringOps.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::unique_ptr<Foam::objectRegistry>
Foam::adiosFoam::adiosReader::getCloud
(
    const word& regName,
    const word& cloudName
)
{
    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    std::unique_ptr<objectRegistry> obrPtr
    (
        new objectRegistry
        (
            IOobject
            (
                cloudName,
                *dummyTimePtr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    auto& obr = *obrPtr;

    const fileName prefixCloud(adiosFoam::cloudPath(regName, cloudName));

    const std::map<std::string, adios2::Params>& vars = availableVars_;

    forAllConstIters(vars, iter)
    {
        const fileName varName(iter->first);

        if (varName.ends_with("_0"))
        {
            // The "_0" fields are addressed separately
            continue;
        }

        const auto nameCmpts = stringOps::split(varName, '/');

        if (nameCmpts.size() == 4 && varName.startsWith(prefixCloud))
        {
            // Or from pre-read?
            const word clsName(getStringAttribute(varName/"class"));

            const bool ok =
            (
                readIntoRegistry<label>(obr, varName, clsName)
             || readIntoRegistry<scalar>(obr, varName, clsName)
             || readIntoRegistry<vector>(obr, varName, clsName)
             || readIntoRegistry<sphericalTensor>(obr, varName, clsName)
             || readIntoRegistry<symmTensor>(obr, varName, clsName)
             || readIntoRegistry<tensor>(obr, varName, clsName)
            );

            if (!ok)
            {
                Info<<"Could not read into registry: " << varName << nl;
            }
            else
            {
                Info<<"Read into registry: " << varName << nl;
            }
        }
    }

    return obrPtr;
}


std::unique_ptr<Foam::objectRegistry>
Foam::adiosFoam::adiosReader::getCloud
(
    const adiosFoam::cloudInfo& cldInfo
)
{
    return getCloud(cldInfo.regionName(), cldInfo.name());
}


// ************************************************************************* //
