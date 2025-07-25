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

\*---------------------------------------------------------------------------*/

#include "cloudInfo.H"
#include "cloud.H"
#include "kinematicCloud.H"
#include "dictionary.H"
#include "mathematicalConstants.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cloudInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeTabbed(os, "nParcels");
    writeTabbed(os, "mass");
    writeTabbed(os, "Dmax");
    writeTabbed(os, "D10");
    writeTabbed(os, "D32");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::cloudInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    functionObjects::logFiles(obr_, name, dict),
    verbose_(false),
    onExecute_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cloudInfo::read(const dictionary& dict)
{
    parcelSelect_.clear();
    verbose_ = false;
    onExecute_ = false;

    if (regionFunctionObject::read(dict) && logFiles::read(dict))
    {
        logFiles::resetNames(dict.get<wordList>("clouds"));

        Info<< type() << " " << name() << ": ";
        if (names().size())
        {
            Info<< "applying to clouds:" << nl;
            for (const word& cldName : names())
            {
                Info<< "    " << cldName << nl;
            }
            Info<< endl;

            // Actions to define selection
            parcelSelect_ = dict.subOrEmptyDict("selection");

            verbose_ = dict.getOrDefault("verbose", false);
            onExecute_ = dict.getOrDefault("sampleOnExecute", false);
        }
        else
        {
            Info<< "no clouds to be processed" << nl << endl;
        }

        if (writeToFile())
        {
            forAll(names(), cloudi)
            {
                writeFileHeader(files(cloudi));
            }
        }
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::performAction(unsigned request)
{
    if (!request || names().empty())
    {
        return true;
    }

    forAll(names(), cloudi)
    {
        // The reported quantities
        label nTotParcels = 0;
        scalar totMass = 0, Dmax = 0, D10 = 0, D32 = 0;
        bool applyFilter = false;

        const word& cloudName = names()[cloudi];

        const auto* kinCloudPtr = obr_.cfindObject<kinematicCloud>(cloudName);

        if (!kinCloudPtr)
        {
            // Safety
            continue;
        }

        const auto& kinCloud = *kinCloudPtr;
        const auto* plainCloudPtr = isA<cloud>(kinCloud);

        if (!parcelSelect_.empty() && plainCloudPtr)
        {
            const auto& plainCloud = *plainCloudPtr;

            // Filtering - simply use cloud methods

            objectRegistry obrTmp
            (
                IOobject
                (
                    "tmp::cloudInfo::" + cloudName,
                    obr_.time().constant(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            plainCloud.writeObjects(obrTmp);

            // Apply output filter (for the current cloud)
            applyFilter = calculateFilter(obrTmp, log);

            // Expected/required fields
            const auto* diamFldPtr = obrTmp.cfindObject<IOField<scalar>>("d");
            const auto* rhoFldPtr = obrTmp.cfindObject<IOField<scalar>>("rho");
            const auto* nParticleFldPtr =
                obrTmp.cfindObject<IOField<scalar>>("nParticle");

            do
            {
                #undef  doLocalCode
                #define doLocalCode(FldPtr, FldName)                          \
                if (applyFilter && !FldPtr)                                   \
                {                                                             \
                    WarningInFunction                                         \
                        << "Missing \"" << #FldName                           \
                        << "\" field - disabling filter" << nl;               \
                    applyFilter = false;                                      \
                    break;                                                    \
                }

                doLocalCode(diamFldPtr, d);
                doLocalCode(rhoFldPtr, rho);
                doLocalCode(nParticleFldPtr, nParticle);
                #undef doLocalCode
            }
            while (false);

            if (applyFilter)
            {
                // Filtered. Need to do everything by hand!

                const auto& diams = *diamFldPtr;
                const auto& rhos = *rhoFldPtr;
                const auto& nParts = *nParticleFldPtr;

                FixedList<scalar, 4> Dsums(Zero);

                for (const label particlei : parcelAddr_)
                {
                    ++nTotParcels;

                    const scalar d = diams[particlei];
                    const scalar rho = rhos[particlei];
                    const scalar np = nParts[particlei];

                    totMass += np*rho*pow3(d);
                    Dmax = max(Dmax, d);

                    Dsums[0] += np;
                    Dsums[1] += np*(d);
                    Dsums[2] += np*(sqr(d));
                    Dsums[3] += np*(pow3(d));
                }

                reduce(nTotParcels, sumOp<label>());
                reduce(totMass, sumOp<scalar>());
                reduce(Dmax, maxOp<scalar>());
                reduce(Dsums, sumOp<scalar>());

                totMass *= (constant::mathematical::pi/6.0);

                Dmax = max(0, Dmax);
                D10 = Dsums[1]/(max(Dsums[0], VSMALL));
                D32 = Dsums[3]/(max(Dsums[2], VSMALL));
            }
        }

        if (!applyFilter)
        {
            // No filter - use regular cloud methods
            nTotParcels = returnReduce(kinCloud.nParcels(), sumOp<label>());
            totMass = returnReduce(kinCloud.massInSystem(), sumOp<scalar>());

            Dmax = kinCloud.Dmax();
            D10 = kinCloud.Dij(1, 0);
            D32 = kinCloud.Dij(3, 2);
        }

        Log << type() << " " << name() <<  " write:" << nl
            << "    number of parcels : " << nTotParcels << nl
            << "    mass in system    : " << totMass << nl
            << "    maximum diameter  : " << Dmax << nl
            << "    D10 diameter      : " << D10 << nl
            << "    D32 diameter      : " << D32 << nl
            << endl;

        if ((request & ACTION_WRITE) && writeToFile())
        {
            auto& os = files(cloudi);

            writeCurrentTime(os);
            os
                << token::TAB << nTotParcels
                << token::TAB << totMass
                << token::TAB << Dmax
                << token::TAB << D10
                << token::TAB << D32
                << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::write()
{
    return performAction(ACTION_ALL);
}


// ************************************************************************* //
