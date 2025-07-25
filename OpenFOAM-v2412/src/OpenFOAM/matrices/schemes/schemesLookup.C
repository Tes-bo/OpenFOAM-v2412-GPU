/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "schemesLookup.H"
#include "Switch.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::schemesLookup::debug(Foam::debug::debugSwitch("schemesLookup", 0));


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::schemesLookup::clear()
{
    ddtSchemes_.clear();
    d2dt2Schemes_.clear();
    interpSchemes_.clear();
    divSchemes_.clear();        // optional
    gradSchemes_.clear();       // optional
    lnGradSchemes_.clear();
    snGradSchemes_.clear();
    laplacianSchemes_.clear();  // optional

    // Do not clear fluxRequired settings
}


void Foam::schemesLookup::checkSteady()
{
    ITstream& is = ddtSchemes_.fallback();

    word schemeName;
    if (is.peek().isWord())
    {
        is >> schemeName;
    }

    // OR: steady_ = schemeName.starts_with("steady");
    steady_ =
    (
        schemeName == "steady"
     || schemeName == "steadyState"
    );
}


void Foam::schemesLookup::read(const dictionary& dict)
{
    ddtSchemes_.populate(dict, "none");
    d2dt2Schemes_.populate(dict, "none");
    interpSchemes_.populate(dict, "linear");
    divSchemes_.populate(dict, "", true);           // Mandatory entry
    gradSchemes_.populate(dict, "", true);          // Mandatory entry
    lnGradSchemes_.populate(dict, "corrected");     // (finiteArea)
    snGradSchemes_.populate(dict, "corrected");     // (finiteVolume)
    laplacianSchemes_.populate(dict, "", true);     // Mandatory entry

    const dictionary* fluxDictPtr = dict.findDict("fluxRequired");
    if (fluxDictPtr)
    {
        fluxRequired_.merge(*fluxDictPtr);

        if (fluxRequired_.found("default"))
        {
            Switch sw(fluxRequired_.lookup("default").peek());

            if (sw.good() && sw.type() != Switch::NONE)
            {
                fluxRequiredDefault_ = bool(sw);
            }
        }
    }

    checkSteady();
}


const Foam::dictionary& Foam::schemesLookup::selectedDict() const
{
    word select;

    if (readIfPresent("select", select, keyType::LITERAL))
    {
        return subDict(select);
    }

    return *this;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::schemesLookup::schemesLookup
(
    const objectRegistry& obr,
    IOobjectOption::readOption rOpt,
    const word& dictName,
    const dictionary* fallback
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            rOpt,
            IOobjectOption::NO_WRITE
        ),
        fallback
    ),

    // Named, but empty dictionaries and default schemes

    ddtSchemes_("ddtSchemes", objectPath()),
    d2dt2Schemes_("d2dt2Schemes", objectPath()),
    interpSchemes_("interpolationSchemes", objectPath()),
    divSchemes_("divSchemes", objectPath()),
    gradSchemes_("gradSchemes", objectPath()),
    lnGradSchemes_("lnGradSchemes", objectPath()),
    snGradSchemes_("snGradSchemes", objectPath()),
    laplacianSchemes_("laplacianSchemes", objectPath()),

    fluxRequired_(objectPath() / "fluxRequired"),
    fluxRequiredDefault_(false),
    steady_(false)
{
    // Treat as READ_MODIFIED whenever possible
    if
    (
        readOpt() == IOobjectOption::MUST_READ
     || (isReadOptional() && headerOk())
    )
    {
        readOpt(IOobjectOption::READ_MODIFIED);
        addWatch();
    }

    // Update: from values read or copied in
    if
    (
        readOpt() == IOobjectOption::READ_MODIFIED
     || !dictionary::empty()
    )
    {
        read(selectedDict());
    }
}


Foam::schemesLookup::schemesLookup
(
    const objectRegistry& obr,
    const word& dictName,
    const dictionary* fallback
)
:
    schemesLookup(obr, obr.readOpt(), dictName, fallback)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::schemesLookup::read()
{
    if (regIOobject::read())
    {
        clear();  // Clear current settings except fluxRequired

        read(selectedDict());

        return true;
    }

    return false;
}


Foam::ITstream& Foam::schemesLookup::ddtScheme(const word& name) const
{
    DebugInfo<< "Lookup ddt scheme for " << name << endl;
    return ddtSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::d2dt2Scheme(const word& name) const
{
    DebugInfo<< "Lookup d2dt2 scheme for " << name << endl;
    return d2dt2Schemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::interpolationScheme(const word& name) const
{
    DebugInfo<< "Lookup interpolation scheme for " << name << endl;
    return interpSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::divScheme(const word& name) const
{
    DebugInfo<< "Lookup div scheme for " << name << endl;
    return divSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::gradScheme(const word& name) const
{
    DebugInfo<< "Lookup grad scheme for " << name << endl;
    return gradSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::lnGradScheme(const word& name) const
{
    DebugInfo<< "Lookup lnGrad scheme for " << name << endl;
    return lnGradSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::snGradScheme(const word& name) const
{
    DebugInfo<< "Lookup snGrad scheme for " << name << endl;
    return snGradSchemes_.lookup(name);
}


Foam::ITstream& Foam::schemesLookup::laplacianScheme(const word& name) const
{
    DebugInfo<< "Lookup laplacian scheme for " << name << endl;
    return laplacianSchemes_.lookup(name);
}


void Foam::schemesLookup::setFluxRequired(const word& name) const
{
    DebugInfo<< "Setting fluxRequired for " << name << endl;
    fluxRequired_.add(name, true, true);
}


bool Foam::schemesLookup::fluxRequired(const word& name) const
{
    DebugInfo<< "Lookup fluxRequired for " << name << endl;
    return (fluxRequired_.found(name) || fluxRequiredDefault_);
}


const Foam::dictionary& Foam::schemesLookup::schemesDict() const
{
    return selectedDict();
}


void Foam::schemesLookup::writeDicts(Ostream& os) const
{
    ddtSchemes_.writeEntryOptional(os);
    d2dt2Schemes_.writeEntryOptional(os);
    interpSchemes_.writeEntryOptional(os);
    divSchemes_.writeEntry(os);  // Mandatory entry
    gradSchemes_.writeEntry(os); // Mandatory entry
    lnGradSchemes_.writeEntryOptional(os);  // (finiteArea)
    snGradSchemes_.writeEntryOptional(os);  // (finiteVolume)
    laplacianSchemes_.writeEntry(os); // Mandatory entry

    if (!fluxRequired_.empty())
    {
        fluxRequired_.writeEntry(os);
    }
}


// ************************************************************************* //
