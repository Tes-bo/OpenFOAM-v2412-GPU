/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::schemesLookup

Description
    Selector class for finite area/finite volume differencing schemes.

    The file will contain these types of entries:
    \table
      Name         | Description                        | Type | Reqd | Default
      ddtSchemes   | ddt                                | dict | no   | none
      d2dt2Schemes | d2dt2 (currently only finite-area) | dict | no   | none
      interpolationSchemes | interpolation              | dict | no   | linear
      divSchemes   | div                                | dict | yes  | -
      gradSchemes  | grad                               | dict | yes  | -
      lnGradSchemes | line-normal grad (finite-area)    | dict | no | corrected
      snGradSchemes | surface-normal grad               | dict | no | corrected
      laplacianSchemes | laplacian                      | dict | yes  | -
      fluxRequired |                                    | dict | no | -
    \endtable

    \par select
    This is a rarely used feature to select between different groups of
    settings within the file. If unspecified, the file contents are used
    directly.

Note
    The schemes data are treated as \c READ_MODIFIED even if the
    requested \p readOption is nominally MUST_READ or READ_IF_PRESENT.
    This allows run-time modifications to behave as expected. <br>
    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SourceFiles
    schemesLookup.C
    schemesLookupDetail.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_schemesLookup_H
#define Foam_schemesLookup_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class schemesLookup Declaration
\*---------------------------------------------------------------------------*/

class schemesLookup
:
    public IOdictionary
{
    // Private Class

        //- Lightweight grouping of scheme dictionary and default
        struct lookupDetail
        {
            word name_;
            dictionary dict_;
            ITstream default_;


        // Constructors

            //- Construct empty with given sub-dictionary name
            lookupDetail
            (
                const word& dictName,
                const fileName& parentDictPath
            );


        // Member Functions

            //- Clear dictionary and default scheme
            void clear();

            //- Return the default scheme (if any)
            ITstream& fallback() const;

            //- Lookup named scheme from dictionary, or return default
            ITstream& lookup(const word& name) const;

            //- Populate dictionary and/or default
            void populate
            (
                const dictionary& dict,
                const word& defaultName,
                const bool mandatory = false
            );

            //- Write dictionary entry
            void writeEntry(Ostream& os) const;

            //- Write dictionary entry if non-empty
            void writeEntryOptional(Ostream& os) const;
        };


    // Private Data

        //- ddt
        lookupDetail ddtSchemes_;

        //- d2dt2
        lookupDetail d2dt2Schemes_;

        //- interpolation
        lookupDetail interpSchemes_;

        //- div
        lookupDetail divSchemes_;

        //- grad
        lookupDetail gradSchemes_;

        //- lnGrad (finiteArea)
        lookupDetail lnGradSchemes_;

        //- snGrad (finiteVolume)
        lookupDetail snGradSchemes_;

        //- laplacian
        lookupDetail laplacianSchemes_;

        //- flux
        mutable dictionary fluxRequired_;
        bool fluxRequiredDefault_;

        //- True if default ddtScheme is steady-state
        bool steady_;


    // Private Member Functions

        //- Clear dictionaries and streams before reading
        void clear();

        //- Check if default ddtScheme is steady-state
        void checkSteady();

        //- Read settings from the dictionary
        void read(const dictionary& dict);

        //- The entire dictionary or the optional "select" sub-dictionary.
        const dictionary& selectedDict() const;


    // Generated Methods

        //- No copy construct
        schemesLookup(const schemesLookup&) = delete;

        //- No copy assignment
        void operator=(const schemesLookup&) = delete;


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        schemesLookup
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const word& dictName,
            const dictionary* fallback = nullptr
        );

        //- Construct for objectRegistry, (system) dictionary name
        //- using the readOption from the registry.
        schemesLookup
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        );


    // Member Functions

        //- True if default ddt scheme is steady-state
        bool steady() const noexcept { return steady_; }

        //- True if default ddt scheme is not steady-state
        bool transient() const noexcept { return !steady_; }

        //- The entire dictionary or the optional "select" sub-dictionary.
        const dictionary& schemesDict() const;


    // Lookup Access

        //- Get ddt scheme for given name, or default
        ITstream& ddtScheme(const word& name) const;

        //- Get d2dt2 scheme for given name, or default
        ITstream& d2dt2Scheme(const word& name) const;

        //- Get interpolation scheme for given name, or default
        ITstream& interpolationScheme(const word& name) const;

        //- Get div scheme for given name, or default
        ITstream& divScheme(const word& name) const;

        //- Get grad scheme for given name, or default
        ITstream& gradScheme(const word& name) const;

        //- Get (finiteArea) lnGrad scheme for given name, or default
        ITstream& lnGradScheme(const word& name) const;

        //- Get (finiteVolume) snGrad scheme for given name, or default
        ITstream& snGradScheme(const word& name) const;

        //- Get laplacian scheme for given name, or default
        ITstream& laplacianScheme(const word& name) const;

        //- Set flux-required for given name (mutable)
        void setFluxRequired(const word& name) const;

        //- Get flux-required for given name, or default
        bool fluxRequired(const word& name) const;


    // Read Access

        //- Access ddt schemes dictionary
        const dictionary& ddtSchemes() const noexcept
        {
            return ddtSchemes_.dict_;
        }

        //- Access d2dt2 schemes dictionary
        const dictionary& d2dt2Schemes() const noexcept
        {
            return d2dt2Schemes_.dict_;
        }

        //- Access interpolation schemes dictionary
        const dictionary& interpolationSchemes() const noexcept
        {
            return interpSchemes_.dict_;
        }

        //- Access div schemes dictionary
        const dictionary& divSchemes() const noexcept
        {
            return divSchemes_.dict_;
        }

        //- Access grad schemes dictionary
        const dictionary& gradSchemes() const noexcept
        {
            return gradSchemes_.dict_;
        }

        //- Access lnGrad schemes dictionary (finiteArea)
        const dictionary& lnGradSchemes() const noexcept
        {
            return lnGradSchemes_.dict_;
        }

        //- Access snGrad schemes dictionary (finiteVolume)
        const dictionary& snGradSchemes() const noexcept
        {
            return snGradSchemes_.dict_;
        }

        //- Access laplacian schemes dictionary
        const dictionary& laplacianSchemes() const noexcept
        {
            return laplacianSchemes_.dict_;
        }

        //- Access flux-required dictionary
        const dictionary& fluxRequired() const noexcept
        {
            return fluxRequired_;
        }


    // Edit Access

        //- Access ddt schemes dictionary
        dictionary& ddtSchemes() noexcept
        {
            return ddtSchemes_.dict_;
        }

        //- Access d2dt2 schemes dictionary
        dictionary& d2dt2Schemes() noexcept
        {
            return d2dt2Schemes_.dict_;
        }

        //- Access interpolation schemes dictionary
        dictionary& interpolationSchemes() noexcept
        {
            return interpSchemes_.dict_;
        }

        //- Access div schemes dictionary
        dictionary& divSchemes() noexcept
        {
            return divSchemes_.dict_;
        }

        //- Access grad schemes dictionary
        dictionary& gradSchemes() noexcept
        {
            return gradSchemes_.dict_;
        }

        //- Access lnGrad schemes dictionary (finiteArea)
        dictionary& lnGradSchemes() noexcept
        {
            return lnGradSchemes_.dict_;
        }

        //- Access snGrad schemes dictionary (finiteVolume)
        dictionary& snGradSchemes() noexcept
        {
            return snGradSchemes_.dict_;
        }

        //- Access laplacian schemes dictionary
        dictionary& laplacianSchemes() noexcept
        {
            return laplacianSchemes_.dict_;
        }

        //- Access flux-required dictionary
        dictionary& fluxRequired() noexcept
        {
            return fluxRequired_;
        }


    // Read

        //- Read schemes from IOdictionary, respects the "select" keyword
        bool read();


    // Write

        //- Write dictionary (possibly modified) settings
        void writeDicts(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for schemesLookup - same content for all ranks
template<>
struct is_globalIOobject<schemesLookup> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
