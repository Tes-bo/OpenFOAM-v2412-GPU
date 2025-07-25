/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::argList

Description
    Extract command arguments and options from the supplied
    \a argc and \a argv parameters.

    Sequences with "(" ... ")" are transformed into a stringList.
    For example,
    \verbatim
        program -listFiles \( *.txt \)
    \endverbatim
    would create a stringList:
    \verbatim
        ( "file1.txt" "file2.txt" ... "fileN.txt" )
    \endverbatim
    The backslash-escaping is required to avoid interpretation by the shell.

    Default command-line options:
      - \par -case \<dir\>
        Select a case directory instead of the current working directory
      - \par -decomposeParDict \<file\>
        Read decomposePar dictionary from specified location
      - \par -parallel
        Specify case as a parallel job
      - \par -doc
        Display the documentation in browser
      - \par -srcDoc
        Display the source documentation in browser
      - \par -help
        Print the usage

    Additionally, the \b -noFunctionObjects and \b -postProcess options
    may be present for some solvers or utilities.

    Environment variables set by argList or by Time:
      - \par FOAM_API
        The value of foamVersion::api
      - \par FOAM_CASE
        The path of the global case.
        It is the same for serial and parallel jobs.
      - \par FOAM_CASENAME
        The name of the global case.
      - \par FOAM_EXECUTABLE
        If not already present in the calling environment,
        it is set to the \a name portion of the calling executable.
      - \par FOAM_APPLICATION
        If not already present in the calling environment,
        it is set to the value of the \c application entry
        (from \c controlDict) if that entry is present.

    The value of the \b FOAM_APPLICATION may be inconsistent if the value of
    the \c application entry is adjusted during runtime.

Note
    - The document browser used is defined by the \b FOAM_DOC_BROWSER
      environment variable or the <tt>Documentation/docBrowser</tt> entry
      in the <tt>\<etc\>/controlDict</tt> file.
      The \%f token is used as a placeholder for the file name.
    - The valid (mandatory) arguments can be adjusted
      via the addArgument static method instead of directly
      manipulating the argList::validArgs static member.
    - The valid options can be adjusted
      via the addOption/removeOption static methods instead of directly
      manipulating the argList::validOptions static member.

SourceFiles
    argList.C
    argListI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_argList_H
#define Foam_argList_H

#include "stringList.H"
#include "SLList.H"
#include "HashSet.H"
#include "fileName.H"
#include "parRun.H"  // "ParRunControl"
#include "ITstream.H"
#include "dlLibraryTable.H"
#include "OSspecific.H"

// Transitional features - older style access (including 1712 release)
#define Foam_argList_1712

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class argList Declaration
\*---------------------------------------------------------------------------*/

class argList
{
    // Private Data

        //- Command arguments are mandatory (default) or optional
        static bool argsMandatory_;

        //- Check presence of processor directories (default: on)
        static bool checkProcessorDirectories_;

        //- MPI threads are desired for the application (default: off)
        static bool parallelThreads_;

        //- Switch on/off parallel mode, dry-run etc.
        //  Construct first so destructor is done last.
        ParRunControl runControl_;

        //- The arguments after removing known options
        stringList args_;

        //- The extracted options
        HashTable<string> options_;

        //- Additional libraries
        mutable dlLibraryTable libs_;

        word executable_;
        fileName rootPath_;
        fileName globalCase_;
        fileName case_;

        //- The command line options and arguments concatenated as a string
        string commandLine_;


    // Private Member Functions

        //- Helper for resolving aliases for -options within validOptionsCompat
        static word optionCompat(const word& optName);

        //- Helper for resolving ignored options
        static int optionIgnore(const word& optName);

        //- Check after reading if the input token stream has unconsumed
        //- tokens remaining or if there were no tokens in the first place.
        //  Emits Warning
        static void checkITstream(const ITstream& is, const label index);

        //- Check after reading if the input token stream has unconsumed
        //- tokens remaining or if there were no tokens in the first place.
        //  Emits Warning
        static void checkITstream(const ITstream& is, const word& optName);

        //- Read a List of values from ITstream,
        //- treating a single entry like a list of size 1.
        template<class T>
        static inline void readList(ITstream& is, List<T>& list);

        //- Trigger FatalError for given option
        void raiseBadInput(const word& optName) const;

        //- Set rootPath_, globalCase_, case_ from one of the following forms
        //   * [-case dir]
        //   * cwd
        //
        // Exports FOAM_CASE and FOAM_CASENAME env variables so they can
        // be used immediately (eg, in decomposeParDict).
        // Exports FOAM_EXECUTABLE env variable.
        //
        // Detects -dry-run option
        void setCasePaths();

        //- Transcribe argv into internal args_.
        //  Transform sequences with "(" ... ")" into string lists
        //  return true if any "(" ... ")" sequences were captured
        bool regroupArgv(int& argc, char**& argv);

        //- Print option compatibility
        void printManCompat() const;


public:

    // Static Data Members

        //- A list of valid (mandatory) arguments
        static SLList<string> validArgs;

        //- The "advanced" options are shown with -help-full (not with --help)
        static HashSet<string> advancedOptions;

        //- A list of valid options
        static HashTable<string> validOptions;

        //- A list of valid parallel options
        static HashTable<string> validParOptions;

        //- A list of aliases for options.
        //  Stored as (alias = canonical, version)
        static HashTable<std::pair<word,int>> validOptionsCompat;

        //- A list of options to ignore.
        //  Stored as (option = bool, version)
        static HashTable<std::pair<bool,int>> ignoreOptionsCompat;

        //- Short description for program arguments
        static HashTable<string, label, Hash<label>> argUsage;

        //- Short description for validOptions
        static HashTable<string> optionUsage;

        //- General usage notes
        static SLList<string> notes;

        //- Min indentation when displaying usage (default: 20)
        static std::string::size_type usageMin;

        //- Max screen width when displaying usage (default: 80)
        static std::string::size_type usageMax;

        //- Standard name for the post-processing option
        static word postProcessOptionName;

        //! \cond internalClass
        //  The constructor populates the standard options
        struct initValidTables
        {
            initValidTables();
        };
        //! \endcond


    // Constructors

        //- Construct from argc and argv
        //- checking the arguments and options as requested.
        //
        //  By default, the argument check respects the value of
        //  argsMandatory() to decide if the arguments should be checked
        //  (when they are mandatory) or not (when they are optional)
        argList
        (
            int& argc,
            char**& argv,
            bool checkArgs = argList::argsMandatory(),
            bool checkOpts = true,
            bool initialise = true
        );

        //- Construct copy with new options
        argList
        (
            const argList& args,
            const HashTable<string>& options,
            bool checkArgs = true,
            bool checkOpts = true,
            bool initialise = true
        );


    //- Destructor
    virtual ~argList();


    // Member Functions

    // Environment

        //- Name of the executable from environment variable
        //
        //  Returns the contents of the \c FOAM_EXECUTABLE variable,
        //  which has previously been set by argList.
        //
        //  This will normally be identical to the value of executable(),
        //  but obtained from the environment.
        static word envExecutable();

        //- Global case (directory) from environment variable
        //
        //  Returns the contents of the \c FOAM_CASE variable,
        //  which has previously been set by argList or by Time.
        //
        //  This will normally be identical to the value of globalPath(),
        //  but obtained via the environment.
        static fileName envGlobalPath();

        //- Return the input relative to the globalPath by stripping off
        //- a leading value of the envGlobalPath
        //
        //  \param input the directory or filename to make case-relative
        //  \param caseTag replace globalPath with \<case\> for later
        //      use with expand(), or prefix \<case\> if the file name was
        //      not an absolute location
        static fileName envRelativePath
        (
            const fileName& input,
            const bool caseTag = false
        );


    // Low-level

        //- Scan for -help, -doc options etc prior to checking the validity
        //- of other args/opts and finally initialising.
        void parse(bool checkArgs, bool checkOpts, bool initialise);


    // Access

        //- Name of executable without the path
        inline const word& executable() const noexcept;

        //- The command line options and arguments concatenated as a string
        inline const string& commandLine() const noexcept;

        //- Return root path
        inline const fileName& rootPath() const noexcept;

        //- Return case name (parallel run) or global case (serial run)
        inline const fileName& caseName() const noexcept;

        //- Return global case name
        inline const fileName& globalCaseName() const noexcept;

        //- Return the full path to the (processor local) case
        //  \note This is guaranteed to be an absolute path
        inline fileName path() const;

        //- Return the full path to the global case
        //  \note This is guaranteed to be an absolute path
        inline fileName globalPath() const;

        //- Return the input relative to the globalPath by stripping off
        //- a leading value of the globalPath
        //
        //  \param input the directory or filename to make case-relative
        //  \param caseTag replace globalPath with \<case\> for later
        //      use with expand(), or prefix \<case\> if the file name was
        //      not an absolute location
        inline fileName relativePath
        (
            const fileName& input,
            const bool caseTag = false
        ) const;

        //- Return the run control (parallel, dry-run etc)
        inline const ParRunControl& runControl() const noexcept;

        //- Return distributed flag
        //- (i.e. are rootPaths different on different machines)
        inline bool distributed() const noexcept;

        //- Return the dry-run flag
        inline int dryRun() const noexcept;

        //- Modify the dry-run flag
        inline int dryRun(const int level) noexcept;

        //- Return the verbose flag
        inline int verbose() const noexcept;

        //- Modify the verbose flag
        inline int verbose(const int level) noexcept;

        //- The controlDict 'functions' entry is allowed to be used.
        //  There must be a command-line option to enable/disable
        //  (-withFunctionObjects, -noFunctionObjects)
        //  and has not been explicitly disabled
        bool allowFunctionObjects() const;

        //- The controlDict 'libs' entry is allowed to be used.
        //- (eg, has not been disabled by the -no-libs option)
        bool allowLibs() const;

        //- Mutable access to the loaded dynamic libraries
        inline dlLibraryTable& libs() const noexcept;

        //- The number of arguments
        inline label size() const noexcept;

        //- Return arguments
        inline const stringList& args() const noexcept;

        //- Non-const access to the command arguments (non-options)
        inline stringList& args() noexcept;

        //- Return options
        inline const HashTable<string>& options() const noexcept;

        //- Return non-const access to the command options
        inline HashTable<string>& options() noexcept;

        //- Return true if the named option is found
        inline bool found(const word& optName) const;

        //- Return how many of the specified options were used
        label count(const UList<word>& optionNames) const;

        //- Return how many of the specified options were used
        label count(std::initializer_list<word> optionNames) const;

        //- Return an input stream from the named option
        inline ITstream lookup(const word& optName) const;

        //- Get a value from the argument at index.
        //  Index 1 is the first (non-option) argument.
        //  For fileName type, invokes fileName::validate()
        template<class T>
        inline T get(const label index) const;

        //- Get a List of values from the argument at index.
        //  Index 1 is the first (non-option) argument.
        template<class T>
        inline List<T> getList(const label index) const;

        //- Get a value from the named option
        //  The default template parameter is string (ie, no conversion).
        //  For fileName type, invokes fileName::validate()
        template<class T=string>
        inline T get(const word& optName) const;

        //- Get a value from the named option if present, or return default.
        template<class T>
        inline T getOrDefault(const word& optName, const T& deflt) const;

        //- Get a List of values from the named option,
        //- treating a single entry like a list of size 1.
        //  \param optName  the option name to read from
        //  \param mandatory if the option is non-mandatory, the behaviour
        //      is similar to readListIfPresent().
        template<class T>
        inline List<T> getList(const word& optName, bool mandatory=true) const;

        //- Read a value from the named option if present.
        //  \return true if the named option was found.
        template<class T>
        inline bool readIfPresent(const word& optName, T& val) const;

        //- Read a value from the named option if present.
        //  \return true if the named option was found, otherwise
        //  use the supplied default and return false.
        template<class T>
        inline bool readIfPresent
        (
            const word& optName,
            T& val,
            const T& deflt
        ) const;

        //- If named option is present, get a List of values
        //- treating a single entry like a list of size 1.
        //  \return true if the named option was found.
        template<class T>
        inline bool readListIfPresent(const word& optName, List<T>& list) const;

        //- Read the named option and check its validity.
        //  FatalError if mandatory and not found, or if the predicate check
        //  failed.
        //
        //  \param optName the option name
        //  \param val the value to read into
        //  \param pred the value check predicate
        //
        //  \return true if the entry was found.
        template<class T, class Predicate>
        inline bool readCheck
        (
            const word& optName,
            T& val,
            const Predicate& pred,
            bool mandatory = true
        ) const;

        //- Read the named option if present and check its validity.
        //  FatalError if found and the predicate check failed.
        //
        //  \param optName the option name
        //  \param val the value to read into
        //  \param pred the value check predicate
        //
        //  \return true if the entry was found.
        template<class T, class Predicate>
        inline bool readCheckIfPresent
        (
            const word& optName,
            T& val,
            const Predicate& pred
        ) const;

        //- Get a value from the named option with additional checking.
        //  FatalError if the predicate check failed.
        //
        //  \param optName the option name
        //  \param pred the value check predicate
        template<class T, class Predicate>
        T getCheck
        (
            const word& optName,
            const Predicate& pred
        ) const;

        //- Get a value from the named option with additional checking
        //- (if present), or return default.
        //  FatalError if the predicate check on the retrieved value failed.
        //
        //  \param optName the option name
        //  \param deflt the default return value
        //  \param pred the value check predicate
        template<class T, class Predicate>
        T getCheckOrDefault
        (
            const word& optName,
            const T& deflt,
            const Predicate& pred
        ) const;


    // Edit

        //- Append a (mandatory) argument to validArgs
        static void addArgument
        (
            const string& argName,
            const string& usage = ""
        );

        //- Add a bool option to validOptions with usage information
        static void addBoolOption
        (
            const word& optName,
            const string& usage = "",
            bool advanced = false
        );

        //- Add an option to validOptions with usage information
        //  An option with an empty param is a bool option
        static void addOption
        (
            const word& optName,
            const string& param = "",
            const string& usage = "",
            bool advanced = false
        );

        //- Set an existing option as being 'advanced' or normal
        static void setAdvanced(const word& optName, bool advanced = true);

        //- Specify an alias for the option name.
        //
        //  \param optName the currently used option name
        //  \param compat alias name and the last OpenFOAM version (YYMM)
        //      when the alias was not needed.
        //      Setting a zero or negative version suppresses warnings about
        //      the alias.
        static void addOptionCompat
        (
            const word& optName,
            std::pair<const char*, int> compat
        );

        //- Specify an option to be ignored.
        //
        //  \param compat optName and the last OpenFOAM version (YYMM)
        //      when the option was directly supported.
        //      Setting a zero or negative version suppresses warnings about
        //      the alias.
        //  \param expectArg the option is non-bool
        static void ignoreOptionCompat
        (
            std::pair<const char*,int> compat,
            bool expectArg
        );

        //- Add option usage information to optionUsage
        static void addUsage
        (
            const word& optName,
            const string& usage
        );

        //- Add extra notes for the usage information
        //  This string is used "as-is" without additional formatting
        static void addNote(const string& note);

        //- Remove option from validOptions and from optionUsage
        static void removeOption(const word& optName);

        //- Flag command arguments as being optional (non-mandatory)
        static void noMandatoryArgs();

        //- Command arguments type (optional/mandatory).
        static bool argsMandatory();

        //- Disable emitting the banner information.
        //  Adjusts the Foam::infoDetailLevel flag.
        static void noBanner();

        //- Banner status (enabled/disabled).
        //  Queries the Foam::infoDetailLevel flag.
        static bool bannerEnabled();

        //- Enable a 'dry-run' bool option, with usage information
        static void addDryRunOption
        (
            const string& usage,  //! usage information (expected)
            bool advanced = false
        );

        //- Enable a 'verbose' bool option, with usage information
        static void addVerboseOption
        (
            const string& usage = "", //! usage information (optional)
            bool advanced = false
        );

        //- Remove '-noFunctionObjects' option and ignore any occurrences.
        //  Optionally add a '-withFunctionObjects' option instead
        static void noFunctionObjects(bool addWithOption = false);

        //- Suppress JobInfo, overriding controlDict setting
        static void noJobInfo();

        //- Add the '-no-libs' command line option
        static void noLibs();

        //- Remove the parallel options
        static void noParallel();

        //- Disable checking of processor directories
        static void noCheckProcessorDirectories();

        //- MPI threads are desired for the application
        static void parallelThreads_on();

        //- Set option directly (use with caution)
        //  An option with an empty param is a bool option.
        //  Not all valid options can also be set: eg, -case, -roots, ...
        //  Return true if the existing option value needed changing,
        //  or if the option did not previously exist.
        bool setOption(const word& optName, const string& param = "");

        //- Unset option directly (use with caution)
        //  Not all valid options can also be unset: eg, -case, -roots ...
        //  Return true if the option existed before being unset.
        bool unsetOption(const word& optName);


    // Helpers

        //- True if the post-processing option is found in the \c argv list
        static bool postProcess(int argc, char *argv[]);

        //- The number of times -verbose is found in the \c argv list
        static int verbose(int argc, char *argv[]);


    // Print

        //- Print option compatibility
        void printCompat() const;

        //- Print notes (if any)
        void printNotes() const;

        //- Print usage
        void printUsage(bool full=true) const;

        //- Print usage as nroff-man format (Experimental)
        void printMan() const;

        //- Display documentation in browser
        //  Optionally display the application source code
        void displayDoc(bool source=false) const;


    // Check

        //- Check the parsed command-line for mandatory arguments and
        //- that all the options are correct.
        //
        //  By default, the argument check respects the value of
        //  argsMandatory() to decide if the arguments should be checked
        //  (when they are mandatory) or not (when they are optional)
        bool check
        (
            bool checkArgs = argList::argsMandatory(),
            bool checkOpts = true
        ) const;

        //- Check root path and case path
        bool checkRootCase() const;


    // Member Operators

        //- The string corresponding to the argument index.
        //  Index 0 is the executable.
        //  Index 1 is the first (non-option) argument.
        inline const string& operator[](const label index) const;

        //- The string associated with the named option
        inline const string& operator[](const word& optName) const;


    // Housekeeping

        //- Deprecated(2020-05) identical to get(const word& optName)
        //  \deprecated(2020-05) - use get() method
        template<class T=string>
        FOAM_DEPRECATED_STRICT(2020-06, "get()")
        T opt(const word& optName) const
        {
            return this->get<T>(optName);
        }

        //- Deprecated(2020-05) identical to getOrDefault(...)
        //  \deprecated(2020-05) - use getOrDefault() method
        template<class T>
        FOAM_DEPRECATED_STRICT(2020-06, "getOrDefault()")
        T opt(const word& optName, const T& deflt) const
        {
            return this->getOrDefault<T>(optName, deflt);
        }

        //- Deprecated(2020-05) identical to getOrDefault(...)
        //  \deprecated(2020-05) - use getOrDefault() method
        template<class T>
        FOAM_DEPRECATED_STRICT(2020-06, "getOrDefault()")
        T get(const word& optName, const T& deflt) const
        {
            return this->getOrDefault<T>(optName, deflt);
        }

        //- Deprecated(2020-05) identical to getOrDefault(...)
        //  \deprecated(2020-05) - use getOrDefault() method
        template<class T>
        FOAM_DEPRECATED_STRICT(2020-06, "getOrDefault()")
        T lookupOrDefault(const word& optName, const T& deflt) const
        {
            return this->getOrDefault<T>(optName, deflt);
        }

        //- Same as runControl() - v2106 and earlier
        const ParRunControl& parRunControl() const { return runControl_; }


    // Older style access (including 1712 release)

    #ifdef Foam_argList_1712

        //- Deprecated(2018-08) read a value from the argument at index.
        //  Index 1 is the first (non-option) argument.
        //  \deprecated(2018-08) - use get() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-08, "get() method")
        T read(const label index) const
        {
            return this->get<T>(index);
        }

        //- Deprecated(2018-01) read a value from the argument at index.
        //  Index 1 is the first (non-option) argument.
        //  \deprecated(2018-01) - use get() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "get() method")
        T argRead(const label index) const
        {
            return this->get<T>(index);
        }

        //- Deprecated(2018-01) return true if the named option is found
        //  \deprecated(2018-01) - use found() method
        FOAM_DEPRECATED_FOR(2018-01, "found() method")
        bool optionFound(const word& optName) const
        {
            return found(optName);
        }

        //- Deprecated(2018-01) return an input stream from the named option
        //  \deprecated(2018-01) - use lookup() method
        FOAM_DEPRECATED_FOR(2018-01, "lookup() method")
        ITstream optionLookup(const word& optName) const
        {
            return lookup(optName);
        }

        //- Deprecated(2018-01) read a value from the named option
        //  \deprecated(2018-01) - use get() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "get() method")
        T optionRead(const word& optName) const
        {
            return get<T>(optName);
        }

        //- Deprecated(2018-01) read a value from the named option if present.
        //  Return true if the named option was found.
        //  \deprecated(2018-01) - use readIfPresent() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "readIfPresent() method")
        bool optionReadIfPresent
        (
            const word& optName,
            T& val
        ) const
        {
            return readIfPresent<T>(optName, val);
        }

        //- Deprecated(2018-01) read a value from the named option if present.
        //  Return true if the named option was found, otherwise
        //  use the supplied default and return false.
        //  \deprecated(2018-01) - use readIfPresent() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "readIfPresent() method")
        bool optionReadIfPresent
        (
            const word& optName,
            T& val,
            const T& deflt
        ) const
        {
            return readIfPresent<T>(optName, val, deflt);
        }

        //- Deprecated(2018-01) read a value from the named option if present.
        //  Return supplied default otherwise.
        //  \deprecated(2018-01) - use getOrDefault() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "getOrDefault() method")
        T optionLookupOrDefault
        (
            const word& optName,
            const T& deflt
        ) const
        {
            return getOrDefault<T>(optName, deflt);
        }

        //- Deprecated(2018-01) read a List of values from the named option
        //  \deprecated(2018-01) - use getList() method
        template<class T>
        FOAM_DEPRECATED_FOR(2018-01, "getList() method")
        List<T> optionReadList(const word& optName) const
        {
            return this->getList<T>(optName);
        }

    #endif  /* Foam_argList_1712 */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "argListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
