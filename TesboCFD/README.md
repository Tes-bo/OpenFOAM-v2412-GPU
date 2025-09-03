## About OpenFOAM
OpenFOAM is a free, open source CFD software [released and developed by OpenCFD Ltd since 2004](http://www.openfoam.com/history/).
It has a large user base across most areas of engineering and science, from both commercial and academic organisations.
OpenFOAM has an extensive range of features to solve anything from complex fluid flows involving chemical reactions, turbulence and heat transfer, to acoustics, solid mechanics and electromagnetics.
[See documentation](http://www.openfoam.com/documentation)

OpenFOAM is professionally released every six months to include
customer sponsored developments and contributions from the community -
individual and group contributors, integrations
(eg, from FOAM-extend and OpenFOAM Foundation Ltd) as well as
[governance guided activities](https://www.openfoam.com/governance/).


## License

OpenFOAM is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.  See the file COPYING in this directory or
[http://www.gnu.org/licenses/](http://www.gnu.org/licenses), for a
description of the GNU General Public License terms under which you
may redistribute files.


## OpenFOAM Trademark

OpenCFD Ltd grants use of its OpenFOAM trademark by Third Parties on a
licence basis. ESI Group and OpenFOAM Foundation Ltd are currently
permitted to use the Name and agreed Domain Name. For information on
trademark use, please refer to the
[trademark policy guidelines][link trademark].

Please [contact OpenCFD](http://www.openfoam.com/contact) if you have
any questions about the use of the OpenFOAM trademark.

Violations of the Trademark are monitored, and will be duly prosecuted.


## Using OpenFOAM

If OpenFOAM has already been compiled on your system, simply source
the appropriate `etc/bashrc` or `etc/cshrc` file and get started.
For example, for the OpenFOAM-v2412 version:
```
source /installation/path/OpenFOAM-v2412/etc/bashrc
```

## Compiling OpenFOAM

If you are compiling OpenFOAM from source, please see the relevant
guides:

| Location    | Readme    | Requirements | Build |
|-------------|-----------|--------------|-------|
| [OpenFOAM][repo openfoam] | [readme][link openfoam-readme] | [system requirements][link openfoam-require] | [build][link openfoam-build] |
| [ThirdParty][repo third] | [readme][link third-readme] | [system requirements][link third-require] | [build][link third-build] |


If you need to modify the versions or locations of ThirdParty
software, please read how the
[OpenFOAM configuration][wiki-config] is structured.


## How do I know which version I am currently using?

The value of the `$WM_PROJECT_DIR` or even `$WM_PROJECT_VERSION` are
not guaranteed to have any correspondence to the OpenFOAM release
(API) value. If OpenFOAM has already been compiled, the build-time
information is embedded into each application. For example, as
displayed from `blockMesh -help`:
```
Using: OpenFOAM-com (2012) - visit www.openfoam.com
Build: b830beb5ea-20210429 (patch=210414)
Arch:  LSB;label=32;scalar=64
```
This output contains all of the more interesting information that we need:

| item                  | value         |
|-----------------------|---------------|
| version               | com  (eg, local development branch) |
| api                   | 2012          |
| commit                | b830beb5ea    |
| author date           | 20210429      |
| patch-level           | (20)210414    |
| label/scalar size     | 32/64 bits    |

The Arch information may also include the `solveScalar` size
if different than the `scalar` size.

As can be seen in this example, the git build information is
supplemented by the date when the last change was authored, which can
be helpful when the repository contains local changes. If you simply
wish to know the current API and patch levels directly, the
`wmake -build-info` provides the relevant information even
when OpenFOAM has not yet been compiled:
```
$ wmake -build-info
make
    api = 2012
    patch = 210414
    branch = master
    build = 308af39136-20210426
```
Similar information is available with `foamEtcFile`, using the
`-show-api` or `-show-patch` options. For example,
```
$ foamEtcFile -show-api
2012

$ foamEtcFile -show-patch
210414
```
This output will generally be the easiest to parse for scripts.
The `$FOAM_API` convenience environment variable may not reflect the
patching changes made within the currently active environment and
should be used with caution.


## ThirdParty directory

OpenFOAM normally ships with a directory of 3rd-party software and
build scripts for some 3rd-party software that is either necessary or
at least highly useful for OpenFOAM, but which are not necessarily
readily available on every operating system or cluster installation.

These 3rd-party sources are normally located in a directory parallel
to the OpenFOAM directory. For example,
```
/path/parent
|-- OpenFOAM-v2412
\-- ThirdParty-v2412
```
There are, however, many cases where this simple convention is inadequate:

* When no additional 3rd party software is actually required (ie, the
  operating system or cluster installation provides it)

* When we have changed the OpenFOAM directory name to some arbitrary
  directory name, e.g. openfoam-sandbox2412, etc..

* When we would like any additional 3rd party software to be located
  inside of the OpenFOAM directory to ensure that the installation is
  encapsulated within a single directory structure. This can be
  necessary for cluster installations, or may simply be a convenient
  means of performing a software rollout for individual workstations.

* When we have many different OpenFOAM directories for testing or
  developing various different features but wish to use or reuse the
  same 3rd party software for them all.

The solution for these problems is a newer, more intelligent discovery
when locating the ThirdParty directory with the following precedence:

1. PROJECT/ThirdParty
   * for single-directory installations
2. PREFIX/ThirdParty-VERSION
   * this corresponds to the traditional approach
3. PREFIX/ThirdParty-vAPI
   * allows for an updated value of VERSION, *eg*, `v2412-myCustom`,
     without requiring a renamed ThirdParty. The API value would still
     be `2412` and the original `ThirdParty-v2412/` would be found.
4. PREFIX/ThirdParty-API
   * same as the previous example, but using an unadorned API value.
5. PREFIX/ThirdParty-common
   * permits maximum reuse for various versions, for experienced
     users who are aware of potential version incompatibilities

If none of these directories are found to be suitable, it reverts to
using PROJECT/ThirdParty as a dummy location (even if the directory
does not exist). This is a safe fallback value since it is within the
OpenFOAM directory structure and can be trusted to have no negative
side-effects. In the above, the following notation has been used:

| name          | value         | meaning       |
|---------------|---------------|---------------|
| PROJECT       | `$WM_PROJECT_DIR`     | The OpenFOAM directory |
| PREFIX        | `dirname $WM_PROJECT_DIR` | The OpenFOAM parent directory |
| API           | `foamEtcFiles -show-api` |  The api or release version |
| VERSION       | `$WM_PROJECT_VERSION` | The version we have chosen |

To reduce the potential of false positive matches (perhaps some other
software also uses ThirdParty-xxx for its naming), the directory test
is accompanied by a OpenFOAM-specific sanity test. The OpenFOAM
ThirdParty directory will contain either an `Allwmake` file or a
`platforms/` directory.


<!-- OpenFOAM -->

[link trademark]: https://www.openfoam.com/opencfd-limited-trade-mark-policy

[repo openfoam]: https://develop.openfoam.com/Development/openfoam/
[repo third]: https://develop.openfoam.com/Development/ThirdParty-common/

[link openfoam-readme]: https://develop.openfoam.com/Development/openfoam/blob/develop/README.md
[link openfoam-issues]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/BuildIssues.md
[link openfoam-build]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md
[link openfoam-require]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/Requirements.md

[wiki-config]: https://develop.openfoam.com/Development/openfoam/-/wikis/configuring


## Useful Links

- Download [source](https://dl.openfoam.com/source/) and [download and installation instructions](http://www.openfoam.com/download/)
- [Documentation](http://www.openfoam.com/documentation)
- [Reporting bugs/issues/feature requests](http://www.openfoam.com/code/bug-reporting.php)
- [Issue tracker](https://develop.openfoam.com/Development/openfoam/-/issues)
- [Code wiki](https://develop.openfoam.com/Development/openfoam/-/wikis/) and [General wiki](http://wiki.openfoam.com/)
- [Governance](http://www.openfoam.com/governance/), [Governance Projects](https://www.openfoam.com/governance/projects)
- [Contacting OpenCFD](http://www.openfoam.com/contact/)

Copyright 2016-2024 OpenCFD Ltd
