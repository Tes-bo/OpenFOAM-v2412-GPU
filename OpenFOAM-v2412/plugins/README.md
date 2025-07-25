[[_TOC_]]

## OpenFOAM plugins

This directory is a location for additional OpenFOAM components or
tools to placed and have them built as part of the normal OpenFOAM
build process. It is assumed that each subdirectory contain an
appropriate `Allwmake` (or `Allwmake.override`) file.

The primary distinction between `modules` and `plugins` is that `modules` are
mainly maintained and released by OpenFOAM developers, whereas `plugins` are an
open and welcoming area for add-ons, predominantly maintained and driven by
OpenFOAM community members and groups.

### How to use

On the first use, you need to register the submodules, and then update them.
You can execute both steps for all the available submodules (including the
nested ones) as follows while you are at `$WM_PROJECT_DIR`:

```bash
cd $WM_PROJECT_DIR

git submodule update --init --recursive
```

Executing this single-line command clones all the submodules from their
respective repositories and prepares them for compilation. Note that you can
also make only a certain group of submodules ready by explicitly specifying the
requested submodules' names at the end of the command above. For example, if
you would like to use only the `turbulence-community` submodule, you specify:

```bash
git submodule update --init --recursive plugins/turbulence-community
```

You can display information about the status of submodules as follows:

```bash
git submodule status --recursive
```

An easy way to see which submodules are actually in use:

```bash
cat .gitmodules
```

Which will reveal content resembling the following:
```
[submodule "xyz"]
    path = plugins/xyz
    url = ...
...
```

If you need to remove a specific submodule or wish to restart the process,
you can simply carry out the task as follows:

```bash
git submodule deinit plugins/turbulence-community
```

This command deregisters the specified submodule and clears the
`plugins/turbulence-community` directory.

A quick overview of `git submodules` can be found in this
[*blog*][blog git-submodule] with full details in the
[*manpage*][man git-submodule].

### Build locations

Any individual _plugin_ will normally also be able to exist outside of
the plugins directory structure and will typically build into user
locations (`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`).

When compiled from the top-level OpenFOAM `Allwmake` or the
`plugins/Allwmake`, they should build into OpenFOAM project locations
(`$FOAM_APPBIN`, `$FOAM_LIBBIN`). This can be adjusted by
supplying an alternative `-prefix=` to the corresponding Allwmake
command.

| Command    | Install location |
|------------|------------------|
| ./Allwmake -prefix=user | `$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN` |
| ./Allwmake -prefix=group | `$FOAM_SITE_APPBIN`, `$FOAM_SITE_LIBBIN` |
| ./Allwmake -prefix=openfoam | `$FOAM_APPBIN`, `$FOAM_LIBBIN` |
| ./Allwmake -prefix=/some/pathname | `/some/pathname/bin`, `/some/pathname/lib` |

### Documentation (doxygen)

To build the doxygen information for the components, it is also
necessary to link the directories to the doc/ subdirectory.
This is a purely manual operation.

### Developer information

#### Build locations

To accomodate building into various locations, the plugin code should
be adapted with the following changes:

- ***Make/files***
   ```
   ...
   EXE = $(FOAM_MODULE_APPBIN)/someExecutable

   LIB = $(FOAM_MODULE_LIBBIN)/libSomeLibrary
   ```

- `Make/options` should include this
  ```
  include $(GENERAL_RULES)/module-path-user
  ...
  ```

The following changes to `Make/options` are universally applicable
(ie, work with older or other versions of OpenFOAM), but more verbose.

- `Make/options` with the following
  ```
  sinclude $(GENERAL_RULES)/module-path-user

  /* Failsafe - user locations */
  ifeq (,$(FOAM_MODULE_APPBIN))
  FOAM_MODULE_APPBIN = $(FOAM_USER_APPBIN)
  endif
  ifeq (,$(FOAM_MODULE_LIBBIN))
  FOAM_MODULE_LIBBIN = $(FOAM_USER_LIBBIN)
  endif
  ...
  ```

<!-- General Information -->

[man git-submodule]:  https://git-scm.com/docs/git-submodule
[blog git-submodule]: http://blog.joncairns.com/2011/10/how-to-use-git-submodules/

---
