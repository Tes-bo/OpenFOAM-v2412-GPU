# Updating

Most of the code in this repository is managed as part of Kitware's
[fork of VisIt][]. This branch pulls from the `for/paraview` branch in that
repository. Changes to readers and the core code should happen there first.

The code may then be updated using the update tools for the [readers][] and
for the [avt core library][]. Simply run the relevant tools and create a merge
request to this project with the updates.

Once this repository is updated, ParaView needs its `Utilities/VisItBridge`
submodule updated to point to the new version of this repository.

[fork of VisIt]: https://gitlab.kitware.com/third-party/visit
[readers]: scripts/update-databases.sh
[avt core library]: scripts/update-avt.sh
