## OpenFOAM&reg; ThirdParty System Requirements

For building some particular third-party libraries from source,
the normal [OpenFOAM System Requirements][link openfoam-require]
may not be sufficient.

This is most notably the case for ParaView and/or QT compilation.
As duly noted in [BUILD][link third-build] and [README][link third-readme] information,
building ParaView from source tends to be the most difficult part of
any third-party compilation.

For general functionality, the paraview version distributed with
the operating system or a [binary package][download ParaView]
is likely [fully adequate for your needs][FAQ ParaView].


***Please help us with keeping the information here up-to-date and accurate.***

### Ubuntu

The full dependency list for building ParaView can be found from the
corresponding [debian/control][debian control] file.

A ***quick*** way to get these dependencies on to your system is to
use the paraview development package:
```
sudo apt install paraview-dev
```
Depending on your setup, the following subset may also be enough:
```
sudo apt install cmake qtbase5-dev qttools5-dev qttools5-dev-tools libqt5opengl5-dev libqt5x11extras5-dev libxt-dev
```


### openSUSE (eg, Leap-15.3)

The full dependency list for building ParaView can be found from the
corresponding [rpm spec][suse spec] file.

A ***quick*** way to get these dependencies on to your system is to
use the paraview development package:
```
sudo zypper install paraview-devel
```

Depending on your setup, the following subset may also be enough:
```
sudo zypper install Mesa-libEGL-devel
sudo zypper install libqt5-qtbase-devel libqt5-qtsvg-devel libqt5-qttools-devel libqt5-qtx11extras-devel
sudo zypper install libxcb-devel libXt-devel
```

For building paraview with ffmpeg support (components: avformat
avcodec avutil swscale), the packman repository may be required, with a
variety of additional packages:
```
sudo zypper install ffmpeg-3
sudo zypper install ffmpeg-3-libavcodec-devel
sudo zypper install ffmpeg-3-libavformat-devel
sudo zypper install ffmpeg-3-libavresample-devel
sudo zypper install ffmpeg-3-libavutil-devel
sudo zypper install ffmpeg-3-libswresample-devel
sudo zypper install ffmpeg-3-libswscale-devel
```


<!-- Quick links -->

[download ParaView]: https://www.paraview.org/download/
[debian control]: https://salsa.debian.org/science-team/paraview/-/blob/master/debian/control
[suse spec]: https://build.opensuse.org/package/view_file/science/paraview/paraview.spec
[FAQ ParaView]: https://discourse.paraview.org/t/i-want-to-visualize-my-openfoam-simulation-results-with-paraview-but-im-confused-which-version-should-i-use


<!-- OpenFOAM -->

[link openfoam-readme]: https://develop.openfoam.com/Development/openfoam/blob/develop/README.md
[link openfoam-config]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Config.md
[link openfoam-build]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md
[link openfoam-require]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/Requirements.md

---
Copyright 2019-2022 OpenCFD Ltd
