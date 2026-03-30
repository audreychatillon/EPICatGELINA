###Project using the epic plugin with NPTool V4

**epic** is a plugin for **NPTool V4**, inspired (copied) from the fission chamber plugin.
We try to avoid hard-coding.
It can be used in project with several fission chambers, having different number of anodes, with different actinides 


## epic plugin management

For plugin management:
[Plugins Management](https://nptool.in2p3.fr/manual-v4/plugins-management/)
```bash
nptool --install epic
# This will
# - install the epic plugin in ~/.local/nptool/default/
# - copy all EpicXXX.h in ~/.local/nptool/default/include/ 
nptool --update epic
# This will update the last features of the plugin
```

## Create your project with the epic plugin

To create your project:
[Project Management](https://nptool.in2p3.fr/manual-v4/project-management/)
```bash
cd /folder/of/your/nptoolV4/projects/
nptool --new-project EPICproject
# This will
# - update the project.list in ~/.local/nptool/default/
# - create the folder EPICproject in /folder/of/your/nptoolV4/projects/
```

You must compile (after each update of the plugin), in your project folder:
```bash
mkdir build install
cmake -B build -DCMAKE_INSTALL_PREFIX=install ./
make -C build/ install
```

## read raw FASTER data

To read FASTER data, configuration files must be provided:
 - `sample.pid`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/pid_files/sample_EPICproto_run24.pid)
 - `detector/detector.yaml`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/detector/detector_run24.yaml)
 - `ConfigEPIC.dat`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/config_files/ConfigEPIC_run24.dat)


