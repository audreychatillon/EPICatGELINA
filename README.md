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
In the EPICproject folder, edit `project.yaml` and uncomment line 10, adding the flag `--disable-mt` 
See an example: [project.yaml](https://github.com/audreychatillon/EPICatGELINA/blob/main/project.yaml)
In the EPICproject folder, create directory and subdirectories to write analysis root trees
```bash
mkdir output
cd output
mkdir analysis conversion simulation
```

## read raw FASTER data

To read FASTER data, configuration files must be provided:
 - `sample.pid`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/pid_files/sample_EPICproto_run24.pid)
 - `detector/detector.yaml`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/detector/detector_run24.yaml)
 - `ConfigEPIC.dat`[example](https://github.com/audreychatillon/EPICatGELINA/blob/main/config_files/ConfigEPIC_run24.dat)

The command `npconversion` is processing the function `EpicDetector::BuildRawEvent()` to fill data at raw level as defined in `EpicData`
To write raw TTree in output/conversion folder (see project.yaml) 
```bash
npconversion --input faster,sample.pid,/path/to/faster_file_num.fast --output root,EpicRawTree,raw_num.root
```

For on-line monitoring of raw histograms use the flag `--input-raw`.
If you want to monitor in a browser with localhost:8082:
```bash
npconversion --input faster,sample.pid,/path/to/FASTER/data/name_faster_file_num.fast --output root,8081
nponline --input-raw root,localhost:8081 --interface root,8082
```

To convert FASTER data and build physical event by applying calibration parameters, 
you should add in `project.yaml` in the `default flag` line `--calibration calibration.txt` which gives the path of all calibration files.
Then run the command `npanalysis` to process the function `EpicDetector::BuildPhysicalEvent()` to fill data at calibration level as defined in `EpicPhysics`.
To monitore the calibrated spectra use the flag `--input-phy`
```bash
npconversion --input faster,sample.pid,file.fast  --output root,8080
npanalysis --input root,localhost:8080 --output root,8081
nponline --input-raw root,localhost:8080 --input-phy root,localhost:8081 --interface root,8082
```

