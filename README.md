###Project using the epic plugin with NPTool V4

##For plugin management:
[Plugins Management](https://nptool.in2p3.fr/manual-v4/plugins-management/)
```bash
nptool --install epic
# This will
# - install the epic plugin in ~/.local/nptool/default/
# - copy all EpicXXX.h in ~/.local/nptool/default/include/ 
nptool --update epic
# This will update the last features of the plugin
```

## Create and manage your project with the epic plugin

To create your project:
[Project Management](https://nptool.in2p3.fr/manual-v4/project-management/)
```bash
cd /folder/of/your/nptoolV4/projects/
nptool --new-project EPICproject
# This will
# - update the project.list in ~/.local/nptool/default/
# - create the folder EPICproject in /folder/of/your/nptoolV4/projects/
```

