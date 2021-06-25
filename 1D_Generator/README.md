### This requires the library eigen to be in this subdirectory 
Add the submodule by doing 

```
git submodule add https://gitlab.com/libeigen/eigen.git
```
### Build with: 

```
make -j4
```

### Set the datagenerator parameters through the config.ini

The config.ini file follows the python `ConfigParser` library syntax. Utilize that for automating data generation. 
