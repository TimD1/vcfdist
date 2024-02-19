Below are some common errors you may run into when attempting to install or run `vcfdist`. 

If you encounter any problems not listed below, please open a new [issue](https://github.com/TimD1/vcfdist/issues).

### Failure to compile on Mac
On Mac devices, `g++` is aliased to `clang`, which is currently not supported. Please install `g++` through [homebrew](https://brew.sh/).

### Errors with libstdc++fs
Although `std::filesystem` is officially part of the C++17 standard library, g++ < 9.1 and llvm < 9.0 require linking to `-lstdc++fs` during compilation. If you're using g++ < 9.1, please add `-lstdc++fs` to `LDFLAGS` in `src/Makefile` before running `make`.

### Errors while loading shared libraries
You may encounter the following error while trying to run vcfdist:
```bash
vcfdist: error while loading shared libraries: libhts.so.2to3part12: cannot open shared object file: No such file or directory
```
This occurs because the HTSlib library cannot be found. If you have installed HTSlib in `/usr/local`, add `/usr/local/lib` to your linker's library path:
```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
Optionally, you can add this export to your [`.bashrc`](https://www.digitalocean.com/community/tutorials/bashrc-file-in-linux) file so that HTSlib is always discoverable:
```bash
echo "export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
```
Please note the escape `\` character before `$`, and use of `>>` to append to your `.bashrc` (NOT `>`, which will overwrite the file).