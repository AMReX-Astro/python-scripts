These are some standalone python routines for interacting with
FBoxLib/AMReX data.

They require the FBoxLib library:

```
git clone https://github.com/AMReX-Codes/FBoxLib
```

You need to build a python module first to load the data.  Typing
```
make
```
should do this.  It requires the `f2py` tool, which comes with NumPy.

Some documentation for these routines is available in the MAESTRO
User's Guide.

Note: it is recommended that you use `yt` instead of these routines,
since `yt` can work directly with the AMR structure of the data.
