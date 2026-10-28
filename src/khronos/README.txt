Khronos OpenCL C API headers (cl.h, cl_platform.h, cl_version.h), version
v2022.01.04 from https://github.com/KhronosGroup/OpenCL-Headers, licensed under
the Apache License 2.0 (see the header of each file).

They are used on macOS in place of the OpenCL framework headers, which mark
every OpenCL function deprecated; the binaries still link -framework OpenCL.
The only change from the upstream files: the #include <CL/...> lines of cl.h
and cl_platform.h use "..." so that no include path is needed.
