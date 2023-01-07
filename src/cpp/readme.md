# c++
This is a C++ project using CMake and Conan.

## status
This is just getting started.

## summary
This is a modern implementation of MoRPE.  The intent is to match the API and numerical functionality of the C# version.

## setup and build

### windows

#### setup
Install VS 2022 (not free) or Build Tools for VS 2022 (free).
  * TO DO:  Do they need to check boxes for C++ support?
  * TO DO:  Do they need to install a particular Windows SDK?

Install conan, cmake.

```pwsh
# Elevated prompt
choco install --yes conan
choco install --yes cmake --installargs 'ADD_CMAKE_TO_PATH=System'
choco install --yes make
```

#### build
From a prompt of type `Developer PowerShell for VS 2022`:

```pwsh
./build.ps1
```

### linux (fedora 34)

#### setup
TO DO

#### build
TO DO
