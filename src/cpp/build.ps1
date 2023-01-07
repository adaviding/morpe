# The Windows Build:  Execute this from a prompt of type `Developer PowerShell for VS 2022`.

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$conan_profile = "$PSScriptRoot/conan/profiles/win64"

$build_types = @("Debug", "RelWithDebInfo")
$subprojects = @("morpe", "tests/morpe")

foreach ($build_type in $build_types)
{
    $build_folder = "$PSScriptRoot/cmake-build-$($build_type.ToLower() )"

    foreach ($subproject in $subprojects)
    {
        $proj_folder = "$PSScriptRoot/$subproject"
        $proj_build_folder = "$build_folder/$subproject"

        Write-Host "[BEGIN] $build_type $subproject conan install"

        mkdir -Force "$proj_build_folder"
        conan install "$proj_folder" --install-folder "$proj_build_folder" --build=outdated --update --profile $conan_profile -s build_type=$build_type

        Write-Host "[END] $build_type $subproject conan install"
    }

    # ------------------------------------------------------------

    Write-Host "[BEGIN] $build_type configure"

    cmake -DCMAKE_BUILD_TYPE=$build_type -S "$PSScriptRoot" -B "$build_folder"

    Write-Host "[END] $build_type configure"

    # ------------------------------------------------------------

    Write-Host "[BEGIN] $build_type build"

    make --directory="$build_folder"

    Write-Host "[END] $build_type build"
}
