# Compile every program of the repository into the "run" folder.
# Usage (from the repository root):   .\build.ps1
# Requires gfortran in the PATH.

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
New-Item -ItemType Directory -Force -Path "$root\run" | Out-Null
Push-Location "$root\run"

function Build([string]$out, [string[]]$sources) {
    Write-Host "Compiling $out ..."
    $paths = $sources | ForEach-Object { "..\src\$_" }
    & gfortran -O2 -fdefault-real-8 @paths -o $out
    if ($LASTEXITCODE -ne 0) { Pop-Location; throw "Compilation failed: $out" }
}

Build "hydro"   @("hydro_solvers.f90",   "hydro.f90")
Build "hydro2d" @("hydro2d_solvers.f90", "hydro2d.f90")
Build "jet2d"   @("hydro2d_solvers.f90", "jet2d.f90")
Build "mhd"     @("mhd_solvers.f90",     "mhd.f90")
Build "orszag"  @("mhd2d_solvers.f90",   "orszag.f90")
Build "jetmhd"  @("mhd2d_solvers.f90",   "jetmhd.f90")

Pop-Location
Write-Host "Done. Executables are in the run folder."
