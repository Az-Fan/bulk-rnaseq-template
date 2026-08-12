param(
    [Parameter(Mandatory=$true, Position=0)][string]$CountsPath,
    [Parameter(Mandatory=$true, Position=1)][string]$ProjectId
)
$ErrorActionPreference = 'Stop'
$resolved = (Resolve-Path -LiteralPath $CountsPath).Path
$linuxSource = (wsl wslpath -a $resolved).Trim()
if (-not $linuxSource) { throw "Unable to convert the input path for WSL: $resolved" }
$project = "projects/$ProjectId"
wsl bash -lc 'cd "$(dirname "$(readlink -f internal/windows/Import-Counts.ps1 2>/dev/null || pwd)")" >/dev/null 2>&1 || true; root="$HOME/work/bulk-rnaseq-v4"; cd "$root" && "$root/.pixi/bin/python" internal/cli/pipeline.py import-counts --project "$1" --source "$2"' bash $project $linuxSource
if ($LASTEXITCODE -ne 0) { throw "Counts import failed with exit code $LASTEXITCODE" }
Read-Host 'Import completed. Press Enter to close'
