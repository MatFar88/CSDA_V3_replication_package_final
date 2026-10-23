#requires -Version 5.1
<#
Windows x64 reporting bootstrap. No installed Python, py, pip, uv, winget or
administrator account is required. Network access is requested explicitly.
Runtime and cache live under .runtime; scientific code and data are not changed.
The executable download is a pinned uv wheel from PyPI, verified against the
SHA-256 in PyPI's HTTPS metadata. Only uv.exe is extracted, to a fixed filename.
#>
[CmdletBinding()]
param([switch]$Yes)
Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
$Root = Split-Path -Parent $PSScriptRoot
$Runtime = Join-Path $Root '.runtime'
$Venv = Join-Path $Root '.venv'
$Python = Join-Path $Venv 'Scripts\python.exe'
$Probe = Join-Path $PSScriptRoot 'check_python_environment.py'
$Requirements = Join-Path $Root 'requirements.txt'
$UvVersion = '0.12.20'
$Uv = Join-Path $Runtime ('uv-' + $UvVersion + '\uv.exe')

function Invoke-Native {
    param([string]$Exe, [string[]]$ArgumentList)
    & $Exe @ArgumentList
    if ($LASTEXITCODE -ne 0) {
        throw ('Command failed (exit {0}): {1}' -f $LASTEXITCODE, $Exe)
    }
}

function Get-Uv {
    $dir = Split-Path -Parent $Uv
    New-Item -ItemType Directory -Path $dir -Force | Out-Null
    $receipt = Join-Path $dir 'download_receipt.json'
    if ((Test-Path -LiteralPath $Uv -PathType Leaf) -and
        (Test-Path -LiteralPath $receipt -PathType Leaf)) {
        $previous = Get-Content -LiteralPath $receipt -Raw | ConvertFrom-Json
        $actual = (Get-FileHash -LiteralPath $Uv -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($actual -eq $previous.executable_sha256) {
            return
        }
        throw 'The cached setup executable has changed. Rename .runtime and rerun setup.'
    }
    $metadataUrl = 'https://pypi.org/pypi/uv/' + $UvVersion + '/json'
    Write-Host ('Fetching pinned setup tool uv ' + $UvVersion + ' from PyPI ...')
    $meta = Invoke-RestMethod -Uri $metadataUrl -TimeoutSec 120
    $wheels = @($meta.urls | Where-Object {
        $_.packagetype -eq 'bdist_wheel' -and
        $_.filename -like ('uv-' + $UvVersion + '-*-win_amd64.whl')
    })
    if ($wheels.Count -ne 1) {
        throw 'Cannot identify the pinned Windows x64 uv wheel in PyPI metadata.'
    }
    $wheel = $wheels[0]
    $url = [Uri]$wheel.url
    if ($url.Scheme -ne 'https' -or $url.Host -ne 'files.pythonhosted.org') {
        throw 'Unexpected executable download host; download refused.'
    }
    $expected = [string]$wheel.digests.sha256
    if ($expected -notmatch '^[0-9a-fA-F]{64}$') {
        throw 'No valid SHA-256 was supplied by PyPI.'
    }
    $archive = Join-Path $dir 'uv-download.whl'
    Invoke-WebRequest -Uri $url.AbsoluteUri -OutFile $archive -UseBasicParsing -TimeoutSec 300
    $actual = (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash
    if ($actual -ne $expected) {
        Remove-Item -LiteralPath $archive -Force
        throw 'The downloaded setup tool failed its SHA-256 check.'
    }
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip = [System.IO.Compression.ZipFile]::OpenRead($archive)
    try {
        $entries = @($zip.Entries | Where-Object {
            $_.FullName -eq 'uv.exe' -or $_.FullName.EndsWith('/uv.exe')
        })
        if ($entries.Count -ne 1) { throw 'The verified wheel must contain exactly one uv.exe.' }
        $src = $entries[0].Open()
        try {
            $dst = [System.IO.File]::Create($Uv)
            try { $src.CopyTo($dst) } finally { $dst.Dispose() }
        } finally { $src.Dispose() }
    } finally { $zip.Dispose() }
    $info = [ordered]@{
        version = $UvVersion
        metadata_url = $metadataUrl
        wheel_url = $url.AbsoluteUri
        wheel_sha256 = $expected.ToLowerInvariant()
        executable_sha256 = (Get-FileHash -LiteralPath $Uv -Algorithm SHA256).Hash.ToLowerInvariant()
        downloaded_at_utc = [DateTime]::UtcNow.ToString('o')
    }
    $info | ConvertTo-Json | Set-Content -LiteralPath $receipt -Encoding UTF8
    Remove-Item -LiteralPath $archive -Force
    Invoke-Native -Exe $Uv -ArgumentList @('--version')
}

$exitCode = 1
$transcribing = $false
$oldLocation = Get-Location
# These settings are process-local; all are restored even when dot-invoked.
$settings = @{
    'UV_CACHE_DIR' = (Join-Path $Runtime 'cache')
    'UV_PYTHON_INSTALL_DIR' = (Join-Path $Runtime 'python')
    'UV_PYTHON_BIN_DIR' = (Join-Path $Runtime 'bin')
    'UV_PYTHON_INSTALL_BIN' = '0'
    'UV_PYTHON_INSTALL_REGISTRY' = '0'
    'UV_NO_MODIFY_PATH' = '1'
    'UV_NO_CONFIG' = '1'
    'UV_SYSTEM_CERTS' = '1'
    'UV_LINK_MODE' = 'copy'
    'UV_PYTHON_DOWNLOADS' = 'automatic'
    'MPLCONFIGDIR' = (Join-Path $Runtime 'matplotlib')
}
$previousEnv = @{}
$oldTls = [Net.ServicePointManager]::SecurityProtocol
try {
    if ($env:OS -ne 'Windows_NT') { throw 'This bootstrap is for Windows only.' }
    $arch = $env:PROCESSOR_ARCHITECTURE
    if ($env:PROCESSOR_ARCHITEW6432) { $arch = $env:PROCESSOR_ARCHITEW6432 }
    if ($arch -ne 'AMD64') {
        throw 'This bootstrap targets x64 Windows. Other platforms use their documented setup route.'
    }
    foreach ($key in $settings.Keys) {
        $previousEnv[$key] = [Environment]::GetEnvironmentVariable($key, 'Process')
        [Environment]::SetEnvironmentVariable($key, $settings[$key], 'Process')
    }
    Set-Location -LiteralPath $Root
    New-Item -ItemType Directory -Path $Runtime -Force | Out-Null
    $logDir = Join-Path $Root 'output\setup_logs'
    New-Item -ItemType Directory -Path $logDir -Force | Out-Null
    $log = Join-Path $logDir ('setup_' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff') + '.txt')
    try { Start-Transcript -Path $log -ErrorAction Stop | Out-Null; $transcribing = $true }
    catch { Write-Host 'Transcript unavailable; console messages are still shown.' }

    $ready = $false
    if (Test-Path -LiteralPath $Python -PathType Leaf) {
        & $Python -I $Probe --base-only --require-venv
        if ($LASTEXITCODE -ne 0) {
            throw 'Existing .venv is not a usable 64-bit Python 3.13 virtual environment. Rename it and rerun setup; it will not be deleted automatically.'
        }
        & $Python -I $Probe --require-venv
        $ready = ($LASTEXITCODE -eq 0)
    } elseif (Test-Path -LiteralPath $Venv) {
        throw 'An incomplete .venv already exists. Rename that directory and rerun setup; no existing directory is deleted automatically.'
    }

    if (-not $ready) {
        Write-Host ''
        Write-Host 'Setup requires internet access to download missing reporting software.'
        Write-Host 'It uses PyPI and, if needed, an Astral-managed CPython 3.13 runtime.'
        Write-Host 'Downloads stay in .runtime and .venv. No administrator account,'
        Write-Host 'global PATH change, registry registration, or MATLAB fitting is requested.'
        if (-not $Yes) {
            $answer = Read-Host 'Continue with these downloads? [y/N]'
            if ($answer -notmatch '^(?i:y|yes|s|si)$') {
                throw 'Setup cancelled. The saved HTML report can be viewed without Python.'
            }
        }
        [Net.ServicePointManager]::SecurityProtocol = $oldTls -bor [Net.SecurityProtocolType]::Tls12
        Get-Uv
        if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
            Invoke-Native -Exe $Uv -ArgumentList @('python', 'install', '3.13', '--no-bin', '--no-registry', '--no-config')
            Invoke-Native -Exe $Uv -ArgumentList @('venv', '--managed-python', '--python', '3.13', '--no-config', $Venv)
        }
        # Never install into a global interpreter, even after a partial setup.
        Invoke-Native -Exe $Python -ArgumentList @('-I', $Probe, '--base-only', '--require-venv')
        Invoke-Native -Exe $Uv -ArgumentList @('pip', 'install', '--python', $Python,
            '--default-index', 'https://pypi.org/simple', '--only-binary', ':all:',
            '--no-config', '--requirements', $Requirements)
        Invoke-Native -Exe $Python -ArgumentList @('-I', $Probe, '--require-venv')
    } else {
        Write-Host 'Existing .venv is ready. No software is downloaded or reinstalled.'
    }
    Invoke-Native -Exe $Python -ArgumentList @((Join-Path $Root 'run_replication.py'), 'CHECK')
    Write-Host ''
    Write-Host 'SETUP COMPLETE. Use RUN_REPLICATION("REPORT") in MATLAB,'
    Write-Host 'or double-click RUN_REPORT_WINDOWS.cmd. No estimator has been fitted.'
    $exitCode = 0
} catch {
    Write-Host ''
    Write-Host ('SETUP STOPPED: ' + $_.Exception.Message) -ForegroundColor Red
    Write-Host 'No fitting command will run. For download errors check network/proxy'
    Write-Host 'access or institutional software policy; do not disable TLS verification.'
    $exitCode = 1
} finally {
    if ($transcribing) { Stop-Transcript | Out-Null }
    [Net.ServicePointManager]::SecurityProtocol = $oldTls
    foreach ($key in $previousEnv.Keys) {
        [Environment]::SetEnvironmentVariable($key, $previousEnv[$key], 'Process')
    }
    Set-Location -LiteralPath $oldLocation.Path
}
exit $exitCode
