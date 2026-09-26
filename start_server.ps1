param([switch]$NoBrowser)
$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonExe = Join-Path $projectRoot ".venv\Scripts\python.exe"
$pythonWindowless = Join-Path $projectRoot ".venv\Scripts\pythonw.exe"
$requirementsFile = Join-Path $projectRoot "requirements.txt"
$pidFile = Join-Path $projectRoot "workbench.pid"
$stdoutFile = Join-Path $projectRoot "server_stdout.log"
$stderrFile = Join-Path $projectRoot "server_stderr.log"

Set-Location $projectRoot
$createdEnvironment = -not (Test-Path -LiteralPath $pythonExe)
if ($createdEnvironment) {
    py -3.12 -m venv .venv
    if ($LASTEXITCODE -ne 0) { throw "Python 3.12 is required." }
    & $pythonExe -m pip install -r $requirementsFile --disable-pip-version-check
    if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }
}
& $pythonExe -m pip check
if ($LASTEXITCODE -ne 0) { throw "The Python environment is incomplete. Delete .venv and start again." }

$healthy = $false
try {
    $reply = Invoke-WebRequest -UseBasicParsing -Uri "http://127.0.0.1:8765/health" -TimeoutSec 2
    $healthy = $reply.StatusCode -eq 200
} catch {}

if (-not $healthy) {
    # Start-Process can fail on hosts that expose both Path and PATH environment
    # keys. ProcessStartInfo uses the current environment without rebuilding that
    # case-insensitive dictionary. app.py writes runtime errors to workbench.log.
    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    # Use python.exe so the returned PID remains the long-lived server process.
    # pythonw.exe in a virtual environment can spawn the base interpreter and
    # exit, leaving a stale PID file that Stop cannot use.
    $startInfo.FileName = $pythonExe
    $startInfo.Arguments = "app.py"
    $startInfo.WorkingDirectory = $projectRoot
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $process = [System.Diagnostics.Process]::Start($startInfo)
    Set-Content -LiteralPath $pidFile -Value $process.Id
    for ($attempt = 0; $attempt -lt 20 -and -not $healthy; $attempt++) {
        Start-Sleep -Milliseconds 250
        try {
            $reply = Invoke-WebRequest -UseBasicParsing -Uri "http://127.0.0.1:8765/health" -TimeoutSec 2
            $healthy = $reply.StatusCode -eq 200
        } catch {}
    }
}

if (-not $healthy) {
    throw "The solver did not start. Read server_stderr.log in this folder."
}

if (-not $NoBrowser) { Start-Process "http://127.0.0.1:8765/" }
Write-Host "FMM Research Workbench is running in the background."
Write-Host "Open http://127.0.0.1:8765/ at any time."
