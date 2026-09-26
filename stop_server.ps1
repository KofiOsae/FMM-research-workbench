$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pidFile = Join-Path $projectRoot "workbench.pid"
if (-not (Test-Path -LiteralPath $pidFile)) {
    Write-Host "No managed workbench server is recorded."
    exit 0
}
$serverPid = [int](Get-Content -LiteralPath $pidFile -Raw)
$process = Get-Process -Id $serverPid -ErrorAction SilentlyContinue
$expectedPython = Join-Path $projectRoot ".venv\Scripts\python.exe"
if ($process -and $process.Path -eq $expectedPython) {
    Stop-Process -Id $serverPid
    Remove-Item -LiteralPath $pidFile
    Write-Host "FMM Research Workbench stopped."
} else {
    Write-Host "The recorded process is no longer the workbench server."
}
