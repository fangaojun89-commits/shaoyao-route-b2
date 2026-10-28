param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Script,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ScriptArgs
)

$ProjectRoot = Split-Path -Parent $PSScriptRoot
$RHome = Join-Path $PSScriptRoot 'R-4.6.1'
$RScript = Join-Path $RHome 'bin\Rscript.exe'
$RLibrary = Join-Path $PSScriptRoot 'R-library\4.6'
$RenvRoot = Join-Path $PSScriptRoot 'renv-cache'

if (-not (Test-Path -LiteralPath $RScript)) {
    throw "Rscript not found at $RScript"
}

$env:R_HOME = $RHome
$env:R_LIBS_USER = $RLibrary
$env:RENV_PATHS_ROOT = $RenvRoot
$env:RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE'
$env:R_PROFILE = 'NUL'
$env:R_PROFILE_USER = 'NUL'
$env:R_ENVIRON = 'NUL'
$env:R_ENVIRON_USER = 'NUL'
$env:LANG = 'English_United States.utf8'
$env:LC_ALL = 'English_United States.utf8'
$env:LC_CTYPE = 'English_United States.utf8'

& $RScript --vanilla $Script @ScriptArgs
exit $LASTEXITCODE
