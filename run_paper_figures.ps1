param(
    [int[]]$Figures = @(2,3,4,5,6),
    [string]$Matlab = 'C:/Program Files/MATLAB/R2022b/bin/matlab.exe',
    [switch]$RunDimensionScaling
)
$ErrorActionPreference = 'Stop'
$taskRoot = $PSScriptRoot
if (($Figures | Where-Object { $_ -lt 1 -or $_ -gt 6 }).Count -gt 0) {
    throw 'Figure numbers must be between 1 and 6.'
}
if (($Figures -contains 1) -and -not $RunDimensionScaling) {
    throw 'Figure 1 requires -RunDimensionScaling; the default run excludes it.'
}
$taskLogDirectory = Join-Path $taskRoot 'results/logs'
New-Item -ItemType Directory -Path $taskLogDirectory -Force | Out-Null
$taskMatlabRoot = $taskRoot.Replace('\','/').Replace("'","''")
foreach ($taskFigure in $Figures) {
    $taskLog = Join-Path $taskLogDirectory ('paper_fig{0:D2}.log' -f $taskFigure)
    $taskAllow = if ($RunDimensionScaling) { 'true' } else { 'false' }
    $taskExpression = "cd('$taskMatlabRoot'); addpath(genpath(pwd),'-begin'); setenv('LC_ALL',''); setenv('LC_CTYPE',''); setenv('LANG',''); setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE',''); set(groot,'defaultFigureVisible','off'); maxNumCompThreads(1); main_run_paper_figures(struct('figures',$taskFigure,'run_dimension_scaling',$taskAllow,'continue_on_error',false));"
    Write-Output ('Computing paper figure {0}; log: {1}' -f $taskFigure,$taskLog)
    & $Matlab -singleCompThread -batch $taskExpression -logfile $taskLog
    if ($LASTEXITCODE -ne 0) {
        throw ('MATLAB failed for figure {0}; inspect {1}' -f $taskFigure,$taskLog)
    }
}
