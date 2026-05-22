# Install all 15 alignment tools via conda on Windows
# Run this script in PowerShell (as Administrator recommended)

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Installing 14 WGBS alignment tools" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$EnvName = "bs-benchmark"

# Check conda
if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: conda not found. Install Miniconda first." -ForegroundColor Red
    Write-Host "Download: https://docs.conda.io/en/latest/miniconda.html" -ForegroundColor Yellow
    exit 1
}

# Create environment
$envs = conda env list
if ($envs -notcontains $EnvName) {
    Write-Host "Creating conda environment: $EnvName" -ForegroundColor Green
    conda create -y -n $EnvName python=3.9
}

# Activate
conda activate $EnvName

# Python dependencies
pip install click pyyaml pandas numpy biopython matplotlib seaborn rich

# Bioconda tools
Write-Host "Installing aligners via conda..." -ForegroundColor Green
conda install -y -c bioconda -c conda-forge `
    bwameth `
    bsbolt `
    bsmap `
    bismark `
    hisat2 `
    bs-seeker2 `
    samtools `
    bedtools

# hisat-3n (may need to build from source on Windows)
Write-Host "`nNote: hisat-3n may need manual installation on Windows." -ForegroundColor Yellow
Write-Host "See: https://github.com/DaehwanKimLab/hisat-3n" -ForegroundColor Yellow

# Pip-based tools
Write-Host "Installing pip-based tools..." -ForegroundColor Green
pip install walt-bio
pip install batmeth2

Write-Host "`nInstallation complete!" -ForegroundColor Green
Write-Host "Run: conda activate $EnvName" -ForegroundColor Cyan

# Check
Write-Host "`n=== Tool Availability Check ===" -ForegroundColor Cyan
$tools = @("bwameth.py", "bsbolt", "bsmap", "bismark", "bs_seeker2-align.py")
foreach ($tool in $tools) {
    $loc = Get-Command $tool -ErrorAction SilentlyContinue
    if ($loc) {
        Write-Host "[OK] $tool" -ForegroundColor Green
    } else {
        Write-Host "[MISSING] $tool" -ForegroundColor Red
    }
}
