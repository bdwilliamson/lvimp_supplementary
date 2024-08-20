:: run all correlated-data simulations for lvim paper
@echo off
setlocal enabledelayedexpansion

:: get R version for the machine (defaults to the latest version)
For /F "Delims=" %%0 In ('where /r "C:\Program Files\R" Rscript.exe') do set scriptpath="%%~0"
echo Using R executable !scriptpath!

:: set up static args
set B=0
set outcome_type=binary
set p=10
set num_timepoints=4
set nreps=1000
set W=0.5
set add_in=1
set leave_out=1
set simple_model=1

:: send output to directory
set outdir="%cd%\rout"
if not exist %outdir% mkdir %outdir%

:: loop over sample sizes
for %%N in (100 250 500 1000 5000 10000) do (
    echo Running n = %%N, outcome = !outcome_type!, cor between = !B!, cor within = !W!
    set this_outfile=!outdir!\output_!outcome_type!_cb_!B!_cw_!W!_n_%%N_simple_!simple_model!.out
    
    !scriptpath! investigate_cross_sectional_performance.R --outcome-type !outcome_type! --cor-between !B! --cor-within !W! --n %%N --p !p! --num-timepoints !num_timepoints! --nreps-total !nreps! --nreps-per-job !nreps! --simple-model !simple_model! 1>!this_outfile! 2>&1
)