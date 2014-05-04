@echo off

set srcdir=%~dp0\Matlab_source_files

set cmd="CommandLineInterface %*"
octave --traditional -q -p %srcdir% --eval %cmd%