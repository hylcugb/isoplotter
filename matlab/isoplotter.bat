@echo off

set srcdir=%~dp0\Matlab_source_files
octave -p %srcdir% -qf %srcdir%\CommandLineInterface.m %*