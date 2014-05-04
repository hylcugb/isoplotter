@echo off

set srcdir=%~dp0\Matlab_source_files

octave -p %srcdir% -q %srcdir%\CommandLineInterface.m %*