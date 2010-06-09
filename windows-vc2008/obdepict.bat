@echo off
rem Displays up to 64 molecules from any file with a chemical extension.
rem Removes hydrogens. Generates 2D coordinates if none present (-a2 is for cml)
babel %1 %temp%\temp.svg -d -a2 -xN64
firefox %temp%\temp.svg
