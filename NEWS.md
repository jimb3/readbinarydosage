# readbinarydosage 1.0

## Initial package set up

These are the things I find useful to do when setting up a package that
I plan to release to CRAN that use Rcpp

* Added a `NEWS.md` file to track changes to the package.
* Added readbinarydosage with lines to properly handle C++ functions
* Created R function to call C++ code
* Documented R function
* Modified skeleton C++ code to return a value of 1
* Modified DESCRIPTION with appropriate entries
* Deleted NAMESPACE
* Deleted files in man directory
* Ran roxygen2::roxygenise() to created documentation and new NAMESPACE file

After doing the above, devtools::check() reports no errors warnings or notes

## Adding testing

These tasks can be added during the initial set up and run. These are for
packages using Rcpp

* Added AppVeyor support
* Added Travis support
* Added first test - needed for code coverage
* Added support for code coverage

## Adding functionality

### Testing

Once functionality items are added these functions should be run. They will
be required before submitting to CRAN. It is better to have them checked
as needed than when you are ready to release.

* R CMD check
* spell_check
* release_checks
* check_rhub
* check_win_devel
* check_win_release

### Opening the file

Simple routines to open the binary dosage file in C++. This is called
from R. It is surprising to see how many checks are needed to verify that
the file was opened successfully.

* Can the file be opened (does it exist)
* Can the header be read
* Is it a binary dosage file
* Is the format variable correctly stored
* Is the format number valid
* Is the subformat valid (2 tests)

After adding the code to open files the testing produced a couple of notes. One was about the case of the title needing to be in title case. The other was about the license.
Removed LICENSE file and changed case of title. Still get a note about being a new submission but that doesn't seem to be important.