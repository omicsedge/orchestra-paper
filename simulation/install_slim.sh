#!/bin/bash
# This script downloads the source archive of SLiM,
# extracts it, creates a build directory and builds
# the command-line utilities for slim and eidos, and
# also the SLiMgui IDE. It then installs them to
# /usr/bin, and installs the FreeDesktop files to the
# appropriate places for desktop integration.

# Bryce Carson, @bryce-carson on GitHub, 12/27/2020
# This script is public domain.  Use freely.
# Please report issues and submit pull requests
# against the SLiM-Extras GitHub repo, tagging Bryce.

# We need superuser privileges.
if [ "$(id -u)" -ne 0 ]; then
        echo 'This script must be run by root' >&2
	echo "Invoke the script with sudo."
        exit 1
fi

# Test that build requirements are satisfied.
# If any one requirement is unmet we say which.
unset cmakeinstalled; dpkg-query -s cmake 2>/dev/null | grep -q ^"Status: install ok installed"$; cmakeinstalled=$?
#unset qmakeinstalled; dpkg-query -s qt5-qmake 2>/dev/null | grep -q ^"Status: install ok installed"$; qmakeinstalled=$?
#unset qtchooserinstalled; dpkg-query -s qtchooser 2>/dev/null | grep -q ^"Status: install ok installed"$; qtchooserinstalled=$?
#unset qtbase5devinstalled; dpkg-query -s qtbase5-dev 2>/dev/null | grep -q ^"Status: install ok installed"$; qtbase5devinstalled=$?
#unset curlinstalled; dpkg-query -s curl 2>/dev/null | grep -q ^"Status: install ok installed"$; curlinstalled=$?
unset wgetinstalled; dpkg-query -s wget 2>/dev/null | grep -q ^"Status: install ok installed"$; wgetinstalled=$?

#[[ $qmakeinstalled == 0 && $qtchooserinstalled == 0 && $qtbase5devinstalled == 0 ]] || printf "All of: qt5-qmake, qtchooser, and qtbase5-dev must be installed.\nInstall the Qt5 requirements with 'sudo apt install qtbase5-dev qtchooser qt5-qmake'.\nInstalling these packages ensures all build\nand runtime requirements are satisfied.\n"
[[ $cmakeinstalled == 0 ]] || printf "cmake is not installed.\nInstall it with 'sudo apt install cmake'.\n"
#[[ $curlinstalled == 0 || $wgetinstalled == 0 ]] || printf "Neither curl nor wget are installed.\nInstall either with one of:\n'sudo apt install wget', OR\n'sudo apt install curl'.\n"
#[[ $qtchooserinstalled == 0 && $qtbase5devinstalled == 0 && $qmakeinstalled == 0 && $cmakeinstalled == 0 ]] || exit #Exit if qtchooser, qtbase5-dev, qt5-qmake, or cmake are not installed. If neither curl nor wget are installed, we exit later.

cd /tmp || { printf "The Filesystem Hierarchy-standard directory /tmp does not exist.\nResolve the issue by creating that directory; inspect this script, and your system,\nas other issues may exist." ; exit ;}
if [[ -f SLiM.zip ]]; then
	printf "/tmp/SLiM.zip already exists! Please delete the file to use this script.\n"
	exit
fi

if [[ -d /tmp/BUILD/ ]]; then
	printf "/tmp/BUILD/ already exists! Please move that directory to use this script.\n"
	exit
fi

if [[ $wgetinstalled == 0 ]]; then
	{ wget https://github.com/MesserLab/SLiM/releases/download/v3.7.1/SLiM.zip && unzip SLiM.zip ;} || { printf "Failed to download SLiM.zip or unzip it.\n"; exit;}
else { exit;} #Exit if neither curl nor wget is installed.
fi

# Proceed with building and installing if all tests succeeded.
{ mkdir BUILD && cd BUILD ;} || { printf "Root is unable to create /tmp/BUILD. It likely already exists. Try again after deleting it."; exit ;}
{ cmake ../SLiM && make -j"$(nproc)" ;} || { printf "Build failed. Please see the output and make a post on the slim-discuss mailing list.\nThe output from this build is stored in '/var/log/' as SLiM-CMakeOutput-%s.log.\nYou may be asked to upload this file during a support request." "$(date -Is)"; mv /tmp/BUILD/CMakeFiles/CMakeOutput.log /var/log/SLiM-CMakeOutput-"$(date -Is)".log; exit;}

{ mkdir -p /usr/bin /usr/share/applications /usr/share/metainfo/;} || { echo "Some directory necessary for installation was not successfully created. Please see the output and make a post on the slim-discuss mailing list."; exit;}

# Exit if installation unsuccessful.
install slim eidos /usr/bin || { echo "Installation to /usr/bin was unsuccessful. Please see the output and make a post on the slim-discuss mailing list."; exit;}
echo "Installation to /usr/bin was successful. Proceeding with desktop integration.";
	
cd ~ || printf "For some reason could not change to ~ before deleting temporary directories."; rm -Rf /tmp/SLiM/ /tmp/BUILD/ /tmp/SLiM.zip || echo "Could not remove temporary files."
