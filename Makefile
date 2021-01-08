DATA_DIR := data


ENV = "./env"
# assuming eq ($(OS),Windows_NT)
NULL := $$null
PROGRAM_CHECK_SCRIPT := .\scripts\check_required_programs.bat
ENV_SCRIPTS = $(ENV)\Scripts
ENV_ACTIVATE = $(ENV_SCRIPTS)\activate
CP := cp
COPYDIR := $(CP) -r
MKDIR := mkdir -p
RM := rm
RMDIR := $(RM) -r
# Windows doesn't install a python3 binary, just python.
PYTHON = python
# Just use what's on the PATH for make.  Avoids issues with escaping spaces in path.
MAKE := make
# Powershell has been inconsistent for allowing make commands to be
# ignored on failure. Many times if a command writes to std error
# powershell interprets that as a failure and exits. Bash shells are
# widely available on Windows now, especially through git-bash
SHELL := /usr/bin/bash
CONDA := conda.bat
BASHLIKE_SHELL_COMMAND := $(SHELL) -c
.DEFAULT_GOAL := windows_installer
JENKINS_BUILD_SCRIPT := .\scripts\jenkins-build.bat
RM_DATA_DIR := $(RM) $(DATA_DIR)
/ := '\'

REQUIRED_PROGRAMS := make zip $(PYTHON) conda makensis

PIP = $(PYTHON) -m pip
VERSION := $(shell $(PYTHON) setup.py --version)
PYTHON_ARCH := $(shell $(PYTHON) -c "import sys; print('x86' if sys.maxsize <= 2**32 else 'x64')")

GSUTIL := gsutil
SIGNTOOL := SignTool

# Output directory names
DIST_DIR := dist
BUILD_DIR := build

# The fork name and user here are derived from the git path on github.
# The fork name will need to be set manually (e.g. make FORKNAME=natcap/invest)
# if someone wants to build from source outside of git (like if they grabbed
# a zipfile of the source code).
FORKNAME := $(word 2, $(subst :,,$(subst github.com, ,$(shell git remote get-url origin))))
FORKUSER := $(word 1, $(subst /, ,$(FORKNAME)))

# We use these release buckets here in Makefile and also in our release scripts.
# See scripts/release-3-publish.sh.
RELEASES_BUCKET := gs://releases.naturalcapitalproject.org
DEV_BUILD_BUCKET := gs://natcap-dev-build-artifacts

DOWNLOAD_DIR_URL := $(subst gs://,https://storage.googleapis.com/,$(DIST_URL_BASE))


TESTRUNNER := $(PYTHON) -m nose -vsP --with-coverage --cover-package=rangeland_production --cover-erase --with-xunit --cover-tests --cover-html --cover-xml --logging-level=DEBUG --with-timer

# Target names.
BINARIES_DIR := $(DIST_DIR)/rangeland_production

.PHONY: fetch install binaries windows_installer clean help check python_packages jenkins purge deploy signcode

# Very useful for debugging variables!
# $ make print-FORKNAME, for example, would print the value of the variable $(FORKNAME)
print-%:
	@echo "$* = $($*)"

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  check             to verify all needed programs and packages are installed"
	@echo "  env               to create a virtualenv with packages from requirements.txt, requirements-dev.txt"
	@echo "  fetch             to clone all managed repositories"
	@echo "  install           to build and install a wheel of rangeland_production into the active python installation"
	@echo "  binaries          to build pyinstaller binaries"
	@echo "  python_packages   to build rangeland_production wheel and source distributions"
	@echo "  windows_installer to build an NSIS installer for distribution"
	@echo "  clean             to remove temporary directories and files (but not dist/)"
	@echo "  purge             to remove temporary directories, cloned repositories and the built environment."
	@echo "  help              to print this help and exit"

$(BUILD_DIR) $(DATA_DIR) $(DIST_DIR):
	$(MKDIR) $@

clean:
	$(PYTHON) setup.py clean
	-$(RMDIR) $(BUILD_DIR)
	-$(RMDIR) rangeland_production.egg-info
	-$(RMDIR) cover
	-$(RM) coverage.xml

purge: clean
	-$(RMDIR) $(ENV)

check:
	@echo "Checking required applications"
	@$(PROGRAM_CHECK_SCRIPT) $(REQUIRED_PROGRAMS)
	@echo "----------------------------"
	@echo "Checking python packages"
	@$(PIP) freeze --all -r requirements.txt -r requirements-dev.txt > $(NULL)


# Python conda environment management
env:
	@echo "NOTE: requires 'requests' be installed in base Python"
	$(PYTHON) ./scripts/convert-requirements-to-conda-yml.py requirements.txt requirements-dev.txt requirements-gui.txt > requirements-all.yml
	$(CONDA) create -p $(ENV) -y -c conda-forge python=3.8 nomkl
	$(CONDA) env update -p $(ENV) --file requirements-all.yml
	@echo "----------------------------"
	@echo "To finish the conda env install:"
	@echo ">> conda activate $(ENV)"
	@echo ">> make install"

# compatible with pip>=7.0.0
# REQUIRED: Need to remove natcap.invest.egg-info directory so recent versions
# of pip don't think CWD is a valid package.
install: $(DIST_DIR)/rangeland_production%.whl
	-$(RMDIR) rangeland_production.egg-info
	$(PIP) install --isolated --upgrade --only-binary rangeland_production --find-links=dist rangeland_production


# Bulid python packages and put them in dist/
python_packages: $(DIST_DIR)/rangeland_production%.whl $(DIST_DIR)/rangeland_production%.zip
$(DIST_DIR)/rangeland_production%.whl: | $(DIST_DIR)
	$(PYTHON) setup.py bdist_wheel

$(DIST_DIR)/rangeland_production%.zip: | $(DIST_DIR)
	$(PYTHON) setup.py sdist --formats=zip


# Build binaries and put them in dist/rangeland_production
binaries: $(BINARIES_DIR)
$(BINARIES_DIR): | $(DIST_DIR) $(BUILD_DIR)
	-$(RMDIR) $(BUILD_DIR)/pyi-build
	-$(RMDIR) $(BINARIES_DIR)
	$(PYTHON) -m PyInstaller --workpath $(BUILD_DIR)/pyi-build --clean --distpath $(DIST_DIR) exe/rangeland_production.spec
	$(CONDA) list --export > $(BINARIES_DIR)/package_versions.txt


# Installers for each platform.
# Windows (NSIS) installer is written to dist/rangeland_production_<version>_x86_Setup.exe
INSTALLER_NAME_FORKUSER := $(FORKUSER)
WINDOWS_INSTALLER_FILE := $(DIST_DIR)/rangeland_production_$(VERSION)_$(PYTHON_ARCH)_Setup.exe
windows_installer: $(WINDOWS_INSTALLER_FILE)
$(WINDOWS_INSTALLER_FILE): $(BINARIES_DIR) build/vcredist_x86.exe
	-$(RM) $(WINDOWS_INSTALLER_FILE)
	makensis /DVERSION=$(VERSION) /DBINDIR=$(BINARIES_DIR) /DARCHITECTURE=$(PYTHON_ARCH) /DFORKNAME=$(INSTALLER_NAME_FORKUSER) installer\windows\rangeland_production_installer.nsi

build/vcredist_x86.exe: | build
	powershell.exe -Command "Start-BitsTransfer -Source https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x86.exe -Destination build\vcredist_x86.exe"

CERT_FILE := StanfordUniversity.crt
KEY_FILE := Stanford-natcap-code-signing-2019-03-07.key.pem
P12_FILE := Stanford-natcap-code-signing-2019-03-07.p12
KEYCHAIN_NAME := codesign_keychain
signcode:
	$(GSUTIL) cp gs://stanford_cert/$(CERT_FILE) $(BUILD_DIR)/$(CERT_FILE)
	$(GSUTIL) cp gs://stanford_cert/$(KEY_FILE) $(BUILD_DIR)/$(KEY_FILE)
	# On some OS (including our build container), osslsigncode fails with Bus error if we overwrite the binary when signing.
	osslsigncode -certs $(BUILD_DIR)/$(CERT_FILE) -key $(BUILD_DIR)/$(KEY_FILE) -pass $(CERT_KEY_PASS) -in $(BIN_TO_SIGN) -out "signed.exe"
	mv "signed.exe" $(BIN_TO_SIGN)
	$(RM) $(BUILD_DIR)/$(CERT_FILE)
	$(RM) $(BUILD_DIR)/$(KEY_FILE)
	@echo "Installer was signed with osslsigncode"

signcode_windows:
	$(GSUTIL) cp 'gs://stanford_cert/$(P12_FILE)' '$(BUILD_DIR)/$(P12_FILE)'
	powershell.exe "& '$(SIGNTOOL)' sign /f '$(BUILD_DIR)\$(P12_FILE)' /p '$(CERT_KEY_PASS)' '$(BIN_TO_SIGN)'"
	-$(RM) $(BUILD_DIR)/$(P12_FILE)
	@echo "Installer was signed with signtool"

deploy:
	-$(GSUTIL) -m rsync $(DIST_DIR) $(DIST_URL_BASE)
	-$(GSUTIL) -m rsync -r $(DIST_DIR)/data $(DIST_URL_BASE)/data
	-$(GSUTIL) -m rsync -r $(DIST_DIR)/userguide $(DIST_URL_BASE)/userguide
	@echo "Application binaries (if they were created) can be downloaded from:"
	@echo "  * $(DOWNLOAD_DIR_URL)/$(subst $(DIST_DIR)/,,$(WINDOWS_INSTALLER_FILE))"


# Notes on Makefile development
#
# * Use the -drR to show the decision tree (and none of the implicit rules)
#   if a task is (or is not) executing when expected.
# * Use -n to print the actions to be executed instead of actually executing them.
