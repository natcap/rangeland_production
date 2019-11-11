DATA_DIR := data


ENV = env
# assuming eq ($(OS),Windows_NT)
NULL := $$null
PROGRAM_CHECK_SCRIPT := .\scripts\check_required_programs.bat
ENV_SCRIPTS = $(ENV)\Scripts
ENV_ACTIVATE = $(ENV_SCRIPTS)\activate
CP := powershell.exe Copy-Item
COPYDIR := $(CP) -Recurse
MKDIR := powershell.exe mkdir -Force -Path
RM := powershell.exe Remove-Item -Force -Recurse -Path
RMDIR := cmd /C "rmdir /S /Q"
# Windows doesn't install a python2 binary, just python.
PYTHON = python
# Just use what's on the PATH for make.  Avoids issues with escaping spaces in path.
MAKE := make
SHELL := powershell.exe
BASHLIKE_SHELL_COMMAND := cmd.exe /C
.DEFAULT_GOAL := windows_installer
JENKINS_BUILD_SCRIPT := .\scripts\jenkins-build.bat
RM_DATA_DIR := $(RM) $(DATA_DIR)
/ := '\'

REQUIRED_PROGRAMS := make zip $(PYTHON) hg makensis

PIP = $(PYTHON) -m pip
VERSION := $(shell $(PYTHON) setup.py --version)
PYTHON_ARCH := $(shell $(PYTHON) -c "import sys; print('x86' if sys.maxsize <= 2**32 else 'x64')")

# Output directory names
DIST_DIR := dist
BUILD_DIR := build

# The fork name and user here are derived from the mercurial path.
# They will need to be set manually (e.g. make FORKNAME=natcap/invest)
# if someone wants to build from source outside of mercurial (like if
# they grabbed a zipfile of the source code)
# FORKUSER should not need to be set from the CLI.
FORKNAME := $(filter-out ssh: http: https:, $(subst /, ,$(shell hg config paths.default)))
FORKUSER := $(word 2, $(subst /, ,$(FORKNAME)))
ifeq ($(FORKUSER),natcap)
	BUCKET := gs://releases.naturalcapitalproject.org
	DIST_URL_BASE := $(BUCKET)/invest/$(VERSION)
else
	BUCKET := gs://natcap-dev-build-artifacts
	DIST_URL_BASE := $(BUCKET)/invest/$(FORKUSER)/$(VERSION)
endif
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


# Python environment management
env:
	$(PYTHON) -m virtualenv --system-site-packages $(ENV)
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install -r requirements.txt -r requirements-gui.txt"
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install -I -r requirements-dev.txt"
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(MAKE) install"

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
	$(BASHLIKE_SHELL_COMMAND) "$(PYTHON) -m pip freeze --all > $(BINARIES_DIR)/package_versions.txt"


# Installers for each platform.
# Windows (NSIS) installer is written to dist/rangeland_production_<version>_x86_Setup.exe
ifeq ($(FORKUSER), natcap)
	INSTALLER_NAME_FORKUSER :=
else
	INSTALLER_NAME_FORKUSER := $(FORKUSER)
endif
WINDOWS_INSTALLER_FILE := $(DIST_DIR)/rangeland_production_$(INSTALLER_NAME_FORKUSER)$(VERSION)_$(PYTHON_ARCH)_Setup.exe
windows_installer: $(WINDOWS_INSTALLER_FILE)
$(WINDOWS_INSTALLER_FILE): $(BINARIES_DIR) build/vcredist_x86.exe
	-$(RM) $(WINDOWS_INSTALLER_FILE)
	makensis /DVERSION=$(VERSION) /DBINDIR=$(BINARIES_DIR) /DARCHITECTURE=$(PYTHON_ARCH) /DFORKNAME=$(INSTALLER_NAME_FORKUSER) installer\windows\rangeland_production_installer.nsi

build/vcredist_x86.exe: | build
	powershell.exe -Command "Start-BitsTransfer -Source https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x86.exe -Destination build\vcredist_x86.exe"

jenkins:
	$(JENKINS_BUILD_SCRIPT)

jenkins_test_ui: env
	$(MAKE) PYTHON=$(ENV_SCRIPTS)/python test_ui

jenkins_test: env $(GIT_TEST_DATA_REPO_PATH)
	$(MAKE) PYTHON=$(ENV_SCRIPTS)/python test

CERT_FILE := StanfordUniversity.crt
KEY_FILE := Stanford-natcap-code-signing-2019-03-07.key.pem
signcode:
	gsutil cp gs://stanford_cert/$(CERT_FILE) $(BUILD_DIR)/$(CERT_FILE)
	gsutil cp gs://stanford_cert/$(KEY_FILE) $(BUILD_DIR)/$(KEY_FILE)
	# On some OS (including our build container), osslsigncode fails with Bus error if we overwrite the binary when signing.
	osslsigncode -certs $(BUILD_DIR)/$(CERT_FILE) -key $(BUILD_DIR)/$(KEY_FILE) -pass $(CERT_KEY_PASS) -in $(BIN_TO_SIGN) -out "signed.exe"
	mv "signed.exe" $(BIN_TO_SIGN)
	rm $(BUILD_DIR)/$(CERT_FILE)
	rm $(BUILD_DIR)/$(KEY_FILE)
	@echo "Installer was signed"

deploy:
	gsutil -m rsync -r $(DIST_DIR)/userguide $(DIST_URL_BASE)/userguide
	gsutil -m rsync -r $(DIST_DIR)/data $(DIST_URL_BASE)/data
	gsutil -m rsync $(DIST_DIR) $(DIST_URL_BASE)
	@echo "Binaries (if they were created) can be downloaded from:"
	@echo "  * $(DOWNLOAD_DIR_URL)/$(subst $(DIST_DIR)/,,$(WINDOWS_INSTALLER_FILE))"
    ifeq ($(BUCKET),gs://releases.naturalcapitalproject.org)  # ifeq cannot follow TABs, only spaces
		gsutil cp "$(BUCKET)/fragment_id_redirections.json" "$(BUILD_DIR)/fragment_id_redirections.json"
		$(PYTHON) scripts/update_installer_urls.py "$(BUILD_DIR)/fragment_id_redirections.json" $(BUCKET) $(notdir $(WINDOWS_INSTALLER_FILE)) $(notdir $(patsubst "%",%,$(MAC_BINARIES_ZIP_FILE)))
		gsutil cp "$(BUILD_DIR)/fragment_id_redirections.json" "$(BUCKET)/fragment_id_redirections.json"
    endif


# Notes on Makefile development
#
# * Use the -drR to show the decision tree (and none of the implicit rules)
#   if a task is (or is not) executing when expected.
# * Use -n to print the actions to be executed instead of actually executing them.
