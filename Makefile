VENV_NAME = venv
REQUIREMENTS_FILE = requirements.txt
PYTHON ?= python3
MAC = Darwin
ifeq ($(OS),Windows_NT)
    VENV_PATH = .\$(VENV_NAME)\Scripts
    ACTIVATE_VENV = $(VENV_PATH)\activate
    RMDIR = rmdir /s /q
    PIP = pip
else
    UNAME_S := $(shell uname -s)
    VENV_PATH = $(VENV_NAME)/bin/activate
    ACTIVATE_VENV = . $(VENV_PATH)
    RMDIR = rm -rf
    PIP = pip3
endif

.PHONY: setup
setup: ## Makes virtual environment and installs requirements
	$(PYTHON) -m venv $(VENV_NAME)
ifeq ($(OS),Windows_NT)
	$(ACTIVATE_VENV) && \
	$(PIP) install -r $(REQUIREMENTS_FILE)
else
	$(ACTIVATE_VENV) && \
	$(PIP) install -r $(REQUIREMENTS_FILE)
endif

.PHONY: clean
clean: ## Removes Virtual environment and logos
ifeq ($(OS),Windows_NT)
	@if exist "$(VENV_NAME)" ($(RMDIR) $(VENV_NAME))
	@if exist "src\__pycache__" (cd src && $(RMDIR) __pycache__)
	@if exist "logos\*.png" (del /Q logos\*.png)
	@if exist "logos\*.jpg" (del /Q logos\*.jpg)
	@if exist "logos\*.jpeg" (del /Q logos\*.jpeg)
else
	@if [ -d "$(VENV_NAME)" ]; then \
		$(RMDIR) $(VENV_NAME); \
	fi
	@if [ -d "src/__pycache__" ]; then \
		cd src && $(RMDIR) __pycache__; \
	fi
	@if [ -d "logos" ]; then \
		rm -f logos/*.png logos/*.jpg logos/*.jpeg; \
	fi
endif
.PHONY: help
help: ## List all available make commands
ifeq ($(OS),Windows_NT)
	$(ACTIVATE_VENV) && \
	$(PYTHON) list_make_commands.py && \
	deactivate
else
	$(ACTIVATE_VENV) && \
	$(PYTHON) list_make_commands.py && \
	deactivate
endif

# Default FASTA file and parameters
FASTA_FILE = data/test_data.fa
LIMIT = 100000
K = 4
SCRIPT = src/main.py

# Target to run the algorithm with parameters
run:
	$(ACTIVATE_VENV) && $(PYTHON) $(SCRIPT) $(FASTA_FILE) --limit $(LIMIT) --k $(K)

