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
clean: ## Removes Virtual environment
ifeq ($(OS),Windows_NT)
	@if exist "$(VENV_NAME)" ($(RMDIR) $(VENV_NAME))
	@if exist "src\__pycache__" (cd src && $(RMDIR) __pycache__)
else
	@if [ -d "$(VENV_NAME)" ]; then \
		$(RMDIR) $(VENV_NAME); \
	fi
	@if [ -d "src/__pycache__" ]; then \
		cd src && $(RMDIR) __pycache__; \
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
FASTA_FILE = data/output.fa
LIMIT = 10000
K = 5
SCRIPT = src/main.py

# Target to run the algorithm with parameters
run:
	$(PYTHON) $(SCRIPT) $(FASTA_FILE) --limit $(LIMIT) --k $(K) --gaps True


