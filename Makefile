.DEFAULT_GOAL := help
SHELL := /bin/sh

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

install:  ## Create .venv and install everything
	uv sync --all-groups

test:  ## Run the suite - no dependencies at all
	uv run pytest -q

real:  ## Run the pipeline on real NCBI genomes (cached in data/)
	uv run python scripts/real_data.py analyse

lint:  ## Lint
	uv run ruff check src tests
	uv run ruff format --check src tests

fmt:  ## Auto-format
	uv run ruff format src tests
	uv run ruff check --fix src tests

.PHONY: help install test demo real realign lint fmt
