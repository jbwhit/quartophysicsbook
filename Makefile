.PHONY: setup test lint format pdf html all preview clean

all: lint test pdf html  ## Lint, run tests, render PDF and HTML

setup:  ## Enable shared git hooks and sync dependencies (run once after cloning)
	git config core.hooksPath .githooks
	uv sync
	@echo "Hooks enabled (core.hooksPath=.githooks) and deps synced."

test:  ## Run sympy derivation tests
	uv run pytest tests/ -v

lint:  ## Check formatting and lint (matches CI)
	uv run ruff format --check .
	uv run ruff check .

format:  ## Auto-format and apply lint fixes
	uv run ruff format .
	uv run ruff check --fix .

pdf: test  ## Render PDF (runs tests first)
	quarto render --to pdf
	cp _book/Think-Like-A-Physicist.pdf .

html: test  ## Render HTML site (runs tests first)
	quarto render --to html

preview:  ## Live preview (HTML only, no tests)
	quarto preview

clean:  ## Remove build output
	rm -rf _book/ _freeze/ .quarto/

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'
