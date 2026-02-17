.PHONY: test pdf html all preview clean

all: test pdf html  ## Run tests, render PDF and HTML

test:  ## Run sympy derivation tests
	uv run pytest tests/ -v

pdf: test  ## Render PDF (runs tests first)
	quarto render --to pdf

html: test  ## Render HTML site (runs tests first)
	quarto render --to html

preview:  ## Live preview (HTML only, no tests)
	quarto preview

clean:  ## Remove build output
	rm -rf _book/ _freeze/ .quarto/

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'
