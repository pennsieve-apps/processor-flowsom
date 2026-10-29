.PHONY: help image test fixture fixture-demo run clean

SERVICE_NAME ?= processor-flowsom
IMAGE ?= processor-flowsom:dev
PYTHON ?= python3
DEMO_DIR ?= ../workflow-immune-cytometry/data/zenodo-15388817

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test          - R tests on a tiny synthetic fixture (local Rscript with FlowSOM, else inside the image)"
	@echo "make fixture       - data/fixture/demo20k/{fcs,metadata}: 20,000 events per demo file, arcsinh-transformed (needs flowio, numpy)"
	@echo "make run           - build the image and run it in Docker over data/fixture/demo20k/, writing to data/output/"
	@echo "make clean         - remove data/output/"

image:
	docker build --platform linux/amd64 -t $(IMAGE) .

test:
	@if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status = !requireNamespace("FlowSOM", quietly = TRUE))' >/dev/null 2>&1; then \
		Rscript tests/test_main.R; \
	else \
		$(MAKE) image && docker run --rm --platform linux/amd64 --memory 2g --entrypoint Rscript -v "$(CURDIR):/src:ro" $(IMAGE) /src/tests/test_main.R; \
	fi

fixture:
	$(PYTHON) -m tests.make_fixture $(DEMO_DIR) data/fixture/demo20k --events 20000

run: clean
	@test -d data/fixture/demo20k || $(MAKE) fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
