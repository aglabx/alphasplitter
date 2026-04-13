.PHONY: all build release debug check test clean fmt install help

BINDIR := target/release
INPUT  ?= input.10kb.fasta
THREADS ?= 8

all: release

build: release

release:
	cargo build --release

debug:
	cargo build

check:
	cargo check

test:
	cargo test --release

fmt:
	cargo fmt

clean:
	cargo clean

install: release
	cargo install --path . --force

# End-to-end pipeline: discover → cut → annotate.
# Usage: make run INPUT=path/to/arrays.10kb.fasta THREADS=96
run: release
	$(BINDIR)/discover_chains $(INPUT) -o chains.json -t $(THREADS)
	$(BINDIR)/motif_cut $(INPUT) -m chains.json -o monomers.tsv -t $(THREADS)
	$(BINDIR)/annotate_cenpb monomers.tsv annotated.tsv

help:
	@echo "AlphaSplitter — make targets"
	@echo "  make             build release binaries (default)"
	@echo "  make debug       build debug binaries"
	@echo "  make check       cargo check"
	@echo "  make test        run tests"
	@echo "  make fmt         cargo fmt"
	@echo "  make clean       cargo clean"
	@echo "  make install     cargo install --path ."
	@echo "  make run INPUT=arrays.10kb.fasta THREADS=96"
	@echo ""
	@echo "Binaries land in $(BINDIR)/ — e.g. $(BINDIR)/discover_chains --help"
