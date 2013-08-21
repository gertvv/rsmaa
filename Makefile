read_version = $(shell grep 'Version:' $1/DESCRIPTION | sed 's/Version: //')

PKG_NAME := smaa
PKG_VERSION := $(call read_version,$(PKG_NAME))
PACKAGE := $(PKG_NAME)_$(PKG_VERSION).tar.gz

all: $(PACKAGE)

$(PACKAGE): $(PKG_NAME)/src/*.c $(PKG_NAME)/R/*.R $(PKG_NAME)/man/*.Rd $(PKG_NAME)/DESCRIPTION $(PKG_NAME)/NAMESPACE $(PKG_NAME)/inst/extdata/* $(PKG_NAME)/tests/Examples/* $(PKG_NAME)/tests/*.R
	R CMD build $(PKG_NAME)

install: $(PACKAGE)
	R CMD INSTALL $(PACKAGE)

check: $(PACKAGE)
	R CMD check $(PACKAGE)

.PHONY: all install
