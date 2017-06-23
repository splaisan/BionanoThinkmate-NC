


all:
	rm -f *.pyc

VERSION := 3.0
SVNREV := $(shell svn info | grep Revision | sed 's/Revision: //')

INSTALLROOT := $(HOME)/tools/Pipeline
INSTALLDIR := $(INSTALLROOT)/$(VERSION).$(SVNREV)
INSTALLFILES := $(wildcard *.py)
INSTALLFILES += $(wildcard *.sh)
INSTALLFILES += $(wildcard *.xml)
INSTALLFILES += $(wildcard *.dox)
INSTALLFILES += Analysis

install:
	mkdir -p $(INSTALLDIR)
	cp -rp $(INSTALLFILES) $(INSTALLDIR)/.

prod: install
	cd $(INSTALLROOT) ; /bin/rm -f $(VERSION) ; ln -s $(VERSION).$(SVNREV) $(VERSION)

