#
# ${R_HOME}/doc/html/search/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = doc/html/search

include $(top_builddir)/Makeconf

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
SOURCES_JAVA = \
	IndexEntry.java \
	IndexStream.java \
	IndexTable.java \
	SearchEngine.java \
	Tracer.java \
	Value.java
OBJECTS_JAVA = $(SOURCES_JAVA:.java=.class)
INSTFILES = \
	$(OBJECTS_JAVA) \
	SearchEngine-foot.html \
	SearchEngine-head.html \
	SearchObject.html
DISTFILES = Makefile.in $(SOURCES_JAVA) $(INSTFILES)

OBJECTS = SearchEngine.html index.txt

CLEANFILES = $(OBJECTS)
DISTCLEANFILES = Makefile

.SUFFIXES:
.SUFFIXES: .class .java

.java.class:
	(cd $(srcdir) && $(JAVAC) $(<F))

all: Makefile R

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

R: Makefile $(OBJECTS) # $(OBJECTS_JAVA)
	@if test "$(BUILDDIR_IS_SRCDIR)" = no ; then \
	  $(MAKE) rhome=$(top_builddir) install-sources; \
	fi

SearchEngine.html: \
  $(srcdir)/SearchEngine-head.html \
  $(srcdir)/SearchEngine-foot.html \
  $(top_srcdir)/doc/KEYWORDS.db
	@cat $(srcdir)/SearchEngine-head.html > $@
	@-$(PERL) $(top_srcdir)/tools/keywords2html.pl \
	  $(top_srcdir)/doc/KEYWORDS.db >> $@
	@cat $(srcdir)/SearchEngine-foot.html >> $@
## <NOTE>
## This is created via
##   (cd $(top_builddir)/src/library && make docs)
index.txt:
	@touch $@
## so we need to remove it ...
## </NOTE>

install: install-message installdirs install-sources
	@for f in $(OBJECTS); do \
	  $(INSTALL_DATA) $${f} "$(rhome)/$(subdir)"; \
	done
installdirs:
	@$(MKINSTALLDIRS) "$(rhome)/$(subdir)"
install-message:
	@echo "installing $(subdir) ..."
install-sources:
	@for f in $(INSTFILES); do \
	  $(INSTALL_DATA) $(srcdir)/$${f} "$(rhome)/$(subdir)"; \
	done
install-strip: install

uninstall:
	@echo "uninstalling $(subdir) ..."
	@rm -rf "$(rhome)/$(subdir)"

mostlyclean:
clean:
	@if test "$(BUILDDIR_IS_SRCDIR)" = no ; then \
	  rm -f $(INSTFILES); \
	fi
	-@test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)
distclean: clean
	-@test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)
maintainer-clean: distclean

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done
