#
# ${R_HOME}/doc/html/Makefile


srcdir = .
top_srcdir = ../..

top_builddir = ../..
subdir = doc/html

include $(top_builddir)/Makeconf

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
INSTFILES = \
	R.css \
	left.jpg logo.jpg logosm.jpg right.jpg up.jpg \
	about.html \
	packages-head.html \
	template.html \
	thanks.html
NON_SVN_INSTFILES = resources.html
DISTFILES = Makefile.in $(INSTFILES) index-default.html
SUBDIRS = search
OBJECTS = index.html
## <FIXME>
## build-help really should not create packages.html in *this* directory.
CLEANFILES = packages.html packages-foot.html $(OBJECTS)
## </FIXME>
DISTCLEANFILES = Makefile $(NON_SVN_INSTFILES)

all: Makefile R

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

R: Makefile $(OBJECTS) svnonly
	@if test "$(BUILDDIR_IS_SRCDIR)" = no ; then \
	  $(MAKE) rhome=$(top_builddir) install-sources; \
	fi
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@) || exit 1; \
	done

svnonly:
	@if test ! -f "$(srcdir)/resources.html" || test -f non-tarball ; then \
	  (cd ../manual && $(MAKE) html-non-svn) ; \
	  touch non-tarball ; \
	else \
	  if test "$(BUILDDIR_IS_SRCDIR)" = no ; then \
	    for f in $(NON_SVN_INSTFILES); do \
	      $(INSTALL_DATA) $(srcdir)/$${f} . ; \
	    done \
	  fi \
	fi

index.html: FORCE
	@$(PERL) $(top_srcdir)/tools/linkcheck.pl \
	  $(srcdir)/index-default.html > $@ || \
	  $(INSTALL_DATA) $(srcdir)/index-default.html $@
FORCE:

install: svnonly install-message installdirs install-sources
	@for f in $(OBJECTS); do \
	  $(INSTALL_DATA) $${f} "$(rhome)/$(subdir)"; \
	done
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@) || exit 1; \
	done
	@for f in $(NON_SVN_INSTFILES); do \
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
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@); \
	done
	@for f in $(CLEANFILES) $(INSTFILES) $(NON_SVN_INSTFILES) $(OBJECTS); do \
	  rm -f "$(rhome)/$(subdir)/$${f}"; \
	done
	@rmdir "$(rhome)/$(subdir)" 2>/dev/null \
	  || echo "  subdir $(subdir) not removed"

mostlyclean:
clean:
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@); \
	done
	@if test "$(BUILDDIR_IS_SRCDIR)" = no ; then \
	  rm -f $(INSTFILES); \
	fi
	-@test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)
distclean: clean
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@); \
	done
	@if test -f non-tarball ; then \
	  rm -f $(NON_SVN_INSTFILES) non-tarball ; \
	fi
	-@test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)
maintainer-clean: distclean

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done
	@for f in $(NON_SVN_INSTFILES); do \
	  cp -p $${f} $(distdir)/$${f}; \
	done
	@for d in $(SUBDIRS); do \
	  test -d $(distdir)/$${d} \
	    || mkdir $(distdir)/$${d} \
	    || exit 1; \
	  chmod 755 $(distdir)/$${d}; \
	  (cd $${d} && $(MAKE) distdir) \
	    || exit 1; \
	done
