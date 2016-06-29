# $Header$

# put the path to where the Python source was untar'd here (not where it
# is installed)
PYSRCDIR=	/Library/Frameworks/Python.framework/Versions/2.4
MKHOWTO=	PATH=$$PATH:$(HOME)/Desktop/Python/PythonDocTools/tools/ mkhowto
#PAPER=		a4
PAPER=		letter

M_CHAPTERS=	c01-intro.tex c02-samples.tex \
		c03-except.tex c04-utils.tex c05-parse.tex c06-tc.tex \
		c07-writer.tex c08-fault.tex c09-resolve.tex \
		c10-dispatch.tex 
	
G_CHAPTERS=	guide01-intro.tex guide02-wsdl2py.tex

M=		zsi
G=		guide

SOURCES=	$M.tex version.tex $(M_CHAPTERS) 
GSOURCES=	$G.tex version.tex $(G_CHAPTERS)


all:	ps pdf html
pdf:	$M.pdf $G.pdf
ps:	$M.ps $G.ps
html:	$M.html $G.html
.PHONY: all pdf ps html

$M.ps: $(SOURCES)
	$(MKHOWTO) --dvi $M.tex
	dvips -N0 -t $(PAPER) -o $@ $M.dvi
	rm $M.dvi

$G.ps: $(GSOURCES)
	$(MKHOWTO) --dvi $G.tex
	dvips -N0 -t $(PAPER) -o $@ $G.dvi
	rm $G.dvi

$M.pdf: $(SOURCES)
	$(MKHOWTO) --pdf $M.tex

$G.pdf: $(GSOURCES)
	$(MKHOWTO) --pdf $G.tex

$M.html: $(SOURCES)
	$(MKHOWTO) --html --split=1 $M.tex
	mv $M/$M.html $M/$M.css .
	rm -rf $M zsi.l2h

$G.html: $(GSOURCES)
	$(MKHOWTO) --html --split=1 $G.tex
	mv $G/$G.html $G/$G.css .
	rm -rf $G zsi.l2h

#zsixsd.tex: zsi.xsd
#	sed -e '1,/-->/d' zsi.xsd | expand >$@

clean:
	rm -f *.aux *.idx *.ilg *.ind *.log *.toc *.bkm *.syn *.pla *.how
clobber: clean
	rm -f $M.ps $M.pdf $M.html $M.l2h $M.css $G.ps $G.pdf $G.html zsixsd.tex 
realclean: clobber
distclean: clobber
.PHONY: clean clobber realclean distclean
