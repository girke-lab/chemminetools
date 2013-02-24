"""
This file is copied from ChemMineV2 with modification and some addition. It 
deals with the main logic in rendering the tree, guessing compounds from it,
and serving tiles for Google Map API.
"""
import os
import Image
from cStringIO import StringIO
from tempfile import TemporaryFile
import sys
from pdfbasic import process_pdf, increase_width, process_plot, pdf_size
from math import ceil
from simplejson import dumps as jsonize
import md5
from string import Template
from tempfile import NamedTemporaryFile as NTF
from tempfile import mkstemp
from logging import (info, warning, error, debug, critical, root, NOTSET,
	WARNING, INFO)
from django.core.cache import cache
max_resolution = 6

WORK_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'working')

# pregenerate a white image
img = Image.new("RGBA", (256, 256), (256, 256, 256, 256))
storage = StringIO()
img.save(storage, format="png")
white = storage.getvalue()

root.setLevel(WARNING)

# gs command
gs = 'gs -dUseCIEColor -dFirstPage=1 -dLastPage=1 -dUseCropBox -dSAFER'\
		' -dBATCH -dNOPAUSE -r%dx%d -dNOPLATFONTS -sDEVICE=png16m'\
		' -dBackgroundColor=16#ffffff -dTextAlphaBits=4'\
		' -dGraphicsAlphaBits=4 -sOutputFile=%s %s > /dev/null'

R_TMPL = Template("""

check_tree <- function(x, desired, level)
{
  if (attr(x, 'members') == 1) {
    # leaf
    if (attr(x, 'label') == desired) return(list(1, NULL))
    else return(list(0, NULL));
  } else {
    # internal node
    for (i in x) {
      s <- check_tree(i, desired, level)
      if (! is.null(s[[2]])) return(s)        # already find it
      if (s[[1]] > 0) {
        # the leaf is found. check my height
        if (s[[1]] == level) return(list(level, x, attributes(x)))    # got it
        else return (list(s[[1]]+1, NULL))                    # ask parent to check
      }
    }
    # nobody find it. so I can't find it
    return(list(0, NULL))
  }
}

get_leaves <- function(tree)
{
	if (length(attr(tree, 'leaf')) && attr(tree, 'leaf') == TRUE)
		return(attr(tree, 'label'))
	r <- c()
	for (i in tree) {
		r <-c(r, get_leaves(i))
	}
	return(r)
}

treefile <- "${working_dir}/${ref}.dnd"
if (! file.exists(treefile)) {
	library(ape)
	.treefile <- "${working_dir}/${ref}.rawdnd"
	d <- as.dendrogram(as.hclust(read.tree(.treefile)))
	save(d, file=treefile)
} else {
	load(treefile)
}
outfile <- '${output}'


fixtree <- function(x, w=lw)
{
	#attr(x, 'edgePar') <- list(lwd=w*5);
	if (length(attr(x, 'leaf')) == 0) {
		for (i in 1:length(x))
			x[[i]] <- fixtree(x[[i]], attr(x, 'height'));
	} else {
		if (! is.null(attr(x, 'label')))
			attr(x, 'label') <- gsub('_____', ':', attr(x, 'label'))
	}
	x
}

fixsubtree <- function(x, order=1)
# just get orders right, after you cut the subtree out. heatmap.2 depends
# on this process
{
	if (length(attr(x, 'leaf')) == 0) {
		for (i in 1:length(x)) {
			r <- fixsubtree(x[[i]], order)
			x[[i]] <- r[[1]]
			order <- r[[2]]
		}
		return(list(x, order))
	} else {
		a <- attributes(x)
		x <- order
		attributes(x) <- a
		return(list(x, order + 1))
	}
}

d <- fixtree(d)
# extract subtree?
if ("${leaf}" != "") {
	dd <- check_tree(d, '${leaf}', ${level})
	if (! is.null(dd[[2]])) {
		d <- dd[[2]]
		attributes(d) <- dd[[3]]
		class(d) <- 'dendrogram'
		d <- fixsubtree(d)[[1]]
		save(d, file="${working_dir}/${new_ref}.dnd")
	}
}

# size estimation
ncompounds <- attr(d, 'members')
height <- max(ceiling((32 * ncompounds) / 100), 6)
width <- max(ceiling(height / 8), 5)
lw <- height / 64

# try load the datafile, if any
y <- NULL
if ("${datafile}" != "" && file.exists("${working_dir}/${datafile}.userdata"))
{
	y <- read.delim("${working_dir}/${datafile}.userdata", row.names=1)
	y <- as.matrix(y)
	# make sure only relevant data remains
	leaves <- get_leaves(d)
	ncols <- dim(y)[2]
	if (ncols == 1) {
		y <- as.matrix(cbind(y, y))
		ncols <- 2
	}
	y.rownames <- rownames(y)
	y <- matrix(y[y.rownames %in% leaves, ], ncol=ncols)
	missing_leaves <- leaves[! leaves %in% y.rownames]
	if (length(missing_leaves)) {
		nas <- rep(NA, ncols * length(missing_leaves))
		nas <- matrix(nas, nrow=length(missing_leaves))
		rownames(nas) <- missing_leaves
		y <- rbind(y, nas)
	}
	# now make sure the ordering is right
	y <- matrix(y[sort(y.rownames, index.return=T)$$ix,], ncol=ncols)
	ordering <- sort(sort(leaves, index.return=T)$$ix, index.return=T)$$ix
	y <- matrix(y[ordering,], ncol=ncols)
	rownames(y) <- y.rownames
}
pdf(outfile, height=height, width=width, bg='white')
if (is.null(y)) {
	plot(d, horiz=T)
} else {
	library(gplots)
	heatmap.2(as.matrix(y), d, trace='none')
}
dev.off()
q(save='no')
""")

R_TMPL_SCATTERPLOT = Template("""
y.col <- ${col}
x.ref <- '${ref}'
y.ref <- '${datafile}'
workdir <- '${working_dir}/'
outfile <- '${output}'
color.size <- 1000
coordfile <- '${coordfile}'

# get X
x.file <- paste(workdir, x.ref, '.distmat', sep="")
if (! file.exists(x.file)) {
	distmat <- read.delim(paste(workdir, x.ref, ".matrix", sep=""), row.name=1)
	x <- cmdscale(distmat)
	rownames(x) <- rownames(distmat)
	rownames(x) <- gsub('_____', ':', rownames(x))
	save(x, file=x.file)
} else {
	load(x.file)
}

# get Y
colors <- NULL
if (y.ref != "") {
	y.file <- paste(workdir, y.ref, '.ydata', sep="")
	if (! file.exists(y.file)) {
		y <- read.delim(paste(workdir, y.ref, '.userdata', sep=""), row.names=1)
		y <- as.matrix(y)

		# match x and y
		# note that when subslicing y, it can lose dimensions and rownames
		ncols <- dim(y)[2]
		y.rownames <- rownames(y)
		y <- matrix(y[y.rownames %in% rownames(x), ], ncol=ncols)
		missing <-rownames(x)[! rownames(x) %in% y.rownames]
		if (length(missing)) {
			nas <- rep(NA, ncols * length(missing))
			nas <- matrix(nas, ncol=ncols)
			rownames(nas) <- missing
			y <- rbind(y, nas)
		}
		# now make sure the ordering is right
		y <- matrix(y[sort(y.rownames, index.return=T)$$ix,], ncol=ncols)
		ordering <- sort(sort(rownames(x), index.return=T)$$ix, index.return=T)$$ix
		y <- matrix(y[ordering,], ncol=ncols)
		
		save(y, file=y.file)
	} else {
		load(y.file)
	}

# y to colors
	yy <- y[,y.col]
	yy.range <- range(yy[!is.na(yy)])
	yy <- floor((yy - yy.range[1])/(yy.range[2] - yy.range[1]) * 999 + 1)
	colors <- rainbow(color.size)[yy]
	colors[is.na(colors)] <- 'black'
}

write.csv(x, file=coordfile)
pdf.size <- ceiling(sqrt(dim(x)[1]))
pdf(outfile, width=pdf.size, height=pdf.size, bg="white")
plot(range(x[,1]), range(x[,2]), type='n', xlab='x', ylab='y')
if (is.null(colors)) {
	for (i in 1:dim(x)[1]) points(x[i,1], x[i,2])
} else {
	for (i in 1:dim(x)[1]) points(x[i,1], x[i,2], col=colors[i])
}
dev.off()
q(save='no')
""")

def os_run(cmd, msg=None):
	info("invoking: " + cmd)
	if os.system(cmd) != 0:
		if not msg: msg = 'cannot run ' + cmd
		critical(msg)
		raise SupportError

class HTTPBadRequest(Exception): pass
class HTTPFound(Exception): pass
class HTTPNotFound(Exception): pass
class SupportError(Exception): pass

def index(session, ref, mode="", subl=None, subh=None, dref=None,
	col=1, compounds=None, *args, **kargs):
	"""
	prepare and create the google map page for viewing dendrogram image.
	if subl and subh are provided, only view a subtree rooted at subl's
	subh-ancestor. This is done by creating the pdf under a new ref and
	redirect.
	if dref is provided, the data file pointed will be used
	"""
	modes = {"":R_TMPL, "sp":R_TMPL_SCATTERPLOT}
	if mode not in modes: 
		raise HTTPBadRequest("Bad Mode")

	try:
		int(ref, 16)
		if dref: int(dref, 16)
		col = int(col)
	except:
		raise HTTPBadRequest("Bad reference or column ID")

	# build the PDF name and update ref if necessary
	pdffile_prefix = os.path.join(WORK_DIR, mode)
	if dref:
		pdffile_prefix += dref
		pdffile_prefix += str(col)
	pdffile = pdffile_prefix + ref + '.pdf'
	if subl and subh:
		new_ref = md5.md5("%s-%s-%s" % (ref, subl, subh)).hexdigest()
		pdffile = pdffile_prefix + new_ref + '.pdf'
	coordfile = pdffile + '.coord'

	# if not present, get the PDF file by calling R script
	if not os.path.exists(pdffile):
		r_params = dict(working_dir=WORK_DIR,
			ref=ref, output=pdffile, level='0', leaf='', new_ref='',
			datafile='', col=col, coordfile=coordfile)
		if dref:
			r_params['datafile'] = dref
		if subl and subh:
			r_params.update(dict(level=subh, leaf=subl, new_ref=new_ref))
		r_prgm = modes[mode].substitute(r_params)

		t, _ = mkstemp()
		n = os.fdopen(t, 'w')
		n.write(r_prgm)
		n.flush()
		cmd = "R CMD BATCH %s %s.Rout" % (_, _)
		os_run(cmd, 'cannot run R script')
		n.close()
		info(n.name + ".Rout")

	# if subtree is requested, redirect to new ref
	if subl and subh:
		if dref:
			raise HTTPFound(request.path + '?ref=' + new_ref + '&dref=' + dref)
		raise HTTPFound(request.path + '?ref=' + new_ref)

	# processing PDF
	info("PDF at %s" % pdffile)
	if mode == '':
		w, h, pdffile = increase_width(pdffile)
	elif mode == 'sp':
		w,h = pdf_size(pdffile)
	# a session for this <ref> can have multiple images, for example, one for
	# tree viewing and one for MDS style. Therefore, use <img_ref> to refer
	# to a specific image. This is used in serving of tiles
	img_ref = md5.md5(pdffile).hexdigest()
	session[img_ref] = pdffile
	if mode == '':
		leaves, tree = process_pdf(pdffile, compounds)
	elif mode == 'sp':
		leaves = process_plot(pdffile, coordfile, compounds)
		tree = []

	# rendering PDF to PNG
	dir = pdffile + '.dir'
	min_res = 2
	marker_zoom = 5
	while h / 72 * min_res < 128:
		min_res <<= 1
		marker_zoom -= 1
	# see whether folder exists for that pdf
	if not os.path.isdir(dir):
		try: os.unlink(dir)
		except: pass
		os.mkdir(dir)
		# convert pdf
		res = min_res
		for i in range(max_resolution + 1):
			os_run(gs % (res, res, os.path.join(dir, '%d.png' % i),
					pdffile))
			res <<= 1
	tile_height = h / 72 * min_res
	tile_width = w / 72 * min_res
	min_resolution = int(ceil(max(w, h) / 72.0 * 2 / 256) - 1)

	# load data file's friendly name
	if dref:
		if 'drefs' not in session: session['drefs'] = dict()
		if dref not in session['drefs']:
			drefname = 'imported:' + file(os.path.join(WORK_DIR,
						dref + '.desc')).read().strip()
			session['drefs'][dref] = drefname
		drefname = session['drefs'][dref]
		maxcol = len(file(os.path.join(
				WORK_DIR, dref + '.userdata')
				).readline().split('\t')) 
	else:
		drefname = None
		maxcol = 0

	drefs = session.get('drefs', {})
	return dict(img_ref=img_ref, ref=ref, pdfsize=(w, h),
		max_resolution=max_resolution, min_resolution=min_resolution,
		tile_height=tile_height, tile_width=tile_width,
		marker_zoom=marker_zoom, mode=mode,
		leaves=jsonize(leaves), tree=jsonize(tree), 
		drefname=drefname, drefs=session.get('drefs', {}), dref=dref, 
		maxcol=maxcol, col=col)

def data(session, ref, data='', name=''):
	"""
	append data file to a pdf.
	This will simply save the data, create a ref to the data file, and
	then redirect to the viewing page with data argument pointing to
	the data file
	"""
	try:
		int(ref, 16)
	except:
		raise HTTPBadRequest('Bad references or col value')

	content = data.read()

	fname = "0" + md5.md5(content).hexdigest()
	if not name:
		name = 'Unnamed-%s' % fname[:3]
	of = file(os.path.join(WORK_DIR, fname + '.userdata'), 'w')
	of.write(content)
	of.close()
	of = file(os.path.join(WORK_DIR, fname + '.desc'), 'w')
	of.write(name)
	of.close()
	if 'drefs' not in session: session['drefs'] = {}
	session['drefs'][fname] = name
	return fname
	
def full(session, img_ref, zoom):
	"""server the full-size PNG"""
	import pdb
	pdb.set_trace()
	pdffile = session.get(img_ref)
	if not pdffile: raise HTTPNotFound
	dir = pdffile + '.dir'
	zoom = int(zoom)
	if zoom > max_resolution:
		raise HTTPNotFound
	f = os.path.join(dir, str(zoom) + '.png')
	return file(f).read()

def make_tile(img, tx, ty):
	startx = tx * 256
	starty = ty * 256
	bx, by = img.getbbox()[2:]
	if startx >= bx or starty >= by:
		return white
	endx = startx + 256
	endy = starty + 256
	img2 = img.crop((startx, starty, min(endx, bx), min(endy, by)))
	if endx > bx or endy > by:
		img3 = Image.new("RGBA", (256, 256), (256, 256, 256, 256))
		img3.paste(img2, (0, 0, min(endx, bx) - startx, min(endy, by) - starty))
		img2 = img3
	storage = StringIO()
	img2.save(storage, format="png")
	return storage.getvalue()

# serve png tiles
def png(session, img_ref, zoom, tx, ty):
	pdffile = session.get(img_ref)
	if not pdffile: raise HTTPNotFound
	dir = pdffile + '.dir'
	f = os.path.join(dir, str(zoom) + '.png')
	try:
		img = Image.open(f)
		img.load()
		tx = int(tx)
		ty = int(ty)
		zoom = int(zoom)
		assert zoom <= max_resolution
	except:
		raise HTTPNotFound
	pixels = make_tile(img, tx, ty)
	return pixels

# precutting images to tile
# this is new addition to the old ChemMineV2 source, mainly to deal with the
# fact that loaded Image object cannot be cached and therefore serving each
# tile would require the whole PNG be reloaded and parsed

