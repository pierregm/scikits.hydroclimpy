# -*- coding: utf-8 -*-
#
# NumPy documentation build configuration file, created by
# sphinx-quickstart on Wed May  7 09:45:59 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.

import sys, os, fnmatch
from scikits.hydroclimpy import __version__ as hc_version

# Add the plot_directive.py directive to the path
sys.path.extend([os.path.dirname(__file__)])


# Check Sphinx version
import sphinx
if sphinx.__version__ < "0.6":
    raise RuntimeError("Sphinx 0.6 or newer required")



# ----------------------------------------------------------------------------
# General configuration
# ----------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary',
              'sphinx.ext.coverage', 'sphinx.ext.doctest',
              'sphinx.ext.inheritance_diagram', 'sphinx.ext.intersphinx', 
              'sphinx.ext.pngmath', 'sphinx.ext.todo',
              'numpydoc', 
              'plot_directive']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'contents'

# General substitutions.
project = 'HydroClimpy Documentation'
copyright = '2009, Pierre GERARD-MARCHANT'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
#
# The short X.Y version.
version = hc_version

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# List of documents that shouldn't be included in the build.
#unused_docs = []

# List of directories, relative to source directory, that shouldn't be searched
# for source files.
exclude_trees = []

# The reST default role (used for this markup: `text`) to use for all documents.
default_role = "autolink"

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

trim_footnote_reference_space = True

# Show both class-level docstring and __init__ docstring in class
# documentation
autoclass_content = 'class'

# ----------------------------------------------------------------------------
# Options for HTML output
# ----------------------------------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
html_theme = 'tstheme'
# html_style = 'mpl.css'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = ["."]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "%s v%s Reference Guide (DRAFT)" % (project, version)

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "HydroClimpy"

# The name of an image file (within the static path) to place at the top of
# the sidebar.
html_logo = '_templates/scipyshiny_small.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Content template for the index page.
html_index = 'index.html'

# Custom sidebar templates, maps document names to template names.
html_sidebars = {'index': 'indexsidebar.html'}

# Additional templates that should be rendered to pages, maps page names to
# template names.
html_additional_pages = {'index': 'index.html',}

# If false, no module index is generated.
html_use_modindex = False

# If true, the reST sources are included in the HTML build as _sources/<name>.
#html_copy_source = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
html_use_opensearch = 'False'

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = '.xhtml'

# Output file base name for HTML help builder.
htmlhelp_basename = 'HydroClimpyDoc'

# Pngmath should try to align formulas properly
pngmath_use_preview = True
pngmath_latex_preamble = '\usepackage{textcomp}'


# ----------------------------------------------------------------------------
# Options for LaTeX output
# ----------------------------------------------------------------------------

# The paper size ('letter' or 'a4').
latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
_stdauthor = u"Pierre Gérard-Marchant"
latex_documents = [
  ('index', 'HydroClimpy.tex', 'HydroClimpy Reference Guide',
   _stdauthor, 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

latex_preamble = """
    \usepackage{amsmath}
    \usepackage{amsfonts}
    \usepackage{amssymb}
    \usepackage{txfonts}
"""

# Additional stuff for the LaTeX preamble.
latex_elements = {'preamble': r'''
\usepackage{amsmath}

% In the parameters section, place a newline after the Parameters
% header
\usepackage{expdlist}
\let\latexdescription=\description
\def\description{\latexdescription{}{} \breaklabel}

% Make Examples/etc section headers smaller and more compact
\makeatletter
\titleformat{\paragraph}{\normalsize\py@HeaderFamily}%
            {\py@TitleColor}{0em}{\py@TitleColor}{\py@NormalColor}
\titlespacing*{\paragraph}{0pt}{1ex}{0pt}
\makeatother

% Do not quote function signatures: asd() instead of `asd()`
\renewcommand\samp[1]{\code{#1}}

% Large function headers
\let\oldmethodline=\methodline
\renewcommand\methodline[2]{\oldmethodline{\Large #1}{#2}}

% Large method headers
\let\oldfuncline=\funcline
\renewcommand\funcline[2]{\oldfuncline{\Large #1}{#2}}

% Fix footer/header
\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter.\ #1}}{}}
\renewcommand{\sectionmark}[1]{\markright{\MakeUppercase{\thesection.\ #1}}}


\usepackage{textcomp}
'''
}

# Documents to append as an appendix to all manuals.
latex_appendices = []

# If false, no module index is generated.
latex_use_modindex = False


# ----------------------------------------------------------------------------
# Plot directive configuration
# ----------------------------------------------------------------------------
plot_formats = [('png', 80), ('hires.png', 200), ('pdf', 50)]


# ----------------------------------------------------------------------------
# Intersphinx configuration
# ----------------------------------------------------------------------------
intersphinx_mapping = {'http://docs.python.org/dev': None,
                       'http://docs.scipy.org/doc/numpy': None,
                       'http://docs.scipy.org/doc/scipy/reference': None,
                       'http://matplotlib.sourceforge.net': None,
                       'http://matplotlib.sourceforge.net/basemap/doc/html': None,
                       'http://pytseries.sourceforge.net': None
                       }


# ----------------------------------------------------------------------------
# Numpy extensions
# ----------------------------------------------------------------------------

# If we want to do a phantom import from an XML file for all autodocs
#phantom_import_file = 'dump.xml'

# Edit links
#numpydoc_edit_link = '`Edit </pydocweb/doc/%(full_name)s/>`__'
numpydoc_show_class_members = False

# ----------------------------------------------------------------------------
# Coverage checker
# ----------------------------------------------------------------------------
coverage_ignore_modules = r"""
    """.split()
coverage_ignore_functions = r"""
    test($|_) (some|all)true bitwise_not cumproduct pkgload
    generic\.
    """.split()
coverage_ignore_classes = r"""
    """.split()

coverage_c_path = []
coverage_c_regexes = {}
coverage_ignore_c_items = {}

#
#plotting_scripts_directory = "examples/plotting"

# ----------------------------------------------------------------------------
# Autogenerate
# ----------------------------------------------------------------------------
autosummary_generate = []
for (r, d, f) in os.walk("."):
    autosummary_generate.extend(os.path.join(r, _) 
                                for _ in fnmatch.filter(f, '*.rst')) 
