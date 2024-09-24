# import gsMap
project = 'gsMap portal'
copyright = '2024, Liyang, Wenhao'
author = 'Liyang, Wenhao'
# release = gsMap.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
    'nbsphinx',
    'myst_parser',
    'sphinx_charts.charts',
    "sphinxcontrib.jquery",
    "sphinx_inline_tabs",
]

exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# html_theme = 'classic'
# html_theme = 'sphinx_rtd_theme'
# html_theme = "pydata_sphinx_theme"
html_theme = "furo"
html_static_path = ['_static']
templates_path = ['_templates']

html_theme_options = {
    # "light_css_variables": {
    #     "color-brand-primary": "#7C4DFF",
    #     "color-brand-content": "#7C4DFF",
    #     "color-code-background": "#f5f5f5",
    # },
}

# add plotly.js to the build
html_js_files = [
    'https://cdn.plot.ly/plotly-latest.min.js',
]

rst_epilog = "\n.. include:: .special.rst\n"

