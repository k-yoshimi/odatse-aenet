project = 'odatse-aenet'
copyright = '2024-, The University of Tokyo'
author = 'ODAT-SE development team'

version = '1.0'
release = '1.0.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []

language = 'ja'
master_doc = 'index'

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

latex_documents = [
    (master_doc, 'usersguide_odatse-aenet_ja.tex',
     'odatse-aenet Users Guide', author, 'manual'),
]
