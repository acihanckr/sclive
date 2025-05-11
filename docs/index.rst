ScLive Documentation
====================

ScLive is a Python package for plotting single cell sequencing data. Its focus is interactivity and compatibility with web applications.
ScLive provides interactive and highly customizable plots commonly used in single cell analysis.
It also has live_dash module to create single cell analysis dashboards from anndata objects.
The resulting dashboard allows to customize various aspects of the plot using an user interface.

ScLive can be installed using pip:

.. code-block:: bash 

   pip install sclive


Once installed and annotated data is loaded, the following code is enough to create a single cell analysis dashboard:

.. code-block:: python

   from sclive.live_dash import create_dash_app, ScLiveDash
   app = create_dash_app(ScLiveDash(adata))

This code will create a Shiny app which can be deployed on a server supporting Shiny apps or using `Shiny server <https://posit.co/products/open-source/shiny-server/>`_.
It can also be run locally using VS Code `Shiny extension <https://marketplace.visualstudio.com/items?itemName=posit.shiny>`_.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   installation
   examples/index
   api/index
