[![Build Status](https://travis-ci.org/subwaystation/vgbrowser.svg?branch=master)](https://travis-ci.org/subwaystation/vgbrowser)

# vgbrowser
Browser for Graph Genomes built with VG.  Provides visualization for variation within a species of plant or animal. Designed to scale up to thousands of specimens and provide useful visualizations.

**Backend:** Toshiyuki Yokoyama
**Algorithms and Design:** Josiah Seaman


# Developer Instructions
**Environment**: [Anaconda 3.7 ](https://www.anaconda.com/distribution/)  
Ggfapy etc. does not have an anaconda package, so it's necessary to use pip:  
`pip install -r requirements.txt`  


**IDE:**  Pycharm Professional 2019.1  
* Travis CI - automatically runs tests on master and development branches
* Jupyter Notebook - run from the same Anaconda environment.  Notebooks are useful for prototyping, mature code gets moved to .py files for reuse.  They can be matured into a user manual.

#### Branches
**master** - should always run tests without errors.  This is currently our development branch until our first release.  Only use Pull Requests to update master.  We will do code review on master.  
**development** - Future branch once we have our first release point.  Staging for master, stable. Travis tests should pass, if they don't work to fix this first.  Only use Pull Requests to update development branch.  
**other branches** - create your own for developing on a specific issue.  For Issue #3  i3_graph_summarization.  You can merge these manually.  
