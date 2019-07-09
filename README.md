# vgbrowser
Browser for Graph Genomes built with VG.  Provides visualization for variation within a species of plant or animal. Designed to scale up to thousands of specimens and provide useful visualizations.

**Backend:** Toshiyuki Yokoyama
**Algorithms and Design:** Josiah Seaman


# Developer Instructions
**Environment**: [Anaconda 3.7 ](https://www.anaconda.com/distribution/)
`pip install -r requirements.txt`  # gfapy does not have an anaconda package

**IDE:**  Pycharm Professional 2019.1  
* Travis CI
* Jupyter Notebook - prototyping environment, mature code gets moved to .py files for reuse.  Can be matured into a user manual.

#### Branches
**master** - should always without errors.  This is a release ready for consumption.  Only use Pull Requests to update master.  We will do code review on master.  
**development** - staging for master, stable. Travis tests should pass, if they don't work to fix this first.  Only use Pull Requests to update development branch.  
**other branches** - create your own for developing on a specific issue.  For Issue #3  i3_graph_summarization.  You can merge these manually.  
