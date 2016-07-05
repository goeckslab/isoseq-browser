* Bokeh has a webGL option that can make zooming and scrolling much faster.  All you need to do is add webGL=True inside the plot function.  Might be worth using with a Try statement to make things more speedy.
http://bokeh.pydata.org/en/0.10.0/docs/user_guide/webgl.html
Don't know if it is still in version 0.11.1 though.

* The other is allowing the user to export the parameters for each gene.  That way the user can recreate a figure without having to write down all the settings manually.  Maybe just dump all the variables to a json file?  It will be helpful for users when they browse to multiple genes and add their own customizations.  Or if they need to recreate the figure for a paper.
