* Display # of isoforms and ref. Transcripts in plot title.

* Clicking enter in gene name loads gene

* It would be really nice to have autocomplete for gene names, but this is broken in bokeh serve right now: https://github.com/bokeh/bokeh/issues/4870 A search box + auto-populated select for genes is doable but not great because
callbacks on a text (search) box are only done on losing focus or hitting enter; what is needed is a callback
on keypresses.

* Move the application into [directory format](http://bokeh.pydata.org/en/latest/docs/user_guide/server.html#userguide-server-applications-directory) and use callbacks to load annotation and pickle file.
