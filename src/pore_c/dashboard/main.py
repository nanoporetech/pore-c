from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Slider
from bokeh.plotting import figure, curdoc
from bokeh.server.server import Server
from bokeh.themes import Theme

import numpy as np
import holoviews as hv
import panel as pn
import param
import hvplot.pandas
import hvplot.dask

pn.extension()
#import holoviews.plotting.bokeh

import sys
from intake import open_catalog

catalog = sys.argv[1]

cat = open_catalog(catalog)
align_df = cat.align_table.to_dask()
read_df = cat.read_table.read()


class ReadLengthHistogram(param.Parameterized):
    min_bin = param.Integer(0, bounds=(0, 10000), step=10)
    max_bin = param.Integer(100000, bounds=(10, 1000000), step=1000)
    num_bins = param.Integer(100, bounds=(10, 1000), step=10)

    def plot(self):
        return read_df.hvplot.hist(y='read_length', bins=self.num_bins, bin_range=(self.min_bin, self.max_bin))

class AlignmentHistogram(param.Parameterized):
    max_bin = param.Integer(30, bounds=(10, 50), step=1)
    #max_bin = param.Integer(100000, bounds=(10, 1000000), step=1000)
    #num_bins = param.Integer(100, bounds=(10, 1000), step=10)

    def plot(self):
        df = (
            read_df.num_pass_aligns.clip(upper=self.max_bin-1).value_counts().to_frame().sort_index()
        )
        return df.hvplot.bar()
        #return read_df.hvplot.hist(y='read_length', bins=self.num_bins, bin_range=(self.min_bin, self.max_bin))



read_length_hist = ReadLengthHistogram()
alignment_hist = AlignmentHistogram()
pane = pn.Column(
    "# Alignments",
    pn.Row(read_length_hist.param, read_length_hist.plot),
    pn.Row(alignment_hist.param, alignment_hist.plot)
)
pane.servable()
