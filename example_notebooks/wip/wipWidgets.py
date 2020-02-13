import matplotlib as plt
from matplotlib.pyplot import subplots
import ipywidgets as widgets
from ipywidgets import Layout

from geenuff.base import types

from .wipDrawableClasses import DrawableSuperLocus, DrawableGeenuffCollection
from .wipClasses import GeenuffCollectionFilter


class OverviewWidget:

    def __init__(self, collection):
        _opts = [i.given_name for i in collection.super_loci]

        self.selectwidget = widgets.SelectMultiple(
            options=_opts,
            value=(),
            # rows=10,
            description='Superloci:',
            disabled=False
        )

        def _on_button_clicked(_):
            res = GeenuffCollectionFilter.filter(collection, filter_by='given_name', to_match=self.selectwidget.value)
            #print(res)
            #for r in res.super_loci:
                #print(r)
            #    DrawableSuperLocus(
            #        super_locus=r, coordinate_piece=collection.coordinate_piece
            #    ).draw()
            DrawableGeenuffCollection.from_geenuff_collection(res).draw()
            #print(res)
            #res.draw()

        self.button = widgets.Button(description='Show')
        self.button.on_click(_on_button_clicked)
        self.container = widgets.VBox([self.selectwidget, self.button])

    def display(self):
        return self.container



class DummyWidget:

    def __init__(self, drawableLocus):
        valid = [t.value for t in types.GeenuffFeature]
        t = [e.given_name for e in drawableLocus.super_locus.transcripts]
        self.menu = widgets.Dropdown(
            options=t,
            value=t[0],
            description='Transcripts:')

        self.zoom_slider = widgets.IntRangeSlider(
            value=[drawableLocus.pos_min, drawableLocus.pos_max],
            min=drawableLocus.pos_min,
            max=drawableLocus.pos_max,
            step=1,
            description='Zoom Pos:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='d',
            layout=Layout(width='50%', height='80px')
        )
        self.button = widgets.Button(description='Plot')
        self.out = widgets.Output()

        def on_button_clicked(_):
            result = []
            # "linking function with output"
            with self.out:
                # what happens when we press the button
                self.out.clear_output()
               # print(self.menu.value)
                res = DrawableSuperLocus.filter_by_given_name(super_locus=drawableLocus.super_locus,
                                                              coordinate_piece=drawableLocus.coordinate_piece,
                                                              given_names=[self.menu.value])
                #print(res)
                with self.out:
                    #print("sth", res)
                    res.draw(zoom_coordinates=self.zoom_slider.value)
                    #drawableLocus.draw(zoom_coordinates=self.zoom_slider.value)
                    #DrawableSuperLocus(r).draw()

        self.button.on_click(on_button_clicked)

        self.box = widgets.VBox([self.zoom_slider, self.menu, self.button, self.out])

    def display(self):
        return self.box
