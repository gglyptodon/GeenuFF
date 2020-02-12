import matplotlib as plt
from matplotlib.pyplot import subplots
import ipywidgets as widgets
from ipywidgets import Layout

from geenuff.base import types

from .wipDrawableClasses import DrawableSuperLocus


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
