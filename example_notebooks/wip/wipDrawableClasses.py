import matplotlib as plt
from matplotlib.pyplot import subplots


from geenuff.base import types
from dna_features_viewer import GraphicFeature, GraphicRecord

from .wipDrawableExceptions import InvalidFileExtensionException
from .wipClasses import Feature, Transcript, SuperLocus, _repr


class DrawableExtensions:
    ALLOWED_EXTENSIONS = ('.png', '.jpg', '.pdf', '.svg')


def color(geenuff_feature_name):
    valid = [t.value for t in types.GeenuffFeature]
    if geenuff_feature_name not in valid:
        raise Exception("{} is not a valid GeenuffFeature name".format(geenuff_feature_name))
    colmap = {
        'geenuff_transcript': '#00bfff',
        'geenuff_cds': '#ba55d3',
        'geenuff_intron': '#ffdab9',
        'missing_utr_5p': '#dc143c',
        'missing_utr_3p': '#dc143c',
        'empty_super_locus': '#dc143c',
        'missing_start_codon': '#dc143c',
        'missing_stop_codon': '#dc143c',
        'wrong_starting_phase': '#dc143c',
        'mismatched_ending_phase': '#dc143c',
        'overlapping_exons': '#dc143c',
        'too_short_intron': '#dc143c'
    }
    return colmap.get(geenuff_feature_name, '#c0c0c0')


def convert_strand_info(geenuff_is_plus_strand):  # todo
    # print(geenuff_is_plus_strand)
    if geenuff_is_plus_strand:
        return +1
    else:
        return -1


class DrawableSuperLocus:
    @classmethod
    def filter_by_given_name(cls, super_locus, coordinate_piece, given_names): # todo: under construction
        #print("called")
        new_transcripts = []
        for t in super_locus.transcripts:
            if t.given_name not in given_names:
                pass
                # print("filtered out {}".format(t.given_name))
            else:
                new_transcripts.append(t)
                # print("appended {}".format(t.given_name))

        return cls(super_locus=SuperLocus(is_fully_contained=super_locus.is_fully_contained,
                                          transcripts=new_transcripts,
                                          type=super_locus.type,
                                          id=super_locus.id,
                                          overlaps=super_locus.overlaps,
                                          given_name=super_locus.given_name),
                   coordinate_piece=coordinate_piece)

    def __init__(self, super_locus, coordinate_piece):
        self.super_locus = super_locus
        self.coordinate_piece = coordinate_piece
        self.graphic_features = []
        self.pos_min = float('inf')
        self.pos_max = float('-inf')
        try:
            for t in self.super_locus.transcripts:
                for feature in t.features:
                    if feature.start < self.pos_min:
                        self.pos_min = feature.start
                    if feature.end > self.pos_max:
                        self.pos_max = feature.end

                    if feature.type in ['geenuff_transcript']:
                        self.graphic_features.append(
                            GraphicFeature(
                                start=feature.start,
                                end=feature.end,
                                strand=convert_strand_info(feature.is_plus_strand),
                                color=color(feature.type), label=t.given_name),

                        )
                    else:
                        self.graphic_features.append(
                            GraphicFeature(
                                start=feature.start,
                                end=feature.end,
                                strand=convert_strand_info(feature.is_plus_strand),
                                color=color(feature.type), label=t.given_name),
                        )
                self.graphic_record = GraphicRecord(sequence=coordinate_piece.sequence, features=self.graphic_features)
        except Exception as e:
            print(e)
            raise Exception  # todo: Exception anticipation and handling

    def draw(self, zoom_coordinates=None, save_to=None, **kwargs):
        """zoom_coordinates expects a tuple (start,end)"""
        if zoom_coordinates is None:
            gr = GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence[self.pos_min:self.pos_max], first_index=self.pos_min)
            ax, _ = gr.plot(**kwargs)
            if save_to:
                print("saving to {}".format(save_to))
                if save_to.endswith(DrawableExtensions.ALLOWED_EXTENSIONS):
                    try:
                        #gr.plot(**kwargs)#figure_width=15)
                        ax.figure.savefig(save_to, bbox_inches='tight')
                    except Exception as e:
                        print(e)
                        raise Exception
                else:
                    raise InvalidFileExtensionException(
                        "Valid file extensions are: {}".format(",".join(DrawableExtensions.ALLOWED_EXTENSIONS)))
        else:
            record = GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence[self.pos_min:self.pos_max], first_index=self.pos_min)
            cropped_record = record.crop(zoom_coordinates)
            #fig, (ax1, ax2) = plt.pyplot.subplots(1, 2, figsize=(14, 2)) #todo: this is weird
            #ax1.set_title("Superlocus " + str(self.super_locus.given_name), loc='left', weight='bold')
            #record.plot(ax=ax1)
            #cropped_record.plot_translation(ax=ax2, location=(400, 400), fontdict={'weight': 'bold'})
            ax, _ = cropped_record.plot(plot_sequence=True, **kwargs)
            #ax2.set_title("Sequence detail " + str(self.super_locus.given_name), loc='left', weight='bold')
            if save_to:
                print("saving to {}".format(save_to))
                if save_to.endswith(DrawableExtensions.ALLOWED_EXTENSIONS):
                    try:
                        #ax, _ = cropped_record.plot(plot_sequence=True, **kwargs)#figure_width=15)
                        ax.figure.savefig(save_to, bbox_inches='tight')
                    except Exception as e:
                        print(e)
                        raise Exception
                else:
                    raise InvalidFileExtensionException(
                        "Valid file extensions are: {}".format(",".join(DrawableExtensions.ALLOWED_EXTENSIONS)))

    def __repr__(self):
        return _repr(self)


class DrawableGeenuffCollection:
    #todo from 'normal\ geenuff collection
    @classmethod
    def from_geenuff_collection(cls, geenuff_collection):
        tmp = []
        for sl in geenuff_collection:
            tmp.append(DrawableSuperLocus(super_locus=sl, coordinate_piece=geenuff_collection.coordinate_piece))
        return cls(list_of_drawable_super_loci=tmp, coordinate_piece=geenuff_collection.coordinate_piece)

    @classmethod
    def from_set(cls, set_of_super_loci, coordinate_piece): # eg from filter, but maybe rather return new collection from filter?
        tmp = []
        for s in set_of_super_loci:
            tmp.append(DrawableSuperLocus(super_locus=s, coordinate_piece=coordinate_piece))
        return cls(list_of_drawable_super_loci=tmp, coordinate_piece=coordinate_piece)

    def __init__(self, list_of_drawable_super_loci, coordinate_piece):
        def _init_graphic_features():
            res = []
            for sl in list_of_drawable_super_loci:
                res += sl.graphic_features
            return res

        self.coordinate_piece = coordinate_piece
        self.drawable_super_loci = list_of_drawable_super_loci
        self.graphic_features = _init_graphic_features()

    def draw(self, zoom_coordinates=None, save_to=None):
        """zoom_coordinates expects a tuple (start,end)"""
        if zoom_coordinates is None:
            GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence).plot(figure_width=10)
        else:
            record = GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence)
            zoom_start, zoom_end = zoom_coordinates
            cropped_record = record.crop((zoom_start, zoom_end))
            fig, (ax1, ax2) = plt.pyplot.subplots(1, 2, figsize=(14, 2)) #todo: this is weird
            ax1.set_title("", loc='left', weight='bold') #+ str(self.super_locus.given_name), loc='left', weight='bold')

            record.plot(ax=ax1)
            cropped_record.plot_translation(ax=ax2, location=(400, 400), fontdict={'weight': 'bold'})
            cropped_record.plot(ax=ax2, plot_sequence=True)
            ax2.set_title("Sequence detail ", loc='left', weight='bold') #str(self.super_locus.given_name),
            if save_to:
                print("saving to {}".format(save_to))
                if save_to.endswith(DrawableExtensions.ALLOWED_EXTENSIONS):
                    try:
                        fig.savefig(save_to, bbox_inches='tight')
                    except Exception as e:
                        print(e)
                        raise Exception
                else:
                    raise InvalidFileExtensionException(
                        "Valid file extensions are: {}".format(",".join(DrawableExtensions.ALLOWED_EXTENSIONS)))

    def __repr__(self):
        _repr(self)

