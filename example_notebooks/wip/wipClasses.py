import matplotlib as plt
from matplotlib.pyplot import subplots

from geenuff.base import types
from dna_features_viewer import GraphicFeature, GraphicRecord


class NotGeenuffDeserializableException(Exception):
    pass


class DrawableExtensions:
    ALLOWED_EXTENSIONS = ('.png', '.jpg', '.pdf', '.svg')


class InvalidFileExtensionException(Exception):
    pass


class Transcript:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        def _init_features():
            feats = dct["features"]
            res = []
            for f in feats:
                res.append(Feature.from_dct(f))
            return res
        return cls(
            is_fully_contained=dct["is_fully_contained"],
            given_name=dct["given_name"],
            id=dct["id"],
            type=dct["type"],
            overlaps=dct["overlaps"],
            features=_init_features()
        )

    def __init__(self, is_fully_contained, features, type, id, overlaps, given_name):
        self.is_fully_contained = is_fully_contained
        self.features = features
        self.type = type
        self.id = id
        self.overlaps = overlaps
        self.given_name = given_name


class SuperLocus:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        def _init_transcripts():
            transcripts = dct["transcripts"]
            res = []
            for t in transcripts:
                res.append(Transcript.from_dct(t))
            return res
        return cls(
            is_fully_contained=dct["is_fully_contained"],
            given_name=dct["given_name"],
            id=dct["id"],
            type=dct["type"],
            overlaps=dct["overlaps"],
            transcripts=_init_transcripts()
        )

    def __init__(self, is_fully_contained, transcripts, type, id, overlaps, given_name):
        self.is_fully_contained = is_fully_contained
        self.transcripts = transcripts
        self.type = type
        self.id = id
        self.overlaps = overlaps
        self.given_name = given_name


class CoordinatePiece:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        return CoordinatePiece(
            seqid=dct["seqid"],
            id=dct["id"],
            sequence=dct["sequence"],
            start=dct["start"],
            end=dct["end"])

    def __init__(self, seqid, id, sequence, start, end):
        self.seqid = seqid
        self.id = id
        self.sequence = sequence
        self.start = start
        self.end = end

    def __repr__(self):
        try:
            return "seqid:{},id:{},sequence:{}...{},start:{},end:{}".format(self.seqid, self.id, self.sequence[:6], self.sequence[-6:], self.start, self.end)
        except Exception as e:
            print(e)
            return "seqid:{},id:{},sequence:...,start:{},end:{}".format(self.seqid, self.id, self.start, self.end)


class Feature:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        return cls(
            end=dct["end"],
            score=dct["score"],
            is_plus_strand=dct["is_plus_strand"],
            start_is_biological_start=dct["start_is_biological_start"],
            start=dct["start"],
            phase=dct["phase"],
            is_fully_contained=dct["is_fully_contained"],
            given_name=dct["given_name"],
            id=dct["id"],
            protein_id=dct["protein_id"],
            type=dct["type"],
            end_is_biological_end=dct["end_is_biological_end"],
            overlaps=dct["overlaps"],
            source=dct["source"]
        )

    def __init__(self, end, score, is_plus_strand, start_is_biological_start, start, phase, is_fully_contained,
                 given_name, id, protein_id, type, end_is_biological_end, overlaps, source):
        if type not in [t.value for t in types.GeenuffFeature]:
            raise NotGeenuffDeserializableException

        self.end = end
        self.score = score
        self.is_plus_strand = is_plus_strand
        self.start_is_biological_start = start_is_biological_start
        self.start = start
        self.phase = phase
        self.is_fully_contained = is_fully_contained
        self.given_name = given_name
        self.id = id
        self.protein_id = protein_id
        self.type = type
        self.end_is_biological_end = end_is_biological_end
        self.overlaps = overlaps
        self.source = source

    def __repr__(self):
        return "{},{},{}".format(self.start, self.end, self.type)


class GeenuffCollection:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        def _init_super_loci():
            #print(dct)
            super_loci = dct["super_loci"]
            res = []
            for sl in super_loci:
                res.append(SuperLocus.from_dct(sl))
            return res

        return cls(
            super_loci=_init_super_loci(),
            coordinate_piece=CoordinatePiece.from_dct(dct=dct["coordinate_piece"])
        )

    def __init__(self, super_loci, coordinate_piece):
        self.super_loci = super_loci
        self.coordinate_piece = coordinate_piece


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
                        print("updated min:", self.pos_min, feature.start)
                        self.pos_min = feature.start
                    if feature.end > self.pos_max:
                        print("updated max:", self.pos_max, feature.end)
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

    def draw(self, zoom_coordinates=None, save_to=None):
        """zoom_coordinates expects a tuple (start,end)"""
        if zoom_coordinates is None:
            #print(self.coordinate_piece.sequence[self.pos_min:self.pos_max], self.pos_min, self.pos_max)

            #GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence, sequence_length=500).plot(figure_width=10)
            GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence[self.pos_min:self.pos_max], first_index=self.pos_min).plot(figure_width=10)

        else:
            record = GraphicRecord(features=self.graphic_features, sequence=self.coordinate_piece.sequence)
            zoom_start, zoom_end = zoom_coordinates
            cropped_record = record.crop((zoom_start, zoom_end))
            fig, (ax1, ax2) = plt.pyplot.subplots(1, 2, figsize=(14, 2)) #todo: this is weird
            ax1.set_title("Whole sequence Superlocus " + str(self.super_locus.given_name), loc='left', weight='bold')
            record.plot(ax=ax1)
            cropped_record.plot_translation(ax=ax2, location=(400, 400), fontdict={'weight': 'bold'})
            cropped_record.plot(ax=ax2, plot_sequence=True)
            ax2.set_title("Sequence detail " + str(self.super_locus.given_name), loc='left', weight='bold')
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


class DrawableGeenuffCollection:
    def __init__(self, list_of_drawable_super_loci, coordinate_piece):
        def _init_graphic_features():
            res = []
            for sl in list_of_drawable_super_loci:
                res += sl.graphic_features
            return res

        self.coordinate_piece = coordinate_piece
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




#
# def plot_superloci(data, zoom_coords_dict_by_superloc_id=None, plot_separately=True, output_directory=".",
#                    output_prefix="tmp_"):
#     for item in data:
#         super_loci = item['super_loci']
#         coordinate_piece = item['coordinate_piece']
#
#         # for vis:
#         seq = coordinate_piece['sequence']
#         features = []
#         for sl in super_loci:
#             if plot_separately:
#                 features = []
#
#             sl_transcripts = sl['transcripts']
#             for sl_t in sl_transcripts:
#                 for feature in sl_t["features"]:
#                     # print(feature['is_plus_strand'])
#                     if feature['type'] in ['geenuff_intron']:
#                         pass
#                     elif feature['type'] in ['geenuff_transcript']:
#                         features.append(
#                             GraphicFeature(
#                                 start=feature["start"],
#                                 end=feature["end"],
#                                 strand=convert_strand_info(feature['is_plus_strand']),  # todo
#                                 color=color(feature['type']), label=feature['given_name']),
#                         )
#                     else:
#                         features.append(
#                             GraphicFeature(
#                                 start=feature["start"],
#                                 end=feature["end"],
#                                 strand=convert_strand_info(feature['is_plus_strand']),  # todo
#                                 color=color(feature['type']), label=feature['type']),
#                         )
#             record = GraphicRecord(sequence=seq, features=features)
#             if zoom_coords_dict_by_superloc_id:
#                 zoom = zoom_coords_dict_by_superloc_id
#                 zoom_start, zoom_end = zoom[str(sl["id"])]  # coordinates of the "detail"
#                 cropped_record = record.crop((zoom_start, zoom_end))
#                 fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 2))
#                 ax1.set_title("Whole sequence Superlocus " + str(sl["id"]), loc='left', weight='bold')
#                 record.plot(ax=ax1)
#
#                 cropped_record.plot_translation(ax=ax2, location=(400, 400),
#                                                 fontdict={'weight': 'bold'})
#                 cropped_record.plot(ax=ax2, plot_sequence=True)
#                 ax2.set_title("Sequence detail", loc='left', weight='bold')
#
#                 fig.savefig(output_prefix + str(sl["id"]) + '.png', bbox_inches='tight')
#             else:
#                 pass
#

