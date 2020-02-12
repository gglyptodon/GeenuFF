from geenuff.base import types as geenufftypes

from .wipDrawableExceptions import NotGeenuffDeserializableException, InvalidGeenuffTypeException


def _repr(obj):  # todo
    res = "{}".format(obj.__class__)
    for it in zip(obj.__dict__.keys(), obj.__dict__.values()):
        if it[0] in ["features", "transcripts"]:
            res += "[...]"
        else:
            res += str(it)
    return res


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

    def __repr__(self):
        return _repr(self)


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

    def __repr__(self):
        return _repr(self)


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
        res = "{}".format(self.__class__)
        try:
            res += "seqid:{},id:{},sequence:{}...{},start:{},end:{}".format(self.seqid, self.id, self.sequence[:6], self.sequence[-6:], self.start, self.end)
        except Exception as e:
            print(e)
            res += "seqid:{},id:{},sequence:...,start:{},end:{}".format(self.seqid, self.id, self.start, self.end)
        return res


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
        if type not in [t.value for t in geenufftypes.GeenuffFeature]:
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
        return _repr(self)


class GeenuffCollection:
    """ todo: docstr"""
    @classmethod
    def from_dct(cls, dct):
        def _init_super_loci():
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

    # def get_filtered_by_type(self, geenuff_types):  # TODO
    #     """for debugging, naive approach"""
    #     res = []
    #     for sl in self.super_loci:
    #         for t in sl.transcripts:
    #             for f in t.features:
    #                 if f.type in geenuff_types:
    #                     res.append(sl)
    #     return set(res)
    #
    # def get_filtered_by_given_name(self, given_name):  # TODO
    #     """for debugging, naive approach"""
    #     res = []
    #     for sl in self.super_loci:
    #         for t in sl.transcripts:
    #             for f in t.features:
    #                 if f.given_name in given_name:
    #                     res.append(sl)
    #     return set(res)

    def __repr__(self):
        return _repr(self)


class GeenuffCollectionFilter:
    FILTERS = ['given_name', 'feature_type']
    @staticmethod
    def filter(collection, filter_by, to_match):
        filter_func = _get_filter(filter_by)
        return filter_func(collection, to_match)


def _get_filter(filter_by):
    if filter_by == 'given_name':
        return _filter_by_given_name
    elif filter_by == 'feature_type':
        return _filter_by_feature_type
    else:
        raise ValueError("Available: {}".format(GeenuffCollectionFilter.FILTERS))


def _filter_by_given_name(collection, to_match):
    res = []
    for sl in collection.super_loci:
        for t in sl.transcripts:
            for f in t.features:
                if f.given_name in to_match:
                    res.append(sl)
    return set(res)


def _filter_by_feature_type(collection, to_match):
    res = []
    valid = [e.value for e in geenufftypes.GeenuffFeature]
    print(to_match, valid)
    err = [m for m in to_match if m not in valid]
    if err:
        raise InvalidGeenuffTypeException("{} not in {}".format(err, valid))
    for sl in collection.super_loci:
        for t in sl.transcripts:
            for f in t.features:
                if f.type in to_match:
                    res.append(sl)
    return set(res)



