from geenuff.base import types


class NotGeenuffDeserializableException(Exception):
    pass


class Transcript:
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

    def __init__(self, is_fully_contained=None, transcripts=None, type=None, id=None, overlaps=None, given_name=None):
        self.is_fully_contained = is_fully_contained
        self.transcripts = transcripts
        self.type = type
        self.id = id
        self.overlaps = overlaps
        self.given_name = given_name


class CoordinatePiece:
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
    @classmethod
    def from_dct(cls, dct):
        def _init_super_loci():
            print(dct)
            super_loci = dct["super_loci"]
            res = []
            for sl in super_loci:
                res.append(SuperLocus.from_dct(sl))
            return res

        return cls(
            super_loci=_init_super_loci(),
            coordinate_piece=CoordinatePiece.from_dct(dct=dct["coordinate_piece"])
        )

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

