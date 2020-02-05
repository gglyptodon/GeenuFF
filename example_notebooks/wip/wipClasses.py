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
    # todo: classmethods for overloading


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
    def __init__(self, seqid, id, sequence, start, end):
        self.seqid = seqid
        self.id = id
        self.sequence = sequence
        self.start = start
        self.end = end


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


