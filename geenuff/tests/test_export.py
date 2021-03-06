import os
import pytest
from ..applications.importer import ImportController
from ..applications.exporters.sequence import FastaExportController
from ..applications.exporters.lengths import LengthExportController
from ..applications.exporters.json import JsonExportController, FeatureJsonable, TranscriptJsonable, SuperLocusJsonable
from ..applications.exporter import MODES
from geenuff.base import orm, types
import json

EXPORTING_DB = 'testdata/exporting.sqlite3'


@pytest.fixture(scope="module", autouse=True)
def prepare_and_cleanup():
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')

    if not os.path.exists(EXPORTING_DB):
        controller = ImportController(database_path='sqlite:///' + EXPORTING_DB)
        controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True,
                              genome_args={'species': 'dummy'})
    yield
    os.remove(EXPORTING_DB)


def seq_len_controllers(mode, longest=False):
    econtroller = FastaExportController(db_path_in='sqlite:///' + EXPORTING_DB, longest=longest)
    econtroller.prep_ranges(range_function=MODES[mode], genomes=None, exclude=None)

    lcontroller = LengthExportController(db_path_in='sqlite:///' + EXPORTING_DB, longest=longest)
    lcontroller.prep_ranges(range_function=MODES[mode], genomes=None, exclude=None)
    return econtroller, lcontroller


def compare2controllers(expect, econtroller, lcontroller):
    iexpect = iter(expect)
    ilengths = iter(lcontroller.export_ranges)
    for i, grp in enumerate(econtroller.export_ranges):
        exp = next(iexpect)
        lgrp = next(ilengths)
        seq = econtroller.get_seq(grp)
        length = lcontroller.get_length(lgrp)
        try:
            assert ''.join(seq) == exp[2]
            assert length == exp[1]
        except AssertionError as e:
            print(exp[0])
            print(i)
            raise e


def test_get_intron_seqs():
    """checks that expected intron sequences are produced"""
    expect = [("Chr1:195000-199000:1543-2114", 572,
               "GTACTTCCAAATCTTCAATTTTGATTCTAAAGATTGGTCCTTTTACTCTGTTTCTCAATT"
               "TGAGTTTTAGGTATTCTTTGATTTTGTATTGGTTTCATTCTAAATATTCATCCTTTACTC"
               "AACTTCTAGATAAGGGATTTAGGTATTCTCAAATTTCCGATTTGATTCCTTTACTCGTTT"
               "CTAGATTGGGGTTTTAGGAATTACCAGTTGGGGGTTTTGCAATTTGCGTAATCAAAGAAT"
               "TTTATTTGTTGTATTGCTTGGTATTGAAGTTTGTCTCTGTTTCTCTACCTCGTCATGTAA"
               "TGTGCTTAGATCCATTAAGTAAATGCTTGTGGATATTTATGTAGATGGTTAAGAGTGATC"
               "GTGATCAGAGTCCTTCTCTTATTTAACTGCATTGCCTGTGAGTTGTGGTCCTGAAGGTTG"
               "TTGTTATTATTGAATTCTATGTATGTATAGATTATGTCATTGGTCTCATGTGGTTTTTAT"
               "GGGTAACGTCTTTACTAATAATAGCACTATGCTTCTGGATTTTGATCTATGTGATCTGTA"
               "ACATTTCTAGTTGGTGTGTCTTTGATTGCCAG"),
              ("Chr1:195000-199000:2220-2299", 80,
               "GTATATATACCGCTGCTCGTATCTCTTTTCCGGTGTTACAAAAGCGATGTCGTGACCTAA"
               "TGCTGGGTTCGTTACTATAG"),
              ("Chr1:195000-199000:2423-2500", 78,
               "GTAAGTCTGGAATAGCTTTTGAGTTGTCCTCTATGTTTATAAGCTATTGTTGTGTGTAAA"
               "CCTTTGTTATATCTGTAG"),
              ("Chr1:195000-199000:2672-2774", 103,
               "GTAAACTATTAAACTCATTAACTCTCTCCTGCAATCTGCAAGGCAGTCTTTAGGAATGTG"
               "AATATTAGGAAATAACTTTTACTTTGTGGGTTGATTTGTTTAG"),
              ("Chr1:195000-199000:2905-2973", 69,
               "GTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCG"
               "TGGCTGTAG"),
              ("Chr1:195000-199000:3184-3266", 83,
               "GTATACAAAACAATTTGCCTTTACGTTTTTACATTTCTTTAAGAGTTTGAAACATGTCTA"
               "AAGCTGGGATAATATTTTTGCAG")]
    expect += expect[:4] + expect[5:6]

    # test data has been simplified from an augustus run that previously resulted in erroneous masks
    # and from a partial gene model in the Rcommumnis genome
    econtroller, lcontroller = seq_len_controllers('introns')
    # expect [(fa_id (samtools faidx), length, sequence)]
    # and the fa_id is just for knowing where it went wrong if it does
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 11
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('introns', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 6
    compare2controllers(expect[:6], econtroller, lcontroller)


def test_get_exons():
    """checks that individual expected exon sequences are produced"""
    expect = [("Chr1:195000-199000:780-1542", 763,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATT"
               "TGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGA"
               "AAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCT"
               "AAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTT"
               "CGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGAT"
               "CTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATC"
               "GAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAA"
               "CAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGT"
               "CGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGC"
               "TTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAA"
               "CTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAG"),
              ("Chr1:195000-199000:2115-2219", 105,
               "GGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAA"
               "GAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAG"),
              ("Chr1:195000-199000:2300-2422", 123,
               "GATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGT"
               "TGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAA"
               "AGG"),
              ("Chr1:195000-199000:2501-2671", 171,
               "ATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAG"
               "CAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCT"
               "CTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGG"),
              ("Chr1:195000-199000:2775-2904", 130,
               "TTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
               "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCC"
               "TCAAGGAGAG"),
              ("Chr1:195000-199000:2974-3183", 210,
               "ACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTGTT"
               "GATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCT"
               "TTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATATG"
               "CTTGAAGCCGAAGATCTACTCTATCGCGAT"),
              ("Chr1:195000-199000:3267-3684", 418,
               "GAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCT"
               "GCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAA"
               "AAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCG"
               "GTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTT"
               "ACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACT"
               "TATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGT"
               "CGATGGATCAAACTTATCGCTTTGGACGAAATGGACCACAATGATTCTTTTTTAGCTC"),
              ("Chr1:195000-199000:812-1542", 731,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATG"
               "CTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAA"
               "GAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTAC"
               "GACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTA"
               "TGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCC"
               "CTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCC"
               "ATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGAT"
               "CGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCC"
               "GGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTG"
               "AGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGT"
               "TACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTG"
               "CTTAACAATAG")]
    expect += expect[1:4]
    expect += [("Chr1:195000-199000:2775-3183", 409,
                "TTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
                "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCC"
                "TCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTA"
                "CTATGTGTCGTGGCTGTAGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCG"
                "AAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAA"
                "ACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAAT"
                "GGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGAT"),
               ("Chr1:195000-199000:3267-3635", 369,
                "GAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCT"
                "GCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAA"
                "AAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCG"
                "GTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTT"
                "ACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACT"
                "TATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGT"
                "CGATGGATC"),
               ("27488:1-3000:2777-2962", 186,
                "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCA"
                "CCTTAAGACTGTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATT"
                "TAGATGGGATGATACTGCAACAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCT"
                "TAATAA")]

    econtroller, lcontroller = seq_len_controllers('exons')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 14
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('exons', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 8
    expect = expect[:7] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_mRNA():
    """tests that exons are propperly spliced into mRNA sequence"""
    expect = [("AT1G01540.2.TAIR10", 1920,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTC"
               "TTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGC"
               "GAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGT"
               "ACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGT"
               "GATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGC"
               "CGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGA"
               "GCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACG"
               "GCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTC"
               "TGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGAT"
               "TGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAGGGGTCAA"
               "GCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTT"
               "TAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCA"
               "ATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATATTATACTGGGG"
               "ATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAA"
               "GCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTC"
               "TGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATACGCTTGCACC"
               "GGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATCACTGGAAGAA"
               "ACCCGGTTGATTATAGTCGGCCTCAAGGAGAGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAA"
               "CCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTG"
               "TTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATA"
               "TGCTTGAAGCCGAAGATCTACTCTATCGCGATGAACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAG"
               "ACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAG"
               "CAAAGATGAAAAAAGAGAGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTC"
               "CGGTATATTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATT"
               "AACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACTTATTAGATATTGTAATAT"
               "GTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGTCGATGGATCAAACTTATCGCTTTGGACG"
               "AAATGGACCACAATGATTCTTTTTTAGCTC"),
              ("AT1G01540.1.TAIR10", 1908,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAA"
               "GCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAA"
               "CGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCG"
               "AAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCA"
               "TCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGC"
               "CTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGA"
               "GTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCC"
               "GGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTT"
               "TGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCA"
               "AAGTCGCCGTCAAGAACTTGCTTAACAATAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGT"
               "CATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATG"
               "CTCGTGTATGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCC"
               "CGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGG"
               "TCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCT"
               "AAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGG"
               "GAACTTTCGGTTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAG"
               "CTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAG"
               "GTATTTGATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGTAGA"
               "CAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAAT"
               "ACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGAT"
               "GCGAACAAGAGACCTAAAATGGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGATG"
               "AACGCCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGA"
               "AAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAAAAGAGAGTCACTTGGGTTAAG"
               "TGATTTCCACACGACCATTATATTTCATTTTTTTTTCCGGTATATTATTGTCTCACTTACATAATATTTT"
               "CATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACG"
               "ATCACATTCACATGTCACTTATTAGATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAG"
               "GTGATTGGTCGATGGATC"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCACCTTAAGACT"
               "GTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATTTAGATGGGATGATACTGCAA"
               "CAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCTTAATAA")]

    econtroller, lcontroller = seq_len_controllers('mRNA')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('mRNA', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_pre_mRNA():
    """tests returning sequence of raw transcript (transcription start-end)"""
    expect = [("Chr1:195000-199000:780-3684", 2905,
               "AAGCCTTTCTCTTTAAATTCGTTATCGTTTTTTTTATTTTATCAATTTAATCTTTTTATT"
               "TGTTTTGTTCTTCCTCGATTCAACACTCGATGCTGTGACAAAGCCCAGATTTTGGCCGGA"
               "AAGTTTTGTGTGTTTCCGGCGAAGAATGCGAAGAAACCGGAAGAATCTGTAACGGAATCT"
               "AAGCTAAAAGTTAAAGTACGATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTT"
               "CGAAACCGACATCGATTTTTGGTCTCCGGCTATGGGTCGTGATCGGAATCTTACTTGGAT"
               "CTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTTAACTTCTCGCCGGAAAAATC"
               "GAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATTTCAAAGGAGA"
               "TTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAA"
               "CAGCGAGTGCAAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGT"
               "CGCATCTTGGATGGGGCCGATGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGC"
               "TTTGTGAAGAGAATGTAATCGGAGAAGGTGGTTACGGGATTGTGTATCGTGGCATTTTAA"
               "CTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAATAGGTACTTCCAAATCTTCA"
               "ATTTTGATTCTAAAGATTGGTCCTTTTACTCTGTTTCTCAATTTGAGTTTTAGGTATTCT"
               "TTGATTTTGTATTGGTTTCATTCTAAATATTCATCCTTTACTCAACTTCTAGATAAGGGA"
               "TTTAGGTATTCTCAAATTTCCGATTTGATTCCTTTACTCGTTTCTAGATTGGGGTTTTAG"
               "GAATTACCAGTTGGGGGTTTTGCAATTTGCGTAATCAAAGAATTTTATTTGTTGTATTGC"
               "TTGGTATTGAAGTTTGTCTCTGTTTCTCTACCTCGTCATGTAATGTGCTTAGATCCATTA"
               "AGTAAATGCTTGTGGATATTTATGTAGATGGTTAAGAGTGATCGTGATCAGAGTCCTTCT"
               "CTTATTTAACTGCATTGCCTGTGAGTTGTGGTCCTGAAGGTTGTTGTTATTATTGAATTC"
               "TATGTATGTATAGATTATGTCATTGGTCTCATGTGGTTTTTATGGGTAACGTCTTTACTA"
               "ATAATAGCACTATGCTTCTGGATTTTGATCTATGTGATCTGTAACATTTCTAGTTGGTGT"
               "GTCTTTGATTGCCAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGG"
               "GCGAGTACGACACAAGAATCTTGTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAG"
               "GTATATATACCGCTGCTCGTATCTCTTTTCCGGTGTTACAAAAGCGATGTCGTGACCTAA"
               "TGCTGGGTTCGTTACTATAGGATGCTCGTGTATGACTTTGTCGACAATGGTAATTTGGAG"
               "CAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAAT"
               "ATTATACTGGGGATGGCCAAAGGGTAAGTCTGGAATAGCTTTTGAGTTGTCCTCTATGTT"
               "TATAAGCTATTGTTGTGTGTAAACCTTTGTTATATCTGTAGATTGGCGTATCTACACGAG"
               "GGTCTTGAGCCAAAAGTTGTTCATCGGGATATTAAATCAAGCAATATCTTACTTGATCGC"
               "CAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGCTCTTGGGGTCTGAGAGCAGT"
               "TATGTGACTACTCGTGTGATGGGAACTTTCGGGTAAACTATTAAACTCATTAACTCTCTC"
               "CTGCAATCTGCAAGGCAGTCTTTAGGAATGTGAATATTAGGAAATAACTTTTACTTTGTG"
               "GGTTGATTTGTTTAGTTATGTAGCACCAGAATACGCTTGCACCGGAATGTTAAACGAGAA"
               "GAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATCACTGGAAGAAACCCGGT"
               "TGATTATAGTCGGCCTCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTCTATGTTT"
               "TTGTAAATGGTCTTACTATGTGTCGTGGCTGTAGACAAATCTAGTGGATTGGCTTAAATC"
               "AATGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATC"
               "CTCAAAGGCTCTTAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAA"
               "CAAGAGACCTAAAATGGGTCATATCATACATATGCTTGAAGCCGAAGATCTACTCTATCG"
               "CGATGTATACAAAACAATTTGCCTTTACGTTTTTACATTTCTTTAAGAGTTTGAAACATG"
               "TCTAAAGCTGGGATAATATTTTTGCAGGAACGCCGAACAACAAGGGACCATGGAAGCCGC"
               "GAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGG"
               "CATCATCAGCAAAAGCAAAGATGAAAAAAGAGAGTCACTTGGGTTAAGTGATTTCCACAC"
               "GACCATTATATTTCATTTTTTTTTCCGGTATATTATTGTCTCACTTACATAATATTTTCA"
               "TTCTTTGGCTCTTCGAGTCTTGGCCTTACGAAATTAACATTCTTCTTCTTTGTAACCTTT"
               "TGTCAACGATCACATTCACATGTCACTTATTAGATATTGTAATATGTAATGTTTGGACCG"
               "ACGTCGAAGCATAAGCAGGTGATTGGTCGATGGATCAAACTTATCGCTTTGGACGAAATG"
               "GACCACAATGATTCTTTTTTAGCTC"),
              ("Chr1:195000-199000:812-3635", 2824,
               "TTTATTTTATCAATTTAATCTTTTTATTTGTTTTGTTCTTCCTCGATTCAACACTCGATG"
               "CTGTGACAAAGCCCAGATTTTGGCCGGAAAGTTTTGTGTGTTTCCGGCGAAGAATGCGAA"
               "GAAACCGGAAGAATCTGTAACGGAATCTAAGCTAAAAGTTAAAGTACGATGTCGGTGTAC"
               "GACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGCTA"
               "TGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCC"
               "CTCTGCTTAACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCC"
               "ATCGCTACACCGCCGATTTCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCT"
               "GTTCCGGCGGAGATCCAGGTCGATATCGGGAAGATCGAGCATCGAGTGGTGTTTTCAGAT"
               "CGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGCAAGTGAAACGGCGTCGTATTCC"
               "GGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGATGGTATACTCTG"
               "AGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTGGT"
               "TACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTG"
               "CTTAACAATAGGTACTTCCAAATCTTCAATTTTGATTCTAAAGATTGGTCCTTTTACTCT"
               "GTTTCTCAATTTGAGTTTTAGGTATTCTTTGATTTTGTATTGGTTTCATTCTAAATATTC"
               "ATCCTTTACTCAACTTCTAGATAAGGGATTTAGGTATTCTCAAATTTCCGATTTGATTCC"
               "TTTACTCGTTTCTAGATTGGGGTTTTAGGAATTACCAGTTGGGGGTTTTGCAATTTGCGT"
               "AATCAAAGAATTTTATTTGTTGTATTGCTTGGTATTGAAGTTTGTCTCTGTTTCTCTACC"
               "TCGTCATGTAATGTGCTTAGATCCATTAAGTAAATGCTTGTGGATATTTATGTAGATGGT"
               "TAAGAGTGATCGTGATCAGAGTCCTTCTCTTATTTAACTGCATTGCCTGTGAGTTGTGGT"
               "CCTGAAGGTTGTTGTTATTATTGAATTCTATGTATGTATAGATTATGTCATTGGTCTCAT"
               "GTGGTTTTTATGGGTAACGTCTTTACTAATAATAGCACTATGCTTCTGGATTTTGATCTA"
               "TGTGATCTGTAACATTTCTAGTTGGTGTGTCTTTGATTGCCAGGGGTCAAGCAGAGAAGG"
               "AATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTTGTTAGGCTTT"
               "TAGGGTATTGCGTGGAAGGTGCATACAGGTATATATACCGCTGCTCGTATCTCTTTTCCG"
               "GTGTTACAAAAGCGATGTCGTGACCTAATGCTGGGTTCGTTACTATAGGATGCTCGTGTA"
               "TGACTTTGTCGACAATGGTAATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAG"
               "CCCGCTAACTTGGGATATACGTATGAATATTATACTGGGGATGGCCAAAGGGTAAGTCTG"
               "GAATAGCTTTTGAGTTGTCCTCTATGTTTATAAGCTATTGTTGTGTGTAAACCTTTGTTA"
               "TATCTGTAGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGATAT"
               "TAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACT"
               "TGCTAAGCTCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGG"
               "GTAAACTATTAAACTCATTAACTCTCTCCTGCAATCTGCAAGGCAGTCTTTAGGAATGTG"
               "AATATTAGGAAATAACTTTTACTTTGTGGGTTGATTTGTTTAGTTATGTAGCACCAGAAT"
               "ACGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCA"
               "TGGAGATAATCACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGGTATTTG"
               "ATAAGCATATTCAATCCTCTCTATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGT"
               "AGACAAATCTAGTGGATTGGCTTAAATCAATGGTGGGAAACCGAAGATCAGAAGAAGTTG"
               "TTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCTTAAACGAGTGTTGCTTGTAG"
               "CTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCATATCATACATA"
               "TGCTTGAAGCCGAAGATCTACTCTATCGCGATGTATACAAAACAATTTGCCTTTACGTTT"
               "TTACATTTCTTTAAGAGTTTGAAACATGTCTAAAGCTGGGATAATATTTTTGCAGGAACG"
               "CCGAACAACAAGGGACCATGGAAGCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGG"
               "TAGTGAAAGTGGTGAGAGCGGTTCACGGCATCATCAGCAAAAGCAAAGATGAAAAAAGAG"
               "AGTCACTTGGGTTAAGTGATTTCCACACGACCATTATATTTCATTTTTTTTTCCGGTATA"
               "TTATTGTCTCACTTACATAATATTTTCATTCTTTGGCTCTTCGAGTCTTGGCCTTACGAA"
               "ATTAACATTCTTCTTCTTTGTAACCTTTTGTCAACGATCACATTCACATGTCACTTATTA"
               "GATATTGTAATATGTAATGTTTGGACCGACGTCGAAGCATAAGCAGGTGATTGGTCGATG"
               "GATC"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCA"
               "CCTTAAGACTGTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATT"
               "TAGATGGGATGATACTGCAACAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCT"
               "TAATAA")]

    econtroller, lcontroller = seq_len_controllers('pre-mRNA')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('pre-mRNA', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_CDS():
    """test production of CDS sequence, ignoring phase"""
    expect = [("AT1G01540.2.TAIR10", 1419,
               "ATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGC"
               "TATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTT"
               "AACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATT"
               "TCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGC"
               "AAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGA"
               "TGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTG"
               "GTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAA"
               "TAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTT"
               "GTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTA"
               "ATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATAT"
               "TATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGAT"
               "ATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGC"
               "TCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATA"
               "CGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATC"
               "ACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGACAAATCTAGTGGATTGGCTTAAATCAA"
               "TGGTGGGAAACCGAAGATCAGAAGAAGTTGTTGATCCGAAAATACCCGAACCACCATCCTCAAAGGCTCT"
               "TAAACGAGTGTTGCTTGTAGCTTTGCGTTGTGTGGATCCTGATGCGAACAAGAGACCTAAAATGGGTCAT"
               "ATCATACATATGCTTGAAGCCGAAGATCTACTCTATCGCGATGAACGCCGAACAACAAGGGACCATGGAA"
               "GCCGCGAGAGACAAGAGACAGCTGTGGTTGCTGCCGGTAGTGAAAGTGGTGAGAGCGGTTCACGGCATCA"
               "TCAGCAAAAGCAAAGATGA"),
              ("AT1G01540.1.TAIR10", 1161,
               "ATGTCGGTGTACGACGCTGCTTTCCTTAATACAGAGCTTTCGAAACCGACATCGATTTTTGGTCTCCGGC"
               "TATGGGTCGTGATCGGAATCTTACTTGGATCTCTAATTGTCATCGCACTCTTTCTTCTCTCCCTCTGCTT"
               "AACTTCTCGCCGGAAAAATCGAAAGCCGAGAGCCGATTTCGCCTCCGCCGCCATCGCTACACCGCCGATT"
               "TCAAAGGAGATTAAAGAGATCGTTCCGGCGCAGAATCAGTCTGTTCCGGCGGAGATCCAGGTCGATATCG"
               "GGAAGATCGAGCATCGAGTGGTGTTTTCAGATCGAGTGTCGAGTGGTGAGAGTAGAGGAACAGCGAGTGC"
               "AAGTGAAACGGCGTCGTATTCCGGTAGCGGGAATTGTGGGCCGGAGGTGTCGCATCTTGGATGGGGCCGA"
               "TGGTATACTCTGAGAGAGCTTGAAGCGGCCACGAATGGGCTTTGTGAAGAGAATGTAATCGGAGAAGGTG"
               "GTTACGGGATTGTGTATCGTGGCATTTTAACTGATGGAACCAAAGTCGCCGTCAAGAACTTGCTTAACAA"
               "TAGGGGTCAAGCAGAGAAGGAATTCAAAGTAGAAGTGGAAGTCATTGGGCGAGTACGACACAAGAATCTT"
               "GTTAGGCTTTTAGGGTATTGCGTGGAAGGTGCATACAGGATGCTCGTGTATGACTTTGTCGACAATGGTA"
               "ATTTGGAGCAATGGATTCACGGTGATGTTGGCGATGTCAGCCCGCTAACTTGGGATATACGTATGAATAT"
               "TATACTGGGGATGGCCAAAGGATTGGCGTATCTACACGAGGGTCTTGAGCCAAAAGTTGTTCATCGGGAT"
               "ATTAAATCAAGCAATATCTTACTTGATCGCCAATGGAATGCTAAGGTTTCGGATTTTGGACTTGCTAAGC"
               "TCTTGGGGTCTGAGAGCAGTTATGTGACTACTCGTGTGATGGGAACTTTCGGTTATGTAGCACCAGAATA"
               "CGCTTGCACCGGAATGTTAAACGAGAAGAGTGATATCTATAGCTTCGGAATACTAATCATGGAGATAATC"
               "ACTGGAAGAAACCCGGTTGATTATAGTCGGCCTCAAGGAGAGGTATTTGATAAGCATATTCAATCCTCTC"
               "TATGTTTTTGTAAATGGTCTTACTATGTGTCGTGGCTGTAG"),
              ("27488.m000034.TIGRR0.1", 186,
               "TGGTTATGAAAATAATAAGAAATTTAAAGTTAAATGTCTACATGAACATGTTGAAAATCACCTTAAGACT"
               "GTAAAAAATGAACGAGCAACTATTTCAACACTTCAAGCTAAAAATGGATTTAGATGGGATGATACTGCAA"
               "CAAAATATGTTTATGATGAACAAATGATGGTAAGCATTCTTAATAA"
               )]
    econtroller, lcontroller = seq_len_controllers('CDS')
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 3
    compare2controllers(expect, econtroller, lcontroller)
    # and now just the longest (1st & 3rd)
    econtroller, lcontroller = seq_len_controllers('CDS', longest=True)
    assert len(econtroller.export_ranges) == len(lcontroller.export_ranges) == 2
    expect = expect[:1] + expect[-1:]
    compare2controllers(expect, econtroller, lcontroller)


def test_get_json_feature():
    controller = JsonExportController(db_path_in='sqlite:///' + EXPORTING_DB)
    f = controller.session.query(orm.Feature).filter(orm.Feature.type == types.GEENUFF_CDS).first()
    print(f.type)
    t = f.transcript_pieces[0].transcript
    coord = controller.session.query(orm.Coordinate).first()
    fh = FeatureJsonable(data=f)
    th = TranscriptJsonable(data=t)
    print(fh.to_jsonable(f, coord, 790, 3000, True, transcript=t))
    #print(json.dumps(th.to_jsonable(t, coord, 790, 3000, True), indent=2))
    print('------')
    sls = controller.session.query(orm.SuperLocus).all()
    #print(len(sls))
    for sl in sls:
        slh = SuperLocusJsonable(sl)
        #print(json.dumps(slh.to_jsonable(sl, coord, 790, 3000, True), indent=2))
    #print(controller.session.query(orm.Genome).all())
    meh = controller.coordinate_range_to_jsonable('dummy', seqid='Chr1:195000-199000', start=1, end=3900, is_plus_strand=True)
    assert len(meh) == 1
    assert meh[0]['coordinate_piece']['seqid'] == 'Chr1:195000-199000'
    assert len(meh[0]['coordinate_piece']['sequence']) == 3899
    assert len(meh[0]['super_loci']) == 1
    print(json.dumps(meh, indent=2))
    # todo, slightly more thorough testing
