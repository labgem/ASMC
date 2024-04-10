import sys
from pathlib import Path

asmc_path = Path(__file__).parents[1]
sys.path.append(str(asmc_path))

from asmc import utils

import pytest

class TestPDB:
    
    @pytest.fixture(scope="class")
    def models(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        m = d / "models.txt"
        
        f1 = d / "id1a.pdb"
        data_f1 = "\n".join([
            "REMARK   1",
            "ATOM      1  N   MET A   1      29.812",
            "ATOM      2  CA  MET A   1      29.297",
            "ATOM      9  N   THR A   2      27.031",
            "ATOM     16  N   ALA A   3      25.078",
            "HETATM   21  N   THR A   4      22.594",
            "ATOM     27  OG1 AAA A   5      23.703",
            "HETATM 4118  O   HOH A 525       9.180"
            ])
        
        f1.write_text(data_f1)
        
        f2 = d / "id2b.pdb"
        data_f2 = "\n".join([
            "ATOM      1  N   TYR A   1      29.812",
            "ATOM      9  N   GLN A   2      27.031",
            "ATOM     16  N   ILE A   3      25.078",
            "ATOM     21  N   MET A   4      22.594",
            "ATOM     27  N   THR A   5      25.078",
            "ATOM   4118  N   ALA A 525       9.180"
            ])

        f2.write_text(data_f2)
        
        m.write_text(f"{f1}\n{f2}")
        
        yield m
    
    def test_get_seq_from_pdb(self, models):
        
        f = Path(models.read_text().split("\n")[0])
        
        seq = utils.get_seq_from_pdb(f)
        assert seq == "MTATX"
        
    def test_read_models(self, models):
        
        d = utils.read_models(models)
        assert d == {"id1a":"MTATX", "id2b":"YQIMTA"}

class TestSeq:
    
    @pytest.fixture(scope="class")
    def fasta(self, tmp_path_factory):
    
        data = ">id2b\nYQIMTA\n>id3c\nMSARGGK\n>id4d\nWQIMTARGL"
        
        d = tmp_path_factory.mktemp("data")
        f = d / "seqs.fasta"
        f.write_text(data)
        
        yield f
    
    def test_read_multi_fasta(self, fasta):
        
        d = utils.read_multi_fasta(fasta)
        assert d == {"id2b":"YQIMTA","id3c":"MSARGGK", "id4d":"WQIMTARGL"}
        
    def test_get_identity(self, fasta):
        
        d = utils.read_multi_fasta(fasta)
        target = d["id4d"]
        del d["id4d"]
        
        best_ref, perc_id = utils.get_identity(d, target)
        assert best_ref == "id2b"
        assert round(perc_id, 2) == 55.56
    
class TestCompareActiveSite:
    
    @pytest.fixture(scope='class')
    def groups_tsv(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        
        f1 = d / "groups1.tsv"
        data_f1 = "\n".join(["idA\tAZERTY\t5",
                             "idB\tQWERTY\t5",
                             "idC\tPEPONV\t6",
                             "idD\tPENOVR\t6",
                             "idE\tPEPOKR\t6"])

        f2 = d / "groups2.tsv"
        data_f2 = "\n".join(["idA\tAZERTY\t5",
                             "idB\tQWERTY\t5",
                             "idC\tPEPOFR\t6",
                             "idE\tPEPOKR\t6",
                             "idF\tQWERTZ\t5"])
        
        f1.write_text(data_f1)
        f2.write_text(data_f2)
        
        yield (f1, f2)
        
    def test_build_comparison_data(self, groups_tsv):
        
        f1 = groups_tsv[0]
        d = utils.build_comparison_data({}, f1)
        assert d["idA"]["f1"] == "AZERTY"
        assert d["idA"]["f2"] is None
        assert d["idA"]["d"] is None
        assert d["idA"]["ref_id"] is None
        assert d["idA"]["ref_seq"] is None
        assert d["idA"]["ref_d1"] is None
        assert d["idA"]["ref_d2"] is None
        assert d["idA"]["ref_pid"] is None
    
        f2 = groups_tsv[1]
        d = utils.build_comparison_data(d, f2, empty=False)

        assert d["idC"]["f1"] == "PEPONV"
        assert d["idC"]["f2"] == "PEPOFR"
        assert d["idD"]["f1"] == "PENOVR"
        assert d["idD"]["f2"] is None
        
    @pytest.fixture(scope="class")
    def identity_data(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        
        f = d / "identity.tsv"
        data = "\n".join(["idB\tidA\t50.0",
                          "idC\tidE\t28.7",
                          "idD\tidE\t31.4"])
        f.write_text(data)
    
        yield f
    
    def test_add_ref_data_to_comparison_data(self, groups_tsv, identity_data):
        
        f1 = groups_tsv[0]
        d = utils.build_comparison_data({}, f1)
        f2 = groups_tsv[1]
        d = utils.build_comparison_data(d, f2, empty=False)
        
        d, s = utils.add_ref_data_to_comparison_data(d, identity_data)
        
        assert s == {"idA", "idE"}
        assert d["idB"]["ref_id"] == "idA"
        assert d["idB"]["ref_pid"] == '50.0'
        assert d["idC"]["ref_id"] == "idE"
        assert d["idC"]["ref_pid"] == '28.7'
        assert d["idF"]["ref_id"] is None
        assert d["idA"]["g1"] == "5"
        assert d["idA"]["g2"] == "5"
    
    def test_LD_two_rows(self):
    
        assert utils.LD_two_rows("sitting", "kiiten") == 3
        
    def test_compute_levenshtein(self, groups_tsv, identity_data):
        
        f1 = groups_tsv[0]
        d = utils.build_comparison_data({}, f1)
        f2 = groups_tsv[1]
        d = utils.build_comparison_data(d, f2, empty=False)
        
        d, s = utils.add_ref_data_to_comparison_data(d, identity_data)
        
        d = utils.compute_levenshtein(d)
        
        assert d["idA"]["d"] == 0
        assert d["idC"]["d"] == 2
        assert d["idD"]["ref_d1"] == 2
        assert d["idD"]["ref_d2"] is None
    
    def test_build_active_site_checking_file(self, groups_tsv, identity_data):
        
        f1 = groups_tsv[0]
        d = utils.build_comparison_data({}, f1)
        f2 = groups_tsv[1]
        d = utils.build_comparison_data(d, f2, empty=False)
        
        d, s = utils.add_ref_data_to_comparison_data(d, identity_data)
        
        d = utils.compute_levenshtein(d)
        
        t = utils.build_active_site_checking_file(d, s)
        t = t.split("\n")
        
        assert t[0] == "ID\tG1\tF1\tG2\tF2\tD\tREF_ID\tPERC_ID\tREF_SEQ\tD_REF_1\tD_REF_2\tNEAR_REF"
        assert t[1] == "idB\t5\tQWERTY\t5\tQWERTY\t0\tidA\t50.0\tAZERTY\t2\t2\tboth"
        assert t[2] == "idC\t6\tPEPONV\t6\tPEPOFR\t2\tidE\t28.7\tPEPOKR\t2\t1\tf2"
        assert t[3] == "idD\t6\tPENOVR\tNone\tNone\tNone\tidE\t31.4\tPEPOKR\t2\tNone\tf1"
        assert t[4] == "idF\tNone\tNone\t5\tQWERTZ\tNone\tNone\tNone\tNone\tNone\tNone\tNone"

class TestExtractAA:
    
    @pytest.fixture
    def data_for_extract_aa(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        f = d / "extract_aa.tsv"
        
        data = "\n".join(["A\tASW\t1",
                          "B\tFTP\t1",
                          "C\tKEL\t2",
                          "D\tRGH\t2",
                          "E\tCYK\t3",
                          "F\tIDN\t4",
                          "G\tMQV\t4"])
        
        f.write_text(data)
        
        bad_f = d / "bad.tsv"
        bad_f.write_text("A\nB\nC")
        
        yield d
        
    def test_extract_aa(self, data_for_extract_aa):
        
        correct_file = data_for_extract_aa / "extract_aa.tsv"
        bad_file = data_for_extract_aa / "bad.tsv"
        
        result = utils.extract_aa(correct_file, 1, "A", None)
        assert result == "A\tASW\t1\n"
        
        result = utils.extract_aa(correct_file, 1, "F", None)
        assert result == "B\tFTP\t1\n"
        
        result = utils.extract_aa(correct_file, 1, "K", None)
        assert result == "C\tKEL\t2\n"
        
        result = utils.extract_aa(correct_file, 3, "N", None)
        assert result == "F\tIDN\t4\n"
        
        result = utils.extract_aa(correct_file, 3, "aromatic", None)
        assert result == "A\tASW\t1\n"
        
        result = utils.extract_aa(correct_file, 2, "acidic", None)
        assert result == "C\tKEL\t2\nF\tIDN\t4\n"
        
        result = utils.extract_aa(correct_file, 3, "basic", None)
        assert result == "D\tRGH\t2\nE\tCYK\t3\n"
        
        result = utils.extract_aa(correct_file, 2, "hydrophobic", None)
        assert result == "D\tRGH\t2\n"
        
        result = utils.extract_aa(correct_file, 1, "polar", None)
        assert result == "E\tCYK\t3\n"
        
        file_format_error = f"'{bad_file}' does not appear to contain at least "
        file_format_error += "2 columns"
        with pytest.raises(utils.FileFormatError, match=file_format_error):
            result = utils.extract_aa(bad_file, 1, "K", None)
            
        position_error = "position must be between 1 and 3, got '5'"
        with pytest.raises(utils.PositionError, match=position_error):
            result = utils.extract_aa(correct_file, 5, "G", None)
            
        amino_acid_type_error = "expected 1-letter amino acid or an amino acid "
        amino_acid_type_error += "type, got 'Z'"
        with pytest.raises(utils.AminoAcidTypeError, match=amino_acid_type_error):
            result = utils.extract_aa(correct_file, 1, 'Z', None)
            
class TestGetUnique:
    
    @pytest.fixture
    def data_for_get_unique(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        f = d / "groups.tsv"
        
        data = "\n".join([
            "A\tAZERTY\t1",
            "B\tAZERTY\t1",
            "C\tAZERTY\t1",
            "D\tQWERTY\t2",
            "E\tQWERTY\t2",
            "F\tQWERTZ\t2",
            "G\tQWERTO\t2",
            "H\tPEPOFR\t3",
            "I\tPEBOVR\t3",
            "J\tPENORR\t3",
            "K\tPETOJR\t3"
        ])
        
        f.write_text(data)
        
        bad_file = d / "bad.tsv"
        
        bad_data = "\n".join([
            "Y\tAAAAAA",
            "Z\tAAAAAA"
        ])
    
        bad_file.write_text(bad_data)
        
        yield d
        
    def test_get_unique(self, data_for_get_unique):
        
        correct_file = data_for_get_unique / "groups.tsv"
        bad_file = data_for_get_unique / "bad.tsv"
        
        with pytest.raises(utils.FileFormatError):
            unique, stats = utils.get_unique(bad_file)
            
        unique, stats = utils.get_unique(correct_file)
        
        print(unique)
        
        expected_unique = {'AZERTY': ('1', {'A', 'B', 'C'}),
                           'QWERTY': ('2', {'D', 'E'}),
                           'QWERTZ': ('2', {'F'}),
                           'QWERTO': ('2', {'G'}),
                           'PEPOFR': ('3', {'H'}),
                           'PEBOVR': ('3', {'I'}),
                           'PENORR': ('3', {'J'}),
                           'PETOJR': ('3', {'K'})}
        
        expected_stats = {
            "1":(1, 3, 0.0), "2":(3, 4, 1.0), "3":(4, 4, 2.0)
        }
        
        assert unique == expected_unique
        assert stats == expected_stats