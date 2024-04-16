import sys
from pathlib import Path

asmc_path = Path(__file__).parents[1]
sys.path.append(str(asmc_path))

import asmc

import pytest


'''
fonction ayant besoin du fichier de ref en input :

build_ds
build_pocket_text
build_multiple_alignment

fonction ayant besoin d'un directory de sortie:

build_ds
extract_pocket
build_pocket_text
'''

class TestPocketDetection:
    
    @pytest.fixture(scope='class')
    def data_for_pocket_detection(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        ref = d / "references.txt"
        ref_text = f"{d / 'refA.pdb'}"
        
        ref.write_text(ref_text)
        
        yield (d, ref)
    
    def test_build_ds_sys_exception(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        ref = data_for_pocket_detection[1]
        chain = "all"
        
        ref_pdb = outdir / "refA.pdb"
        
        error_msg = f"An error has occured while reading {ref}:\n{ref_pdb}"
        error_msg += " file not found"
        
        with pytest.raises(Exception, match=error_msg):
            ds, text = asmc.build_ds(ref, outdir, chain)
    
    def test_build_ds(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        ref = data_for_pocket_detection[1]
        chain = "all"
        
        ref_pdb = outdir / "refA.pdb"
        ref_pdb.write_text("")
        
        ds, text = asmc.build_ds(ref, outdir, chain)
        
        expected_text = f"HEADER: protein chains\n\n{ref_pdb} *"
        expected_ds = outdir / "data.ds"
        
        assert ds == expected_ds
        assert text == expected_text
        
    def test_extract_pocket_no_prediction(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        
        error_msg = "No predictions file after running p2rank"
        
        with pytest.raises(RuntimeError, match=error_msg):
            asmc.extract_pocket(outdir)
            
    def test_extract_pocket(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        pred_file = outdir / "refA.pdb_predictions.csv"
        
        pred_data = "\n".join([
            "name,rank,score,probability,sas_points,surf_atoms,center_x,"+
            "center_y,center_z,residue_ids,surf_atom_ids",
            "pocket1,1,26.03,0.910,128,78,-3.0958,-14.5148,38.5031,A_1 B_7 A_3"+
            ",4924 4938 4939 4942 4946 4951 4952 5095 5099 5100 5101",
            "pocket2,2,20.03,0.910,128,78,-3.0958,-14.5148,38.5031,C_4 C_5 C_6"+
            ",4924 4938 4939 4942 4946 4951 4952 5095 5099 5100 5101"
        ])
        
        pred_file.write_text(pred_data)
        
        res_dict = asmc.extract_pocket(outdir)
        assert res_dict == {'C':[4,5,6]}
        
    def test_build_pocket_text_exception(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        ref = data_for_pocket_detection[1]
        ref_dict = {}
        query_chain = "Z"
        
        error_msg = "None results for p2rank, this may be due to an incorrect "
        error_msg += f"query chain value : {query_chain}"
        
        with pytest.raises(Exception, match=error_msg):
            asmc.build_pocket_text(ref, ref_dict, outdir, query_chain)
            
    def test_build_pocket_text(self, data_for_pocket_detection):
        
        outdir = data_for_pocket_detection[0]
        ref = data_for_pocket_detection[1]
        ref_dict = {'C':[4,5,6]}
        query_chain = "C"
        
        output, text = asmc.build_pocket_text(ref, ref_dict, outdir, query_chain)
        
        assert output == output
        assert text == "refA,C,4,5,6"
        

class TestStructuralAlignment:
    
    @pytest.fixture(scope='class')
    def data_for_structural_aln(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        ref_file = d / "reference.txt"
        ref_data  = f"{d / 'refA.pdb'}"
        ref_file.write_text(ref_data)
        
        refA_pdb = d / "refA.pdb"
        refA_pdb_data = "\n".join([
            "REMARK    1",
            "ATOM      1  CA  MET C   1      29.297",
            "ATOM      2  CA  ALA C   2      29.297",
            "ATOM      3  C   ALA C   2      29.297",
            "ATOM      4  CA  GLN C   3      29.297",
            "ATOM      5  CA  GLY C   4      29.297",
            "ATOM      6  CA  GLY C   5      29.297",
            "ATOM      7  CA  SER C   6      29.297",
            "ATOM      8  CA  TRP C   7      29.297",
            "ATOM      9  CA  LYS C   8      29.297"
        ])
        
        refA_pdb.write_text(refA_pdb_data)
        
        pocket_file = d / "pocket.csv"
        pocket_data = "\n".join([
            "refA,C,1,2,3,4,5,6,7,8"
        ])
        
        pocket_file.write_text(pocket_data)
        
        pairwise_dir = tmp_path_factory.mktemp("pairwise")
        pair_file_B = pairwise_dir / "protB_-refA.fasta"
        pair_file_C = pairwise_dir / "protC_-refA.fasta"
        pair_file_D = pairwise_dir / "protD_-refA.fasta"
        
        pair_file_B_data = "\n".join([
            ">refA.pdb\nMA--QG-GS-WK",
            ">protB.pdb\nMGNVQG-GS-WK\n#"
        ])
        
        pair_file_C_data = "\n".join([
            ">refA.pdb\n---MA-Q-G--G--S-WK",
            ">protC.pdb\nGGGMA-Q-GVLKG-TS-WR\n#"
        ])
        
        pair_file_D_data = "\n".join([
            ">refA.pdb\nMAQGGSWK",
            ">protD.pdb\nMFNLVLHH\n#"
        ])
        
        pair_file_B.write_text(pair_file_B_data)
        pair_file_C.write_text(pair_file_C_data)
        pair_file_D.write_text(pair_file_D_data)
        
        yield d, pairwise_dir
        
    def test_renumber_residues_sys_exception(self, data_for_structural_aln):
            
        d = data_for_structural_aln[0]
        refA_pdb = d / "refA.pdb"
            
        ref_list = [refA_pdb, "C", ["1","2","10"]]
        
        error_msg = "An error has occured when renumbering the residues "
        error_msg += "of reference. This may caused by a residue number "
        error_msg += "indicated in the pocket file not found in the "
        error_msg += f"'{refA_pdb}' or a duplicated residue number"
        
        with pytest.raises(asmc.RenumberResiduesError, match=error_msg):
            renum = asmc.renumber_residues(ref_list)
            
    def test_renumber_residues(self, data_for_structural_aln):
        
        d = data_for_structural_aln[0]
        refA_pdb = d / "refA.pdb"
            
        ref_list = [refA_pdb, "C", ["1","2","3"]]
        renum = asmc.renumber_residues(ref_list)
        assert renum == [0, 1, 2]
        
    def test_extract_aligned_pos(self, data_for_structural_aln):
        
        d = data_for_structural_aln[0]
        paiwrise_dir = data_for_structural_aln[1]
        refA_pdb = d / "refA.pdb"
        protB_file = paiwrise_dir / "protB_-refA.fasta"
        
            
        ref_list = [refA_pdb, "C", ["1","4","6","7","8"],
                    [0,3,5,6,7]]
        
        text = asmc.extract_aligned_pos("refA", "protB", ref_list,
                                        protB_file, True)
        
        assert text == ">refA\nMGSWK\n>protB\nMGSWK\n"
        
    def test_build_multiple_alignment(self, data_for_structural_aln):
        
        d = data_for_structural_aln[0]
        pairwise_dir = data_for_structural_aln[1]
        ref_file = d / "reference.txt"
        pocket_file = d / "pocket.csv"
        
        text = asmc.build_multiple_alignment(pairwise_dir, ref_file,
                                             pocket_file)
        
        expected_text = ">refA\nMAQGGSWK\n>protB\nMGQGGSWK\n>protC\nMAQGKT-W\n"
        expected_text += ">protD\nMFNLVLHH\n"
        
        assert text == expected_text
        

class TestMultipleSequencesAlignment:
    
    @pytest.fixture(scope='class')
    def data_for_msa(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        msa_one_ref = d / "msa_1_ref.txt"
        msa_multiple_ref = d / "msa_multiple_ref.txt"
        msa_fasta = d / "msa.fasta"
        msa_identity = d / "identity_target_ref.tsv"
        
        msa_one_ref_data = "\n".join([
            "refA,3,5,7,8,9,12",
            f"{msa_fasta}"
        ])
        
        msa_multiple_ref_data = "\n".join([
            "refA,3,5,7,8,9,12",
            "refB,4,6,8,10,12,13",
            f"{msa_identity}\n"
            f"{msa_fasta}"
        ])
        
        msa_identity_data = "\n".join([
            "idC\trefA\t62.50",
            "idD\trefA\t68.75",
            "idE\trefB\t68.75",
            "idF\trefB\t50.00",
            "idG\trefB\t62.50"
        ])
        
        msa_fasta_data = "\n".join([
            ">refA\nMA-TFGLKHASNKRIPM",
            ">refB\nMAETFNFQFLGDQLWS-",
            ">idC\nMG-TLGLRHASNSHIAM",
            ">idD\nMA-TFVLHHASNSRILS",
            ">idE\nMADTFNLQFLGDILFT-",
            ">idF\nMADTFKLQIAGEQYWG-",
            ">idG\nMAYTFNLKKLGNELWS-"
        ])
        
        msa_one_ref.write_text(msa_one_ref_data)
        msa_multiple_ref.write_text(msa_multiple_ref_data)
        msa_identity.write_text(msa_identity_data)
        msa_fasta.write_text(msa_fasta_data)
        
        yield d
        
    def test_search_active_site_in_msa_one_ref(self, data_for_msa):
        
        d = data_for_msa
        msa_one_ref = d / "msa_1_ref.txt"
        
        expected_text = ">refA\nTGKHAK\n>refB\nTNQFLQ\n>idC\nTGRHAS\n>idD\n"
        expected_text += "TVHHAS\n>idE\nTNQFLI\n>idF\nTKQIAQ\n>idG\nTNKKLE\n"
        
        text = asmc.search_active_site_in_msa(msa_one_ref)
        assert  text == expected_text
        
    def test_search_active_site_in_msa_multiple_ref(self, data_for_msa):
        
        d = data_for_msa
        msa_multiple_ref = d / "msa_multiple_ref.txt"
        
        expected_text = ">refA\nTGKHAK\n>idC\nTGRHAS\n>idD\nTVHHAS\n>refB\n"
        expected_text += "TNQLDQ\n>idE\nTNQLDI\n>idF\nTKQAEQ\n>idG\nTNKLNE\n"
        
        text = asmc.search_active_site_in_msa(msa_multiple_ref)
        assert text == expected_text
        
class TestClustering:
    
    @pytest.fixture(scope='class')
    def data_for_clustering(self, tmp_path_factory):
        
        d = tmp_path_factory.mktemp("data")
        
        active_site_aln = d / 'active_site_alignment.fasta'
        active_site_data = ">refA\nTGKHAK\n>idC\nTGRHAS\n>idD\nTVHHAS\n>refB\n"
        active_site_data += "TNQLDQ\n>idE\nTNQLDI\n>idF\nTKQAEQ\n>idG\nTNKLNE\n"
        active_site_data += ">idH\nTG-HAR\n"
        
        active_site_aln.write_text(active_site_data)
        
        matrix = d / "matrix.tsv"
        matrix_data = ('\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW'
                       '\tY\tV\nA\t0\t14\t5\t11\t6\t7\t10\t2\t10\t5\t4\t10\t4\t'
                       '10\t3\t4\t5\t12\t11\t3\nR\t14\t0\t11\t16\t15\t9\t18\t15'
                       '\t10\t15\t14\t6\t11\t15\t14\t10\t9\t18\t17\t14\nN\t5\t'
                       '11\t0\t12\t10\t4\t7\t4\t10\t5\t5\t8\t7\t13\t5\t5\t5\t13'
                       '\t12\t4\nD\t11\t16\t12\t0\t11\t9\t5\t13\t18\t15\t14\t15'
                       '\t11\t16\t11\t8\t7\t17\t17\t13\nC\t6\t15\t10\t11\t0\t8'
                       '\t13\t9\t14\t7\t8\t15\t5\t9\t8\t6\t6\t11\t10\t7\nQ\t7\t'
                       '9\t4\t9\t8\t0\t9\t9\t11\t8\t7\t7\t4\t10\t7\t4\t4\t11\t'
                       '11\t7\nE\t10\t18\t7\t5\t13\t9\t0\t9\t14\t9\t8\t12\t11\t'
                       '17\t8\t10\t9\t16\t15\t8\nG\t2\t15\t4\t13\t9\t9\t9\t0\t'
                       '11\t5\t4\t12\t7\t13\t3\t5\t6\t14\t14\t3\nH\t10\t10\t10'
                       '\t18\t14\t11\t14\t11\t0\t10\t9\t7\t11\t10\t11\t11\t11\t'
                       '9\t8\t9\nI\t5\t15\t5\t15\t7\t8\t9\t5\t10\t0\t1\t11\t4\t9'
                       '\t6\t9\t8\t9\t11\t2\nL\t4\t14\t5\t14\t8\t7\t8\t4\t9\t1'
                       '\t0\t10\t3\t9\t5\t8\t7\t10\t10\t1\nK\t10\t6\t8\t15\t15'
                       '\t7\t12\t12\t7\t11\t10\t0\t10\t16\t10\t10\t10\t15\t15\t'
                       '10\nM\t4\t11\t7\t11\t5\t4\t11\t7\t11\t4\t3\t10\t0\t6\t7'
                       '\t5\t4\t9\t9\t4\nF\t10\t15\t13\t16\t9\t10\t17\t13\t10\t9'
                       '\t9\t16\t6\t0\t13\t11\t9\t4\t5\t10\nP\t3\t14\t5\t11\t8'
                       '\t7\t8\t3\t11\t6\t5\t10\t7\t13\t0\t5\t6\t14\t13\t4\nS\t'
                       '4\t10\t5\t8\t6\t4\t10\t5\t11\t9\t8\t10\t5\t11\t5\t0\t2'
                       '\t12\t12\t7\nT\t5\t9\t5\t7\t6\t4\t9\t6\t11\t8\t7\t10\t4'
                       '\t9\t6\t2\t0\t11\t11\t6\nW\t12\t18\t13\t17\t11\t11\t16'
                       '\t14\t9\t9\t10\t15\t9\t4\t14\t12\t11\t0\t4\t11\nY\t11\t'
                       '17\t12\t17\t10\t11\t15\t14\t8\t11\t10\t15\t9\t5\t13\t12'
                       '\t11\t4\t0\t11\nV\t3\t14\t4\t13\t7\t7\t8\t3\t9\t2\t1\t10'
                       '\t4\t10\t4\t7\t6\t11\t11\t0')
    
        matrix.write_text(matrix_data)
        
        yield d
        
        
    def test_read_alignment(self, data_for_clustering):
        
        d = data_for_clustering
        active_site_aln = d / 'active_site_alignment.fasta'
        
        seq, removed = asmc.read_alignment(active_site_aln)
        
        expected_seq = {'refA':'TGKHAK','idC':'TGRHAS','idD':'TVHHAS',
                        'refB':'TNQLDQ','idE':'TNQLDI','idF':'TKQAEQ',
                        'idG':'TNKLNE'}
        expected_removed = {'idH':'TG-HAR'}
        
        assert removed == expected_removed
        assert seq == expected_seq
    
    def test_read_matrix(self, tmp_path_factory):
        
        matrix_tmp_dir = tmp_path_factory.mktemp("mat")
        matrix = matrix_tmp_dir / "matrix.tsv"
        
        matrix_data = "\n".join([
            "\tA\tR",
            "A\t0\t14",
            "R\t14\t0"
        ])
        
        matrix.write_text(matrix_data)
        
        res = asmc.read_matrix(matrix)
        
        assert res == {"A":{"A":0, "R":14}, "R":{"A":14, "R":0}}
        
    def test_pairwise_score(self, data_for_clustering):
        
        d = data_for_clustering
        matrix = d / "matrix.tsv"
        
        scoring_dict = asmc.read_matrix(matrix)
        
        testA = asmc.pairwise_score(scoring_dict, "ARN", "ADN", [])
        assert testA == (16, '')
        
        testB = asmc.pairwise_score(scoring_dict, "CQFD", "CQFD", [])
        assert testB == (0, '')
        
        testC = asmc.pairwise_score(scoring_dict, "WE-K", "WXVS", [])
        assert testC == (50, '')
        
        testD = asmc.pairwise_score(scoring_dict, "YML", "-GT", [1,3])
        assert testD == (142, '')
        
        testE = asmc.pairwise_score(scoring_dict, "A", "Z", [])
        assert testE == (20, {'Z'})
        
    def test_dissimilarity(self, data_for_clustering):
        
        d = data_for_clustering
        active_site_aln = d / 'active_site_alignment.fasta'
        matrix = d / "matrix.tsv"
        
        scoring_dict = asmc.read_matrix(matrix)
        seq, removed = asmc.read_alignment(active_site_aln)
        
        key_list, data, warn_set = asmc.dissimilarity(seq, scoring_dict, [])
    
        assert key_list == ['refA', 'idC', 'idD', 'refB', 'idE', 'idF', 'idG']
        assert data.shape == (7,7)
        assert data.min() == 0.0
        assert data.max() == 1.0
        assert warn_set == set()
        
    def test_dbscan_clustering(self, data_for_clustering):
        
        d = data_for_clustering
        active_site_aln = d / 'active_site_alignment.fasta'
        matrix = d / "matrix.tsv"
        
        scoring_dict = asmc.read_matrix(matrix)
        seq, removed = asmc.read_alignment(active_site_aln)
        
        key_list, data, warn_set = asmc.dissimilarity(seq, scoring_dict, [])
        
        labels = asmc.dbscan_clustering(data, 0.4, 1, 1)
        assert labels.shape == (7,)
        assert list(labels) == [0,0,0,1,1,1,2]
        
    def test_formatting_output(self, data_for_clustering):
        
        d = data_for_clustering
        active_site_aln = d / 'active_site_alignment.fasta'
        matrix = d / "matrix.tsv"
        
        scoring_dict = asmc.read_matrix(matrix)
        seq, removed = asmc.read_alignment(active_site_aln)
        
        key_list, data, warn_set = asmc.dissimilarity(seq, scoring_dict, [])
        
        labels = asmc.dbscan_clustering(data, 0.4, 1, 1)
        
        expected_G = [('refA', 'TGKHAK', 0), ('idC', 'TGRHAS', 0),
                      ('idD', 'TVHHAS', 0), ('refB', 'TNQLDQ', 1),
                      ('idE', 'TNQLDI', 1), ('idF', 'TKQAEQ', 1),
                      ('idG', 'TNKLNE', 2)]
        
        G = asmc.formatting_output(seq, key_list, labels)
        assert G == expected_G
        
class TestSequenceLogo:
    
    @pytest.fixture(scope="class")
    def group_data(self):
        
        G = [('refA', 'TGKHAK', 0), ('idC', 'TGRHAS', 0),
             ('idD', 'TVHHAS', 0), ('refB', 'TNQLDQ', 1),
             ('idE', 'TNQLDI', 1), ('idF', 'TKQAEQ', 1),
             ('idG', 'TNKLNE', 2)]
        
        yield G
        
    def test_build_fasta(self, group_data):
        
        text = asmc.build_fasta(group_data)
        expected_text = ('>refA\nTGKHAK\n>idC\nTGRHAS\n>idD\nTVHHAS\n>refB\n'
                         'TNQLDQ\n>idE\nTNQLDI\n>idF\nTKQAEQ\n>idG\nTNKLNE\n')
        
        assert text == expected_text
        