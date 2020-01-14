import os
import json
from indra.statements import Activation, ActivityCondition, Phosphorylation, \
    Agent, Evidence
from emmaa.analyze_tests_results import ModelRound, TestRound, \
    ModelStatsGenerator, TestStatsGenerator
from emmaa.tests.test_model import indra_stmts as previous_stmts


TestRound.__test__ = False
TestStatsGenerator.__test__ = False


path_here = os.path.abspath(os.path.dirname(__file__))
previous_results_file = os.path.join(path_here, 'previous_results.json')
new_results_file = os.path.join(path_here, 'new_results.json')
previous_test_stats_file = os.path.join(path_here, 'previous_stats.json')
previous_model_stats_file = os.path.join(path_here, 'previous_model_stats.json')
with open(previous_results_file, 'r') as f:
    previous_results = json.load(f)
with open(new_results_file, 'r') as f:
    new_results = json.load(f)
with open(previous_test_stats_file, 'r') as f:
    previous_test_stats = json.load(f)
with open(previous_model_stats_file, 'r') as f:
    previous_model_stats = json.load(f)


new_stmts = previous_stmts + [
    Activation(Agent('BRAF', db_refs={'HGNC': '1097'}),
               Agent('AKT', db_refs={'FPLX': 'AKT'}),
               evidence=[Evidence(text='BRAF activates AKT')]),
    Activation(Agent('AKT', db_refs={'FPLX': 'AKT'},
                     activity=ActivityCondition('activity', True)),
               Agent('MTOR', db_refs={"HGNC": "3942"}),
               evidence=[Evidence(text='AKT activate MTOR')])]


def test_model_round():
    mr = ModelRound(previous_stmts, 'test', '2020-01-01-00-00-00')
    assert mr
    assert mr.get_total_statements() == 2
    assert len(mr.get_stmt_hashes()) == 2
    assert mr.get_statement_types() == [('Activation', 2)]
    assert all(agent_tuple in mr.get_agent_distribution() for agent_tuple in
               [('BRAF', 1), ('MAP2K1', 2), ('MAPK1', 1)])
    assert all((stmt_hash, 1) in mr.get_statements_by_evidence() for stmt_hash
               in mr.get_stmt_hashes())
    mr2 = ModelRound(new_stmts, 'test', '2020-01-02-00-00-00')
    assert mr2
    assert mr2.get_total_statements() == 4
    assert len(mr2.get_stmt_hashes()) == 4
    assert mr2.get_statement_types() == [('Activation', 4)]
    assert all(agent_tuple in mr2.get_agent_distribution() for agent_tuple in
               [('BRAF', 2), ('MAP2K1', 2), ('MAPK1', 1), ('MTOR', 1),
                ('AKT', 2)])
    assert len(mr2.find_delta_hashes(mr, 'statements')['added']) == 2


def test_test_round():
    tr = TestRound(previous_results, 'test', '2020-01-01-00-00-00')
    assert tr
    assert tr.get_total_applied_tests() == 1
    assert tr.get_number_passed_tests() == 1
    assert tr.get_applied_test_hashes() == tr.get_passed_test_hashes()
    assert tr.passed_over_total() == 1.0
    tr2 = TestRound(new_results, 'test', '2020-01-02-00-00-00')
    assert tr2
    assert tr2.get_total_applied_tests() == 2
    assert tr2.get_number_passed_tests() == 2
    assert tr2.get_applied_test_hashes() == tr2.get_passed_test_hashes()
    assert tr2.passed_over_total() == 1.0
    assert len(tr2.find_delta_hashes(tr, 'applied_tests')['added']) == 1
    assert len(tr2.find_delta_hashes(tr, 'passed_tests')['added']) == 1
    assert len(tr2.find_delta_hashes(tr, 'paths')['added']) == 1


def test_model_stats_generator():
    latest_round = ModelRound(new_stmts, 'test', '2020-01-02-00-00-00')
    previous_round = ModelRound(previous_stmts, 'test', '2020-01-01-00-00-00')
    sg = ModelStatsGenerator('test', latest_round=latest_round,
                             previous_round=previous_round,
                             previous_json_stats=previous_model_stats)
    sg.make_stats()
    assert sg.json_stats
    model_summary = sg.json_stats['model_summary']
    assert model_summary['model_name'] == 'test'
    assert model_summary['number_of_statements'] == 4
    assert model_summary['stmts_type_distr'] == [('Activation', 4)]
    assert all(agent_tuple in model_summary['agent_distr'] for
               agent_tuple in [('AKT', 2), ('BRAF', 2), ('MAP2K1', 2),
                               ('MTOR', 1), ('MAPK1', 1)])
    assert len(model_summary['stmts_by_evidence']) == 4
    assert len(model_summary['all_stmts']) == 4
    model_delta = sg.json_stats['model_delta']
    assert len(model_delta['statements_hashes_delta']['added']) == 2
    changes = sg.json_stats['changes_over_time']
    assert changes['number_of_statements'] == [2, 4]
    assert len(changes['dates']) == 2


def test_test_stats_generator():
    latest_round = TestRound(new_results, 'test', '2020-01-02-00-00-00')
    previous_round = TestRound(previous_results, 'test', '2020-01-01-00-00-00')
    sg = TestStatsGenerator('test', latest_round=latest_round,
                            previous_round=previous_round,
                            previous_json_stats=previous_test_stats)
    sg.make_stats()
    assert sg.json_stats
    test_round_summary = sg.json_stats['test_round_summary']
    assert test_round_summary['number_applied_tests'] == 2
    assert len(test_round_summary['all_test_results']) == 2
    assert test_round_summary['pysb']['number_passed_tests'] == 2
    assert test_round_summary['pysb']['passed_ratio'] == 1.0
    assert test_round_summary['pybel']['number_passed_tests'] == 2
    assert test_round_summary['pybel']['passed_ratio'] == 1.0
    assert test_round_summary['signed_graph']['number_passed_tests'] == 2
    assert test_round_summary['signed_graph']['passed_ratio'] == 1.0
    assert test_round_summary['unsigned_graph']['number_passed_tests'] == 2
    assert test_round_summary['unsigned_graph']['passed_ratio'] == 1.0
    tests_delta = sg.json_stats['tests_delta']
    assert len(tests_delta['applied_hashes_delta']['added']) == 1
    assert len(tests_delta['pysb']['passed_hashes_delta']['added']) == 1
    assert len(tests_delta['pybel']['passed_hashes_delta']['added']) == 1
    assert len(tests_delta['signed_graph']['passed_hashes_delta']['added']) == 1
    assert len(tests_delta['unsigned_graph']['passed_hashes_delta']['added']) == 1
    changes = sg.json_stats['changes_over_time']
    assert changes['number_applied_tests'] == [1, 2]
    assert len(changes['dates']) == 2
    assert changes['pysb']['number_passed_tests'] == [1, 2]
    assert changes['pysb']['passed_ratio'] == [1, 1]
    assert changes['pybel']['number_passed_tests'] == [1, 2]
    assert changes['pybel']['passed_ratio'] == [1, 1]
    assert changes['signed_graph']['number_passed_tests'] == [1, 2]
    assert changes['signed_graph']['passed_ratio'] == [1, 1]
    assert changes['unsigned_graph']['number_passed_tests'] == [1, 2]
    assert changes['unsigned_graph']['passed_ratio'] == [1, 1]
