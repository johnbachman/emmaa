from os.path import abspath, dirname, join
from datetime import datetime
from nose.plugins.attrib import attr
from emmaa.answer_queries import QueryManager, format_results, \
    is_query_result_diff
from emmaa.queries import Query, DynamicProperty, get_agent_from_trips
from emmaa.model_tests import ModelManager
from emmaa.tests.test_db import _get_test_db
from emmaa.tests.test_model import create_model


test_query = {'type': 'path_property', 'path': {'type': 'Activation',
              'subj': {'type': 'Agent', 'name': 'BRAF',
                       'db_refs': {'HGNC': '1097'}},
              'obj': {'type': 'Agent', 'name': 'MAPK1',
                      'db_refs': {'HGNC': '6871'}},
                      'obj_activity': 'activity'}}
simple_query = 'BRAF activates MAPK1.'
query_object = Query._from_json(test_query)
dyn_ag = get_agent_from_trips('active MAP2K1')
dyn_query = DynamicProperty(dyn_ag, 'eventual_value', 'high')
open_qj = {'type': 'open_search_query',
            'entity': {'type': 'Agent', 'name': 'BRAF',
                       'db_refs': {'HGNC': '1097'}},
            'entity_role': 'subject', 'stmt_type': 'Activation'}
open_query = Query._from_json(open_qj)
test_response = {
    '3801854542': {
        'path': 'BRAF → MAP2K1 → MAPK1',
        'edge_list': [
            {'edge': 'BRAF → MAP2K1',
             'stmts': [
                 ['/evidence?stmt_hash=-23078353002754841&source='
                  'model_statement&model=test&date=2020-01-01',
                  'BRAF activates MAP2K1.', '']]},
            {'edge': 'MAP2K1 → MAPK1',
             'stmts': [
                 ['/evidence?stmt_hash=-34603994586320440&source='
                  'model_statement&model=test&date=2020-01-01',
                  'Active MAP2K1 activates MAPK1.', '']]}]}}
query_not_appl = {'2413475507': 'Query is not applicable for this model'}
fail_response = {
    '521653329': 'No path found that satisfies the test statement'}
# Create a new EmmaaModel and ModelManager for tests instead of depending
# on S3 version
test_model = create_model()
test_mm = ModelManager(test_model)
test_mm.date_str = '2020-01-01-00-00-00'
test_email = 'tester@test.com'


def test_format_results():
    date = datetime.now()
    results = [('test', query_object, 'pysb', test_response, date),
               ('test', query_object, 'signed_graph', fail_response, date),
               ('test', query_object, 'unsigned_graph', test_response, date)]
    formatted_results = format_results(results)
    assert len(formatted_results) == 1
    qh = query_object.get_hash_with_model('test')
    assert qh in formatted_results
    assert formatted_results[qh]['model'] == 'test'
    assert formatted_results[qh]['query'] == simple_query
    assert isinstance(formatted_results[qh]['date'], str)
    assert formatted_results[qh]['pysb'] == [
        'Pass', [test_response['3801854542']]]
    assert formatted_results[qh]['pybel'] == [
        'n_a', 'Model type not supported']
    assert formatted_results[qh]['signed_graph'] == [
        'Fail', fail_response['521653329']]
    assert formatted_results[qh]['unsigned_graph'] == [
        'Pass', [test_response['3801854542']]]


@attr('nonpublic')
def test_answer_immediate_query():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    query_hashes = qm.answer_immediate_query(
        test_email, 1, query_object, ['test'], subscribe=False)[
            'path_property']
    assert query_hashes == [35683418474694258], query_hashes
    results = qm.retrieve_results_from_hashes(query_hashes)
    assert len(results) == 1
    assert query_hashes[0] in results
    result_values = results[query_hashes[0]]
    assert result_values['model'] == 'test'
    assert result_values['query'] == simple_query
    assert isinstance(result_values['date'], str)
    assert result_values['pysb'] == ['Pass', [test_response['3801854542']]], \
        result_values['pysb']
    assert result_values['pybel'] == ['Pass', [test_response['3801854542']]]
    assert result_values['signed_graph'][0] == 'Pass'
    assert result_values['unsigned_graph'][0] == 'Pass'


@attr('nonpublic')
def test_immediate_dynamic():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    query_hashes = qm.answer_immediate_query(
        test_email, 1, dyn_query, ['test'], subscribe=False)[
            'dynamic_property']
    assert query_hashes == [-27775603206605897], query_hashes
    results = qm.retrieve_results_from_hashes(query_hashes, 'dynamic_property')
    assert len(results) == 1, results
    assert query_hashes[0] in results
    result_values = results[query_hashes[0]]
    assert result_values['model'] == 'test'
    assert result_values['query'] == 'Active MAP2K1 is eventually high.'
    assert isinstance(result_values['date'], str)
    assert result_values['result'] == [
        'Pass', 'Satisfaction rate is 100% after 2 simulations.']
    assert isinstance(result_values['image'], str)


def test_immediate_open():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    query_hashes = qm.answer_immediate_query(
        test_email, 1, open_query, ['test'], subscribe=False)[
            'open_search_query']
    assert query_hashes == [-13552944417558866], query_hashes
    results = qm.retrieve_results_from_hashes(query_hashes)
    assert len(results) == 1
    assert query_hashes[0] in results
    result_values = results[query_hashes[0]]
    assert result_values['model'] == 'test'
    assert result_values['query'] == 'What does BRAF activate?'
    assert isinstance(result_values['date'], str)
    for mc_type in ['pysb', 'pybel', 'signed_graph', 'unsigned_graph']:
        assert result_values[mc_type][0] == 'Pass'
        assert isinstance(result_values[mc_type][1], list)
        assert test_response['3801854542']['path'] in [
            res['path'] for res in result_values[mc_type][1]]


@attr('nonpublic')
def test_answer_get_registered_queries():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    # Put all types of queries in db, answer together
    qm.db.put_queries(test_email, 1, query_object, ['test'],
                      subscribe=True)
    qm.db.put_queries(test_email, 1, dyn_query, ['test'],
                      subscribe=True)
    qm.db.put_queries(test_email, 1, open_query, ['test'],
                      subscribe=True)
    qm.answer_registered_queries('test')
    # Retrieve results for path query
    results = qm.get_registered_queries(test_email, 'path_property')
    qh = query_object.get_hash_with_model('test')
    assert qh in results
    assert len(results) == 1
    assert results[qh]['model'] == 'test'
    assert results[qh]['query'] == simple_query
    assert isinstance(results[qh]['date'], str)
    assert results[qh]['pysb'] == ['Pass', [test_response['3801854542']]], \
        (results[qh]['pysb'], test_response['3801854542'])
    assert results[qh]['pybel'] == ['Pass', [test_response['3801854542']]]
    assert results[qh]['signed_graph'][0] == 'Pass'
    assert results[qh]['unsigned_graph'][0] == 'Pass'
    # Retrieve results for dynamic query
    results = qm.get_registered_queries(test_email, 'dynamic_property')
    qh = dyn_query.get_hash_with_model('test')
    assert qh in results
    assert results[qh]['model'] == 'test'
    assert results[qh]['query'] == 'Active MAP2K1 is eventually high.'
    assert isinstance(results[qh]['date'], str)
    assert results[qh]['result'] == [
        'Pass', 'Satisfaction rate is 100% after 2 simulations.']
    assert isinstance(results[qh]['image'], str)
    # Retrieve results for dynamic query
    results = qm.get_registered_queries(test_email, 'open_search_query')
    qh = open_query.get_hash_with_model('test')
    assert qh in results
    assert results[qh]['model'] == 'test'
    assert results[qh]['query'] == 'What does BRAF activate?'
    assert isinstance(results[qh]['date'], str)
    for mc_type in ['pysb', 'pybel', 'signed_graph', 'unsigned_graph']:
        assert results[qh][mc_type][0] == 'Pass'
        assert isinstance(results[qh][mc_type][1], list)
        assert test_response['3801854542']['path'] in [
            res['path'] for res in results[qh][mc_type][1]]


def test_is_diff():
    assert not is_query_result_diff(query_not_appl, query_not_appl)
    assert is_query_result_diff(test_response, query_not_appl)


@attr('nonpublic')
def test_report_one_query():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    # Using results from db
    qm.db.put_queries(test_email, 1, query_object, ['test'],
                      subscribe=True)
    qm.db.put_results('test', [(query_object, 'pysb', test_response),
                               (query_object, 'pysb', query_not_appl)])
    str_msg = qm.get_report_per_query('test', query_object)[0]
    assert str_msg
    assert 'A new result to query' in str_msg
    assert 'Query is not applicable for this model' in str_msg
    assert 'BRAF → MAP2K1 → MAPK1' in str_msg
    # String report given two responses explicitly
    str_msg = qm.make_str_report_one_query(
        'test', query_object, 'pysb', test_response, query_not_appl)
    assert str_msg
    assert 'A new result to query' in str_msg, str_msg
    assert 'Query is not applicable for this model' in str_msg
    assert 'BRAF → MAP2K1 → MAPK1' in str_msg
    assert simple_query in str_msg


@attr('nonpublic')
def test_report_files():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    qm.db.put_queries(test_email, 1, query_object, ['test'],
                      subscribe=True)
    qm.db.put_results('test', [(query_object, 'pysb', query_not_appl)])
    results = qm.db.get_results(test_email, latest_order=1)
    qm.make_str_report_per_user(results,
                                filename='test_query_delta.txt')
    report_file = join(dirname(abspath(__file__)), 'test_query_delta.txt')
    with open(report_file, 'r') as f:
        msg = f.read()
    assert msg
    assert 'This is the first result to query' in msg, msg
    assert 'Query is not applicable for this model' in msg
    qm.db.put_results('test', [(query_object, 'pysb', test_response)])
    results = qm.db.get_results(test_email, latest_order=1)
    qm.make_str_report_per_user(results,
                                filename='new_test_query_delta.txt')
    new_report_file = join(dirname(abspath(__file__)),
                           'new_test_query_delta.txt')
    with open(new_report_file, 'r') as f:
        msg = f.read()
    assert msg
    assert 'A new result to query' in msg
    assert 'BRAF → MAP2K1 → MAPK1' in msg


@attr('nonpublic')
def test_user_query_delta():
    db = _get_test_db()
    qm = QueryManager(db=db, model_managers=[test_mm])
    # Using results from db
    qm.db.put_queries(test_email, 1, query_object, ['test'],
                      subscribe=True)
    qm.db.put_results('test', [(query_object, 'pysb', test_response)])
    qm.db.put_results('test', [(query_object, 'pysb', query_not_appl)])
    str_rep, html_rep = qm.get_user_query_delta(user_email=test_email)
    assert str_rep, print(str_rep)
    assert html_rep, print(html_rep)
    assert '</html>' in html_rep
    assert '</body>' in html_rep
    assert 'pysb' in html_rep
    assert 'unsubscribe' in html_rep
