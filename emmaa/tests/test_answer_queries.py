from datetime import datetime
from nose.plugins.attrib import attr
from emmaa.answer_queries import (
    answer_immediate_query, answer_registered_queries, get_registered_queries,
    format_results, load_model_manager_from_s3)
from emmaa.queries import Query
from emmaa.model_tests import ModelManager
from emmaa.db import get_db
from indra.statements.statements import Activation
from indra.statements.agent import Agent


# Tell nose to not run tests in the imported modules
answer_immediate_query.__test__ = False
answer_registered_queries.__test__ = False
get_registered_queries.__test__ = False
format_results.__test__ = False
load_model_manager_from_s3.__test__ = False
ModelManager.__test__ = False
Activation.__test__ = False
Agent.__test__ = False
get_db.__test__ = False
Query.__test__ = False


test_query = {'type': 'path_property', 'path': {
    'type': 'Activation', 'subj': {'type': 'Agent', 'name': 'BRAF'},
    'obj': {'type': 'Agent', 'name': 'MAPK1'}}}
processed_query = Query._from_json(test_query).to_json()
test_response = {3801854542: [
    ('BRAF activates MAP2K1.',
     'https://db.indra.bio/statements/from_agents?subject=1097@HGNC&object=6840@HGNC&type=Activation&format=html'),
    ('Active MAP2K1 activates MAPK1.',
     'https://db.indra.bio/statements/from_agents?subject=6840@HGNC&object=6871@HGNC&type=Activation&format=html')]}


def test_load_model_manager_from_s3():
    mm = load_model_manager_from_s3('test')
    assert isinstance(mm, ModelManager)


def test_format_results():
    results = [('test', test_query, test_response, datetime.now())]
    formatted_results = format_results(results)
    assert len(formatted_results) == 1
    assert formatted_results[0]['model'] == 'test'
    assert formatted_results[0]['query'] == test_query
    assert isinstance(formatted_results[0]['response'], str)
    assert isinstance(formatted_results[0]['date'], str)


@attr('nonpublic')
def test_answer_immediate_query():
    results = answer_immediate_query('tester@test.com', processed_query, ['test'],
                                     subscribe=False, db_name='test')
    assert len(results) == 1
    assert results[0]['model'] == 'test'
    assert results[0]['query'] == processed_query, (results[0]['query'], processed_query)
    assert isinstance(results[0]['response'], str)
    assert 'BRAF activates MAP2K1.' in results[0]['response'], results[0]['response']
    assert isinstance(results[0]['date'], str)


@attr('nonpublic')
def test_answer_get_registered_queries():
    db = get_db('test')
    db.drop_tables(force=True)
    db.create_tables()
    db.put_queries('tester@test.com', processed_query, ['test'], subscribe=True)
    answer_registered_queries('test', db_name='test')
    results = get_registered_queries('tester@test.com', db_name='test')
    assert len(results) == 1
    assert results[0]['model'] == 'test'
    assert results[0]['query'] == processed_query
    assert isinstance(results[0]['response'], str)
    assert 'BRAF activates MAP2K1.' in results[0]['response'], results[0]['response']
    assert isinstance(results[0]['date'], str)
