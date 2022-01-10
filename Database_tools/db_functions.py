"""
Script to set up DB sessions and other common functions for accessing CTV.db
Python 3.7+
Carmen Sheppard 2020-2021
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from Database_tools.sqlalchemydeclarative import Base


def session_maker(db_path):
    """
    Sets up sessions database MUST be called CTV.db
    :param db_path:  path to database (str)
    :return: DB session object
    """
    engine = create_engine(f'sqlite:///{db_path}/CTV.db')

    # Bind the engine to Base class to allow access through DBSession
    Base.metadata.bind = engine
    # establish db session
    DBSession = sessionmaker(bind=engine)
    session = DBSession()
    return session


def searchlike(text, table, fieldname, session):
    """
    Universal LIKE search function for table/field name,returns all fields in table
    Session must be started prior to running this function
    :param text: search text (string)
    :param table: db table class object
    :param fieldname: dbtable field object
    :param session: Database session object
    :return: query results
    """
    text = '%' + text + '%' # add wildcards for LIKE search to search within string
    results = session.query(table).filter(fieldname\
                                          .like(text)).all()

    # check if query has returned results
    if results:
        return results
    else:
        return ['No results found']


def searchexact(text, table, fieldname, session):
    """
    Universal exact search function for table/field name
    :param text: search text (string)
    :param table: db table class object
    :param fieldname: dbtable field object
    :param session: Database session object
    :return: query results
    """

    results = session.query(table).filter(fieldname == text).all()

    # check if query has returned results
    if results:
        return results
    else:
        return ['No results found']

