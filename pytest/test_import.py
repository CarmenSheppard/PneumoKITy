# Testing script for import
# Python 3.8 using Pytest module
# Carmen Sheppard March 2020

import pytest

import pandas as pd
import os
from Database_tools.import_from_excel import add_variant, add_serotypevariants, add_serotype

@pytest.fixture
def sero_df():
    """
    Creates a serotype dataframe for testing
    """
    df = pd.read_excel(os.path.join("pytest","data","sero_data.xlsx"),
                       dtype={"predicted_pheno": str, "subtypes": bool,
                              "alt_vars": str, "var1": str, "serotype_1": str,
                              "position": float})
    df = df.dropna(how="all")

    # remove leading/trailing whitespaces
    df = df.apply(lambda x: x.str.strip() if x.dtype == object else x)
    return df

@pytest.fixture
def variant_df():
    """
    Creates a variant dataframe for testing
    """
    df = pd.read_excel(os.path.join("pytest","data","variant_data.xlsx"),
                       dtype={"predicted_pheno": str, "subtypes": bool,
                              "alt_vars": str, "var1": str, "serotype_1": str,
                              "position": float})
    df = df.dropna(how="all")

    # remove leading/trailing whitespaces
    df = df.apply(lambda x: x.str.strip() if x.dtype == object else x)
    return df

@pytest.fixture
def serovariant_df():
    """
    Creates a serotype variant dataframe for testing
    """
    df = pd.read_excel(os.path.join("pytest","data","serovar_data.xlsx"),
                       dtype={"predicted_pheno": str, "subtypes": bool,
                              "alt_vars": str, "var1": str, "serotype_1": str,
                              "position": float})
    df = df.dropna(how="all")

    # remove leading/trailing whitespaces
    df = df.apply(lambda x: x.str.strip() if x.dtype == object else x)
    return df

@pytest.mark.dbimport
def test_sero(sero_df, capsys):
    """ test import of serovar"""
    add_serotype(sero_df,"pytest/data/test_db")
    captured = capsys.readouterr()
    assert captured.err == ''
    assert captured.out != ''

@pytest.mark.dbimport
def test_variant(variant_df, capsys):
    """
    test import of variant
    """
    add_variant(variant_df,"pytest/data/test_db")
    # Use capsys built in pytest feature to capture standard out and std error
    captured = capsys.readouterr()
    assert captured.err == ''
    assert captured.out != ''

@pytest.mark.dbimport
def test_serovar(serovariant_df, capsys):
    """ test import of serovar"""
    add_serotypevariants(serovariant_df,"pytest/data/test_db")
    captured = capsys.readouterr()
    assert captured.err == ''
    assert captured.out == ''