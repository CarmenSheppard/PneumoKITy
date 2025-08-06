"""
Python 3.7+
Database SQLalchemy declarative script for ORM.
Creates CTV.db structure
Carmen Sheppard 2020-22
"""

from sqlalchemy import Boolean, Column, String, Integer, ForeignKey
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship, backref

# initiates a base sqlalchemy class for python to SQLite translation
Base = declarative_base()


# Set up tables with no foreign keys
class Group(Base):
    """
    Define Group lookup table to look up group name (human
    friendly)
    """
    __tablename__ = 'genogroup'
    # id assigned by db
    id = Column(Integer, primary_key=True, autoincrement=True)
    # Human friendly group name
    group_name = Column(String, primary_key=False)

    # Set up one to many relationships
    group_id = relationship("Serotype", backref=backref("genogroup", lazy="joined"))
    grp_id = relationship("VariantGroup", backref=backref("genogroup", lazy="joined"))


class Genes(Base):
    """
    Define Genes table - these are the serotype discrimination genes
    gene ids are taken from the operon_genes table
    """
    __tablename__ = 'genes'

    # id assigned by db
    id = Column(Integer, primary_key=True, autoincrement=True)
    # name of gene
    gene_name = Column(String, primary_key=False)

    # relationships- one to many
    gene = relationship("Variants", backref=backref("genes", lazy="joined"))
    gene_id = relationship("OperonGene", backref=backref("genes", lazy="joined"))


# Set up tables dependent on above tables only
class Serotype(Base):
    """
    Define a serotype class for serotype table includes predicted pheno
    and serotype hit.
   """
    __tablename__ = 'serotype'
    # id assigned by db
    id = Column(Integer, primary_key=True, autoincrement=True)
    # Predicted phenotype (Final result)
    predicted_pheno = Column(String, primary_key=False)
    # corresponding hit serotype  from stage 1 screen (may be more than one
    # per predicted pheno)
    serotype_hit = Column(String, primary_key=False)
    # id of the group record from the group table - null for non groups
    group_id = Column(Integer, ForeignKey('genogroup.id'), nullable=True)
    # does this serotype have subtypes (Boolean)
    subtype = Column(Boolean, default=False)

    # Set up one to many relationships
    sero_id = relationship("SerotypeVariants", backref=backref("serotype", lazy="joined"))
    serotype_id = relationship("OperonGene", backref=backref("serotype", lazy="joined"))


class Variants(Base):
    """
    Variants table containing positions and variant/ alternative variants
    and the variant types.
    """
    __tablename__ = 'variants'

    id = Column(Integer, primary_key=True, autoincrement=True)
    var_type = Column(String, primary_key=False)
    gene = Column(Integer, ForeignKey('genes.id'))
    position = Column(Integer, primary_key=False, nullable=True)
    variant = Column(String, primary_key=False)


    # Set up one to many relationships
    variant_id = relationship("SerotypeVariants", backref=backref("variants", lazy="joined"))
    var_id = relationship("VariantGroup", backref=backref("variants", lazy="joined"))


class OperonGene(Base):
    """
    Define operon genes table with reference to serotype_hit in serotype table
    """
    __tablename__ = 'operon_gene'
    # id assigned by db
    id = Column(Integer, primary_key=True, autoincrement=True)
    # id of serotype (1st stage)
    serotype_id = Column(Integer, ForeignKey('serotype.id'))
    # Gene id from genes table
    gene_id = Column(Integer, ForeignKey('genes.id'))

    # Numerical order of the genes in the operon for the hit serotype
    operon_position = Column(Integer, primary_key=False)



# set tables with joins
class VariantGroup(Base):
    """
    Variant/Group join table.
    """
    __tablename__ = 'variant_group'
    id = Column(Integer, primary_key=True, autoincrement=True)
    grp_id = Column('group_id', Integer, ForeignKey('genogroup.id'))
    var_id = Column('variants_id', Integer, ForeignKey('variants.id'))


class SerotypeVariants(Base):
    """
    Define serotype_variant mapping table containing serotypes and their
    associated variant ID.
    """
    __tablename__ = 'serotype_variants'
    id = Column(Integer, primary_key=True, autoincrement=True)
    sero_id = Column(Integer, ForeignKey('serotype.id'))
    variant_id = Column(Integer, ForeignKey('variants.id'))

if __name__ == "__main__":
    #### CREATE DATABASE
    engine = create_engine('sqlite:///../ctvdb/CTV.db')

    # Create all tables
    Base.metadata.create_all(engine)
