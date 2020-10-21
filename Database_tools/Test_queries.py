"""
Test queries to run outside of PneumoCaT.
"""
from Database_tools.db_functions import searchexact, session_maker
from Database_tools.sqlalchemydeclarative import VariantGroup, Variants, Group, Serotype, Genes, SerotypeVariants

# variables for testing
database = "ctvdb"
results = ["19A", "19F", "19AF", "03"]
groups = []
grp_id = 12
# go through results  from mash and find group in DB
variants = {"wzy": 2}
session = session_maker(database)



# find all variants possible in group
serotypes = session.query(SerotypeVariants).join(Variants). \
    filter(Variants.variant_id == ).all()
genotype = int
for target in variants:
    for group_var in group_vars:












out_res = []
for hit in results:
    # get back
    grp = session.query(Serotype).join(Group).filter(Serotype.serotype_hit == hit).all()
    # if hit is in group get group name
    if grp:
        for g in grp:
            out_res.append(g.genogroup.group_name)
    # if hit is type get phenotype name
    else:
        pheno = session.query(Serotype.predicted_pheno).filter(Serotype.serotype_hit == hit) \
            .all()
        out_res.append(pheno[0][0])
out_res = set(out_res)





# go through results  from mash and find group in DB
for i in results:
    db_result = searchexact(i, Serotype, Serotype.serotype_hit, session)
    if db_result != ["No results found"]:
        for record in db_result:
            groups.append(record.group_id)

print("groups are:")
print(groups)

if len(set(groups)) > 1:
# must contain non-group serotypes in addition
    print("mix")

# if only one hit and that hit is not in a group
elif groups == [None] and len(results) == 1:
    stage1_result = results[0]
    print("type")

# if not meeting above criteria must be a group (even if only 1 hit)
else:
    session = session_maker(database)
    # all groups in list are identical so search on first entry
    db_result = searchexact(groups[0], Group, Group.id, session)
    session.close()



# retrieve group variants from group id using join query
variants = []


vars = session.query(Variants).join(VariantGroup).filter(VariantGroup.grp_id == grp_id).all()



alleles = {}
for item in vars:
    print(item.genes.gene_name)
    print(item.var_type)
    if item.var_type == "allele":
        # query for gene name (file name for FASTA)
        # retrieve genes from var id
        records = searchexact(item.gene, Genes, Genes.id, session)
        # there will only by 1 match as searching unique key.
        for record in records:
            gene = record.gene_name
            # add gene name and variant query results to dictionary
            alleles[gene] = item.gene


print("Allele genes are:")
for allele in alleles:
    print(allele)
var_list2 = []

for i in vars:
    var_list2.append(i.var_type)
var_list2 = "/t".join(var_list2)
print("var_list is")
print(var_list2)
session.close()


print(session.query(Serotype).filter(Serotype.serotype_hit == "19A",
    Serotype.predicted_pheno == "19A",
                                Serotype.group_id == "19A_19AF",
                                 Serotype.subtype == True).all())
print(session.query(Group.id).filter(Group.group_name == "19A_19AF").first())