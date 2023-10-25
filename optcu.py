
#_NOT_IN_PAPER [
########### Seed random gen ############
import numpy as np
import random
seed = 8
np.random.seed(seed)
random.seed(seed)
#_NOT_IN_PAPER ]

# PARENT LATTICE - create
if 1==0:
    ########### ParentLattice ############
    from ase.build import fcc111, add_adsorbate

    pristine = fcc111('Cu', a=3.59, size=(1,1,3)) # 3-atomic-layer Cu(111) slab
    add_adsorbate(pristine,'X',1.7,position='fcc') # Hollow fcc vacancy site
    pristine.center(vacuum=10.0, axis=2) # add vacuum along z-axis

    # Note: CELL is imported with the name "clusterx"
    from clusterx.parent_lattice import ParentLattice

    symbols = [['Cu'],['Cu'],['Cu','Pt'],['X','O']]
    platt = ParentLattice(pristine, symbols=symbols)
    #platt.get_sublattice_types(pretty_print=True)
    ######################################

    #_NOT_IN_PAPER [
    platt.serialize(fname = "optcu_PARENT_LATTICE.json")
    #_NOT_IN_PAPER ]    

# SUPER CELL - create   
if 1==0:
    ########### SuperCell ############
    from clusterx.super_cell import SuperCell
    scell = SuperCell(platt,[[4,0],[-2,4]])

    #from clusterx.visualization import juview
    #juview(scell)
    ######################################

    #_NOT_IN_PAPER [
    scell.serialize(fname="optcu_SUPER_CELL.json")
    #_NOT_IN_PAPER ]

# PARENT LATTICE and SUPER CELL - read
if 1==1:
    #_NOT_IN_PAPER [
    from clusterx.parent_lattice import ParentLattice
    from clusterx.super_cell import SuperCell
    platt = ParentLattice(json_db_filepath="optcu_PARENT_LATTICE.json")
    scell = SuperCell(json_db_filepath="optcu_SUPER_CELL.json")
    #_NOT_IN_PAPER ]

# STRUCTURES SET - create
if 1==1:    
    ########### RandomStructures ############
    from clusterx.structures_set import StructuresSet
    sset = StructuresSet(platt)
    for i in range(50):
        rnd_str = scell.gen_random_structure()
        sset.add_structure(rnd_str)
    #juview(sset,n=4)
    ######################################

    #_NOT_IN_PAPER [
    sset.serialize(path="optcu_STRUCTURES_SET_1.json",overwrite=True)
    #_NOT_IN_PAPER ]

# STRUCTURES SET - read
if 1==0:
    #_NOT_IN_PAPER [
    from clusterx.structures_set import StructuresSet
    sset = StructuresSet(json_db_filepath="optcu_STRUCTURES_SET_1.json")
    #_NOT_IN_PAPER ]

# CONCENTRATIONS
if 1==1:
    fc_dict = sset.get_structure(0).get_fractional_concentrations()
    fc_list = sset.get_concentrations(site_type=2,sigma=1)

    #_NOT_IN_PAPER [    
    print(fc_dict)
    print(fc_list)
    print(np.unique(np.array(fc_list),return_counts=True))
    #_NOT_IN_PAPER ]

# TOTAL ENERGY - compute non relaxed
if 1==1:
    ########### Total energy unrelaxed ############
    from ase.calculators.emt import EMT
    sset.set_calculator(EMT()) # Assign EMT() calculator to StructuresSet.

    # The total energy of every structure in sset is evaluated.
    sset.calculate_property(prop_name="e_tot")
    ######################################

    #_NOT_IN_PAPER [
    energies = sset.get_property_values("e_tot")
    print(energies)
    print(np.unique(energies,return_counts=True))


    sset.set_property_values(property_name="e_tot", property_vals=sset.get_property_values("e_tot"))
    #_NOT_IN_PAPER ]

# REFERENCES - compute
if 1==1:
    #_NOT_IN_PAPER [
    from optcu_functions import references
    eCu_bulk, ePt_bulk, eO, eCu_slab, infostr = references()
    with open("optcu_REFERENCES.txt", "w") as f:
        f.write(infostr)
    print(infostr)
    #_NOT_IN_PAPER ]

if 1==1:
    #_NOT_IN_PAPER [
    from optcu_functions import total_energy

    e1, rel_str, e0, structure = total_energy(scell.get_pristine_structure())
    print(e0,e1)
    #_NOT_IN_PAPER [

if 1==1:
    ################ Adsorption energy relaxed ################
    from optcu_functions import ads_energy
    sset.calculate_property(prop_name="e_ads",prop_func=ads_energy)
    ######################################

    #_NOT_IN_PAPER [
    print(sset.get_property_values("e_ads"))

    sset.set_property_values(property_name="e_ads", property_vals=sset.get_property_values("e_ads"))
    #_NOT_IN_PAPER ]
    

if 1==1:
    #_NOT_IN_PAPER [
    ################ Other properties ################
    xO = lambda i, s : s.get_fractional_concentrations()[0][1]
    xPt = lambda i, s : s.get_fractional_concentrations()[2][1]
    sset.calculate_property(prop_name="conc_O",prop_func=xO)
    sset.calculate_property(prop_name="conc_Pt",prop_func=xPt)
    ######################################
    #_NOT_IN_PAPER ]


if 1==1:
    ################ Visualization: e ads rnd ################
    from clusterx.visualization import plot_property_vs_concentration
    plot_property_vs_concentration(
        sset,
        site_type=2,
        property_name ="e_ads",
        show_plot = False,
        #yaxis_label = "$E_{\\textrm{ads}}(\\bm{\\sigma})$ [arb. units]",
        yaxis_label = "Energy of adsorption [a.u.]",
        fig_fname="optcu_PLOT_1.png",
        data_fname="optcu_PLOT_1.txt"
    )
    ######################################

if 1==0:
    ############### Clusters Pool ###############
    from clusterx.clusters.clusters_pool import ClustersPool
    cpool = ClustersPool(platt,
    			 npoints=[1,2,3,4],
    			 radii=[0,-1,-1,4.3],
    			 super_cell=scell, method=1)
    
    cpool.serialize(db_name="optcu_CPOOL.json")
    cpool.print_info()
    ######################################

if 1==1:
    #_NOT_IN_PAPER [
    from clusterx.clusters.clusters_pool import ClustersPool
    cpool = ClustersPool(json_db_filepath="optcu_CPOOL.json")
    #_NOT_IN_PAPER ]

if 1==1:
    from clusterx.correlations import CorrelationsCalculator
    corrcal = CorrelationsCalculator("trigonometric", platt, cpool)

    # Compute full correlations matrix for sset
    x = corrcal.get_correlation_matrix(sset)

    # Extract ab initio property values
    p = sset.get_property_values("e_ads")
    
    #_NOT_IN_PAPER [
    corrcal.serialize("optcu_CCALC.pickle","pickle")
    np.savez('optcu_XP', x=x, p=p)
    
    # Compute correlations for a structure
    #corrs = corrcal.get_cluster_correlations(sset[3])
    #_NOT_IN_PAPER ]

if 1==0:
    #_NOT_IN_PAPER [
    from clusterx.correlations import CorrelationsCalculator

    corrcal = CorrelationsCalculator(filepath="optcu_CCALC.pickle")
    xp = np.load('optcu_XP.npz')
    x = xp["x"]
    p = xp["p"]
    #_NOT_IN_PAPER ]

if 1==1:
    from clusterx.estimators.estimator_factory import EstimatorFactory
    from clusterx.model import Model

    # Build and fit a ridge-regression estimator of scikit-learn
    re = EstimatorFactory.create("skl_Ridge", alpha=1e-8) 
    re.fit(x,p)
    ce_model = Model(corrcal, "e_ads", estimator = re) # Create a CE model
    ce_model.report_errors(sset)

    #_NOT_IN_PAPER [
    ce_model.serialize(filepath="optcu_CE_MODEL_1.pickle", fmt="pickle")
    #_NOT_IN_PAPER ]
    
if 1==0:
    #_NOT_IN_PAPER [
    from clusterx.model import Model
    ce_model = Model(filepath = "optcu_CE_MODEL_1.pickle")
    #_NOT_IN_PAPER ]
    

if 1==1:
    from clusterx.model import ModelBuilder
    mb_2 = ModelBuilder(
        selector_type="subsets_cv", 
        selector_opts={"clusters_sets": "size"}, 
        estimator_type="skl_Ridge", 
        estimator_opts={"alpha":1e-8,"fit_intercept":True})

    ce_model_2 = mb_2.build(sset, cpool, "e_ads")
    ce_model_2.report_errors(sset)

    #_NOT_IN_PAPER [
    ce_model_2.serialize(filepath="optcu_CE_MODEL_2.pickle", fmt="pickle")
    #mb_2.serialize(filepath="optcu_MODELBDR_2.pickle") -- fails with "TypeError: cannot pickle '_io.BufferedWriter' object"
    opt_cpool_2=mb_2.get_opt_cpool()
    opt_cpool_2.serialize(db_name="optcu_CPOOL_OPT_2.json")
    #_NOT_IN_PAPER ]
    
if 1==0:
    #_NOT_IN_PAPER [
    from clusterx.model import Model, ModelBuilder
    ce_model_2 = Model(filepath = "optcu_CE_MODEL_2.pickle")
    mb_2 = ModelBuilder(filepath="optcu_MODELBDR_2.pickle")
    opt_cpool_2 = mb_2.get_opt_cpool()
    #_NOT_IN_PAPER ]
   
if 1==1:
    ################ Visualization: e ads rnd ################
    from clusterx.visualization import plot_property_vs_concentration

    plot_property_vs_concentration(
        sset,
        site_type=2,
        cemodel=ce_model_2,
        property_name ="e_ads",
        yaxis_label = "Energy of adsorption [a.u.]",
        show_plot = False,
        fig_fname="optcu_PLOT_2.png",
        data_fname="optcu_PLOT_2")
    ######################################

if 1==1:
    from clusterx.visualization import plot_optimization_vs_number_of_clusters
    plot_optimization_vs_number_of_clusters(
        mb_2.get_selector(),
        ymin=-0.0005,
        ymax=0.024,
        xmin=0,
        xmax=50,
        yaxis_label = "Energy [a.u.]",
        show_plot = False,
        fig_fname="optcu_PLOT_3.png",
        data_fname="optcu_PLOT_3"
    )

if 1==1:
    from clusterx.model import ModelBuilder
    mb_3 = ModelBuilder(
        selector_type="lasso_cv", 
        selector_opts={'sparsity_max': 1e-2,'sparsity_min': 1e-6}, 
        estimator_type="skl_Ridge", 
        estimator_opts={"alpha":1e-8,"fit_intercept":True}
    )

    ce_model_3 = mb_3.build(sset, cpool, "e_ads")
    ce_model_3.report_errors(sset)

    #_NOT_IN_PAPER [
    ce_model_3.serialize(filepath="optcu_CE_MODEL_3.pickle", fmt="pickle")
    #mb_3.serialize(filepath="optcu_MODELBDR_3.pickle") # Fails with TypeError: cannot pickle '_io.BufferedWriter' object
    opt_cpool_3 = mb_3.get_opt_cpool()
    opt_cpool_3.serialize(db_name="optcu_CPOOL_OPT_3.json")
    #_NOT_IN_PAPER ]
    
if 1==0:
    #_NOT_IN_PAPER [
    from clusterx.model import Model, ModelBuilder
    ce_model_3 = Model(filepath = "optcu_CE_MODEL_3.pickle")
    mb_3 = ModelBuilder(filepath="optcu_MODELBDR_3.pickle")
    opt_cpool_3 = mb_3.get_opt_cpool()
    #_NOT_IN_PAPER

if 1==1:
    from clusterx.visualization import plot_optimization_vs_number_of_clusters
    plot_optimization_vs_number_of_clusters(
        mb_3.get_selector(),
        ymin=-0.001,
        ymax=0.02,
        xmin=0,
        xmax=125,
        yaxis_label = "Energy [a.u.]",
        show_plot = False,
        fig_fname="optcu_PLOT_4.png",
        data_fname="optcu_PLOT_4"
    )

if 1==1:
    from clusterx.visualization import plot_optimization_vs_sparsity
    plot_optimization_vs_sparsity(
        mb_3.get_selector(),
        xmin=1e-6,
        xmax=1e-3,
        ymin=-0.0003,
        ymax=0.009,
        yaxis_label = "Energy [a.u.]",
        show_plot = False,
        fname="optcu_PLOT_5.png"
    )
    
