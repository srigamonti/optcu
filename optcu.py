
import numpy as np
import random
seed = 8
np.random.seed(seed)
random.seed(seed)

def example_ce_optcu(task):

    if task == 1: # PARENT LATTICE - create
        ########### ParentLattice (Listing 1) ############
        from ase.build import fcc111, add_adsorbate

        pristine = fcc111('Cu', a=3.59, size=(1,1,3)) # 3-atomic-layer Cu(111) slab
        add_adsorbate(pristine,'X',1.7,position='fcc') # Hollow fcc vacancy site
        pristine.center(vacuum=10.0, axis=2) # add vacuum along z-axis

        # Note: CELL is imported with the name "clusterx"
        from clusterx.parent_lattice import ParentLattice

        symbols = [['Cu'],['Cu'],['Cu','Pt'],['X','O']]
        platt = ParentLattice(pristine, symbols=symbols)
        #platt.get_sublattice_types(pretty_print=True)
        ##################################################

        platt.serialize(fname = "optcu_PARENT_LATTICE.json")
      
    if task == 2: # SUPER CELL - create 
        from clusterx.parent_lattice import ParentLattice
        platt = ParentLattice(json_db_filepath="optcu_PARENT_LATTICE.json")

        ########### SuperCell (Listing 2) ############
        from clusterx.super_cell import SuperCell
        scell = SuperCell(platt,[[4,0],[-2,4]])

        #from clusterx.visualization import juview
        #juview(scell)
        ##############################################

        scell.serialize(fname="optcu_SUPER_CELL.json")

    if task in [3, 5, 9, 10]: # PARENT LATTICE and SUPER CELL - read
        from clusterx.parent_lattice import ParentLattice
        from clusterx.super_cell import SuperCell
        platt = ParentLattice(json_db_filepath="optcu_PARENT_LATTICE.json")
        scell = SuperCell(json_db_filepath="optcu_SUPER_CELL.json")
    
    if task == 3: # STRUCTURES SET - create
        ########### RandomStructures (Listing 3) ############
        from clusterx.structures_set import StructuresSet
        sset = StructuresSet(platt)
        for i in range(50):
            rnd_str = scell.gen_random_structure()
            sset.add_structure(rnd_str)
        #juview(sset,n=4)
        #####################################################

        sset.serialize(path="optcu_STRUCTURES_SET_1.json",overwrite=True)
    
    if task in [4,5,6,7,8,10,11,12,13,15]: # STRUCTURES SET - read
        from clusterx.structures_set import StructuresSet
        sset = StructuresSet(json_db_filepath="optcu_STRUCTURES_SET_1.json")
    
    if task == 4: # TOTAL ENERGY - compute non relaxed
        ########### Total energy unrelaxed (Listing 4) ############
        from ase.calculators.emt import EMT
        sset.set_calculator(EMT()) # Assign EMT() calculator to StructuresSet.

        # The total energy of every structure in sset is evaluated.
        sset.calculate_property(prop_name="e_tot")
        ###########################################################

        energies = sset.get_property_values("e_tot")
        print(energies)
        print(np.unique(energies,return_counts=True))

        sset.set_property_values(property_name="e_tot", property_vals=sset.get_property_values("e_tot"))

    if task == 5: # Compute reference energies
        from optcu_functions import references
        eCu_bulk, ePt_bulk, eO, eCu_slab, infostr = references(scell.get_pristine_structure())
        with open("optcu_REFERENCES.txt", "w") as f:
            f.write(infostr)
        print(infostr)

    if task == 6:
        ################ Adsorption energy relaxed (Listing 5) ################
        from optcu_functions import ads_energy
        sset.calculate_property(prop_name="e_ads",prop_func=ads_energy)
        #######################################################################

        print(sset.get_property_values("e_ads"))

        sset.set_property_values(property_name="e_ads", property_vals=sset.get_property_values("e_ads"))        

    if task == 7: # Other example properties
        xO = lambda i, s : s.get_fractional_concentrations()[0][1]
        xPt = lambda i, s : s.get_fractional_concentrations()[2][1]
        sset.calculate_property(prop_name="conc_O",prop_func=xO)
        sset.calculate_property(prop_name="conc_Pt",prop_func=xPt)

    if task == 8:
        #### Visualization: Property versus concentration (Listing 6 and Figure 8) ####
        from clusterx.visualization import plot_property_vs_concentration
        plot_property_vs_concentration(
            sset,
            site_type=2,
            property_name ="e_ads",
            show_plot = False,
            yaxis_label = "Energy of adsorption [a.u.]",
            fig_fname="optcu_PLOT_1.png",
            data_fname="optcu_PLOT_1.txt"
        )
        ################################################################################

    if task == 9:
        ############### Clusters Pool (Listing 7) ###############
        from clusterx.clusters.clusters_pool import ClustersPool
        cpool = ClustersPool(platt,
                    npoints=[1,2,3,4],
                    radii=[0,-1,-1,4.3],
                    super_cell=scell, method=1)
        
        cpool.serialize(db_name="optcu_CPOOL.json")
        cpool.print_info()
        ##########################################################

    if task in [10,12,15]: # CLUSTERS POOL - read 
        from clusterx.clusters.clusters_pool import ClustersPool
        cpool = ClustersPool(json_db_filepath="optcu_CPOOL.json")

    if task ==10:
        ############### Input data (Listing 8) ###############
        from clusterx.correlations import CorrelationsCalculator
        corrcal = CorrelationsCalculator("trigonometric", platt, cpool)

        # Compute full correlations matrix for sset
        x = corrcal.get_correlation_matrix(sset)

        # Extract ab initio property values
        p = sset.get_property_values("e_ads")
        #######################################################

        corrcal.serialize("optcu_CCALC.pickle","pickle")
        np.savez('optcu_XP', x=x, p=p)
        
        # Compute correlations for a structure
        #corrs = corrcal.get_cluster_correlations(sset[3])

    if task == 11: # INPUT DATA - read
        from clusterx.correlations import CorrelationsCalculator

        corrcal = CorrelationsCalculator(filepath="optcu_CCALC.pickle")
        xp = np.load('optcu_XP.npz')
        x = xp["x"]
        p = xp["p"]

    if task == 11: # CE MODE with ridge regression - create
        ###### Create CE model with ridge regression (Listing 9) ######
        from clusterx.estimators.estimator_factory import EstimatorFactory
        from clusterx.model import Model

        # Build and fit a ridge-regression estimator of scikit-learn
        re = EstimatorFactory.create("skl_Ridge", alpha=1e-8) 
        re.fit(x,p)
        ce_model = Model(corrcal, "e_ads", estimator = re) # Create a CE model
        ce_model.report_errors(sset)
        ################################################################

        ce_model.serialize(filepath="optcu_CE_MODEL_1.pickle", fmt="pickle")
        
    if task == -1: # CE MODE with ridge regression - read
        from clusterx.model import Model
        ce_model = Model(filepath = "optcu_CE_MODEL_1.pickle")
        

    if task == 12: # CE MODEL with subset selection - create
        ###### Create CE model with subset selection (Listing 10) ######
        from clusterx.model import ModelBuilder
        mb_2 = ModelBuilder(
            selector_type="subsets_cv", 
            selector_opts={"clusters_sets": "size"}, 
            estimator_type="skl_Ridge", 
            estimator_opts={"alpha":1e-8,"fit_intercept":True})

        ce_model_2 = mb_2.build(sset, cpool, "e_ads")
        ce_model_2.report_errors(sset)
        ################################################################

        ce_model_2.serialize(filepath="optcu_CE_MODEL_2.pickle", fmt="pickle")
        mb_2.serialize(filepath="optcu_MODELBDR_2.pickle")
        opt_cpool_2=mb_2.get_opt_cpool()
        opt_cpool_2.serialize(db_name="optcu_CPOOL_OPT_2.json")
        
    if task in [13,14]: # CE MODEL with subset selection - read
        from clusterx.model import Model, ModelBuilder
        ce_model_2 = Model(filepath = "optcu_CE_MODEL_2.pickle")
        mb_2 = ModelBuilder(filepath="optcu_MODELBDR_2.pickle")
        opt_cpool_2 = mb_2.get_opt_cpool()
    
    if task == 13:
        #### Visualization: Property versus concentration (Figure 13) ####
        from clusterx.visualization import plot_property_vs_concentration

        plot_property_vs_concentration(
            sset,
            site_type=2,
            cemodel=ce_model_2,
            property_name ="e_ads",
            yaxis_label = "Energy of adsorption [a.u.]",
            show_plot = False,
            fig_fname="optcu_PLOT_2.png",
            data_fname="optcu_PLOT_2"
        )
        ###################################################################

    if task == 14:
        #### Visualization: Optimization by cluster selector (Figure 14) ####
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
        ###################################################################

    if task == 15:
        ###### Create CE model with LASSO (Listing 11) ######
        from clusterx.model import ModelBuilder
        mb_3 = ModelBuilder(
            selector_type="lasso_cv", 
            selector_opts={'sparsity_max': 1e-2,'sparsity_min': 1e-6}, 
            estimator_type="skl_Ridge", 
            estimator_opts={"alpha":1e-8,"fit_intercept":True}
        )

        ce_model_3 = mb_3.build(sset, cpool, "e_ads")
        ce_model_3.report_errors(sset)
        ######################################################

        ce_model_3.serialize(filepath="optcu_CE_MODEL_3.pickle", fmt="pickle")
        mb_3.serialize(filepath="optcu_MODELBDR_3.pickle")
        opt_cpool_3 = mb_3.get_opt_cpool()
        opt_cpool_3.serialize(db_name="optcu_CPOOL_OPT_3.json")
        
    if task == 16: # Read CE model with LASSO
        from clusterx.model import Model, ModelBuilder
        ce_model_3 = Model(filepath = "optcu_CE_MODEL_3.pickle")
        mb_3 = ModelBuilder(filepath="optcu_MODELBDR_3.pickle")
        opt_cpool_3 = mb_3.get_opt_cpool()

    if task == 16:
        #### Visualization: Optimization by cluster selector (Figure 15) ####
        from clusterx.visualization import plot_optimization_vs_sparsity
        plot_optimization_vs_sparsity(
            mb_3.get_selector(),
            #xmin=1e-6,
            #xmax=1e-3,
            #ymin=-0.0003,
            #ymax=0.009,
            yaxis_label = "Energy [a.u.]",
            show_plot = False,
            fname="optcu_PLOT_5.png"
        )
        ###################################################################

def print_info(task):
    if task == 1: print("Create parent lattice (Listing 1)")
    if task == 2: print("Create super cell (Listing 2)")
    if task == 3: print("Create structures set (Listing 3)")
    if task == 4: print("Compute total energy, non relaxed (Listing 4)")
    if task == 5: print("Compute reference energies")
    if task == 6: print("Compute adsorption energy, relaxed (Listing 5)")
    if task == 7: print("Consider other example properties")
    if task == 8: print("Visualization: Property versus concentration (Listing 6 and Figure 8)")
    if task == 9: print("Create clusters pool (Listing 7)")
    if task == 10: print("Create input matrix X and vector of properties P (Listing 8)")
    if task == 11: print("Create CE model with ridge regression (Listing 9)")
    if task == 12: print("Create CE model with subset selection (Listing 10)")
    if task == 13: print("Visualization: Property versus concentration (Figure 13)")
    if task == 14: print("Visualization: Optimization by cluster selector (Figure 14)")
    if task == 15: print("Create CE model with LASSO (Listing 11)")
    if task == 16: print("Visualization: Optimization by cluster selector (Figure 15)")

if __name__ == "__main__":
    """Example CE of O-Pt/Cu(111)

    task:
        1: Create parent lattice (Listing 1)
        2: Create super cell (Listing 2)
        3: Create structures set (Listing 3)
        4: Compute total energy, non relaxed (Listing 4)
        5: Compute reference energies
        6: Compute adsorption energy, relaxed (Listing 5)
        7: Consider other example properties
        8: Visualization: Property versus concentration (Listing 6 and Figure 8)
        9: Create clusters pool (Listing 7)
        10: Create input matrix X and vector of properties P (Listing 8)
        11: Create CE model with ridge regression (Listing 9)
        12: Create CE model with subset selection (Listing 10)
        13: Visualization: Property versus concentration (Figure 13)
        14: Visualization: Optimization by cluster selector (Figure 14)
        15: Create CE model with LASSO (Listing 11)
        16: Visualization: Optimization by cluster selector (Figure 15)
    """
    import sys
    task = sys.argv[1]
    print_info(int(task))
    example_ce_optcu(int(task))     
