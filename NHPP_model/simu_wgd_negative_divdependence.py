from scipy.special import gammaln
from tqdm import trange
from scipy.optimize import differential_evolution
from Bio import Phylo
import click
from rich.logging import RichHandler
import logging
import numpy as np
import pandas as pd

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    This script is to simulate WGD on a timetree under negative diversity-dependence
    Usage: python simu_wgd_negative_divdependence.py simulatewgdontree Maintree.txt WGD_dates_info.tsv -o Simu_WGDates_1000000_Dependence.tsv
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('tree', type=click.Path(exists=True))
@click.argument('date', type=click.Path(exists=True))
@click.option('--output', '-o', default='output_file', show_default=True, help="file name of output file")
def simulatewgdontree(tree,date,output):
    tree = Phylo.read(tree,'newick')
    df = pd.read_csv(date,header=0,index_col=None,sep='\t')
    df = df.drop_duplicates(subset=['WGD ID'])
    WGD_IDs_Ages = {wgdid:age for wgdid,age in zip(df["WGD ID"],df["Consensus Mean"])}
    observed_WGD_dates = df['Consensus Mean'].to_list()
    for indice,i in enumerate(tree.get_nonterminals()): i.name = "Node_{}".format(indice)
    y = lambda x:" ".join(x.split("_")[1:])
    yy = lambda x:{key.replace("α","a").replace("β","b"):value for key,value in x.items()}
    WGD_clades_dics = getwgdlocations(tree)
    WGD_clades_dics = yy(WGD_clades_dics)
    WGDIDs_treenodeids = {y(key):value.name for key,value in WGD_clades_dics.items()}
    Tree_nodes_WGDids = getTree_nodes_WGDIDs(WGDIDs_treenodeids)
    Ages_nodes = get_nodes_ages(tree)
    Ages_nodes = {key:list(value)+[len(Tree_nodes_WGDids.get(key,[]))] for key,value in Ages_nodes.items()}
    branches = nodestobranches(Ages_nodes)
    Num_WGD = df.shape[0]
    a,b,c = fit_nhpp_counts_10000median(branches, Num_WGD)
    sim_WGDs = {}
    simulated_WGD_dates_ = []
    y2 = lambda x:", ".join([str(i) for i in x])
    y3 = lambda x:[-i for i in x]
    Iterations = 1000000
    for time in trange(Iterations):
        simulated_WGD_Dates = []
        for key,value in Ages_nodes.items():
            end_time,start_time = value[1],value[2]
            simulated_WGD_dates = y3(simulation_WGD_dependence(-start_time, -end_time, a, b, c, random_seed=None))
            sim_WGDs[key] = [y2(simulated_WGD_dates)] if sim_WGDs.get(key) is None else sim_WGDs[key] + [y2(simulated_WGD_dates)]
    df = pd.DataFrame.from_dict(sim_WGDs)
    df.index = ["Iteration_{}".format(i) for i in range(df.shape[0])]
    df.to_csv(output,header=True,index=True,sep='\t')

def lambda_t(a, b, c, t):
    return a * np.exp(b * t) + c

def integral_lambda(a, b, c, t0, t1):
    return (a / b) * (np.exp(b * t1) - np.exp(b * t0)) + c * (t1 - t0)

def compute_c_from_expectation(a, b, branches, target_total):
    num = target_total
    denom = sum(br['t1'] - br['t0'] for br in branches)
    adj = sum((a/b) * (np.exp(b*br['t1']) - np.exp(b*br['t0'])) for br in branches)
    c = (num - adj) / denom
    return c

def nodestobranches(Ages_nodes):
    br_dic = [{'t0':-startt,'t1':-endt,'duration':duration,'n':nWGD} for duration,endt,startt,nWGD in Ages_nodes.values()]
    return br_dic

def neg_log_likelihood(theta, branches, target_total):
    a, b = np.exp(theta)
    a = -a
    c = compute_c_from_expectation(a, b, branches, target_total)
    t0_arr = np.array([br['t0'] for br in branches])
    t1_arr = np.array([br['t1'] for br in branches])
    n_arr  = np.array([br['n']  for br in branches])
    Lambda = integral_lambda(a, b, c, t0_arr, t1_arr)
    if np.any(Lambda <= 0):
        return 1e20
    if np.any(lambda_t(a, b, c, t0_arr) <= 0) or np.any(lambda_t(a, b, c, t1_arr) <= 0):
        return 1e20
    ll = np.sum(n_arr * np.log(Lambda) - Lambda - gammaln(n_arr+1))
    return -ll

def fit_nhpp_counts_10000median(branches, target_total):
    a_,b_= [],[]
    for sd in trange(10000):
        bounds = [(np.log(1e-6), np.log(10)),(np.log(1e-6), np.log(10))]
        res = differential_evolution(lambda th: neg_log_likelihood(th, branches, target_total),bounds, maxiter=10000,seed=sd)
        alpha, beta = res.x
        a = -np.exp(alpha)
        b = np.exp(beta)
        a_.append(a)
        b_.append(b) 
    a,b = np.median(a_),np.median(b_)
    c = compute_c_from_expectation(a, b, branches, target_total)
    return a,b,c

def simulation_WGD_dependence(t0, t1, a, b, c, random_seed=None):
    if random_seed is not None: np.random.seed(random_seed)
    max_rate = max(lambda_t(a,b,c,t0),lambda_t(a, b, c, t1))
    events = []
    current_time = t0
    while current_time < t1:
        current_time += np.random.exponential(1 / max_rate)
        if current_time > t1:
            break
        if np.random.random() < (lambda_t(a,b,c,current_time) / max_rate):
            events.append(current_time)
    return events

def get_nodes_ages(tree):
    nodes = [i for i in tree.get_nonterminals()]
    Ages_nodes = {}
    for clade in nodes:
        if type(clade.branch_length) is not float:
            continue
        older_bound = max(clade.depths().values())
        younger_bound = older_bound - clade.branch_length
        Ages_nodes[clade.name] = (clade.branch_length*100,younger_bound*100,older_bound*100)
    tips = [i for i in tree.get_terminals()]
    for clade in tips:
        younger_bound = 0
        older_bound = clade.branch_length
        Ages_nodes[clade.name] = (clade.branch_length*100,younger_bound*100,older_bound*100)
    return Ages_nodes

def getwgdlocations(tree):
    WGD_clade1 = {"WGD1_MALE":tree.common_ancestor({'name':'Pyrus_bretschneideri'},{'name':'Crataegus_pinnatifida'}),"WGD2_PAPI":tree.common_ancestor({'name':'Ammopiptanthus_nanus'},{'name':'Vigna_radiata'}),"WGD2_CAES":tree.common_ancestor({'name':'Acacia_pycnantha'},{'name':'Senna_tora'}),"WGD3_GLYC":tree.common_ancestor({'name':'Glycine_soja'},{'name':'Glycine_max'}),"WGD4_JUGL":tree.common_ancestor({'name':'Carya_illinoinensis'},{'name':'Juglans_regia'}),"WGD5_BRCE":tree.common_ancestor({'name':'Sinapis_alba'},{'name':'Brassica_oleracea'}),"WGD6_BRAS_α":tree.common_ancestor({'name':'Aethionema_arabicum'},{'name':'Brassica_oleracea'}),"WGD7_BRAS_β":tree.common_ancestor({'name':'Capparis_spinosa'},{'name':'Brassica_oleracea'}),"WGD9_TARE":tree.common_ancestor({'name':'Tarenaya_hassleriana'},{'name':'Tarenaya_hassleriana'}),"WGD12_MALV":tree.common_ancestor({'name':'Gossypium_barbadense'},{'name':'Durio_zibethinus'}),"WGD13_MAHE":tree.common_ancestor({'name':'Hevea_brasiliensis'},{'name':'Manihot_esculenta'}),"WGD14_POSA":tree.common_ancestor({'name':'Populus_wilsonii'},{'name':'Salix_arbutifolia'}),"WGD15_LINU_α":tree.common_ancestor({'name':'Linum_usitatissimum'},{'name':'Linum_usitatissimum'}),"WGD15_LINU_β":tree.common_ancestor({'name':'Linum_usitatissimum'},{'name':'Linum_usitatissimum'}),"WGD16_MYRT":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Sonneratia_caseolaris'}),"WGD17_LYTH":tree.common_ancestor({'name':'Sonneratia_caseolaris'},{'name':'Lagerstroemia_indica'}),"WGD19_TRAP":tree.common_ancestor({'name':'Trapa_bicornis'},{'name':'Trapa_bicornis'}),"WGD20_SOLA":tree.common_ancestor({'name':'Capsicum_chinense'},{'name':'Petunia_axillaris'}),"WGD21_OLEA_α":tree.common_ancestor({'name':'Syringa_oblata'},{'name':'Olea_europaea'}),"WGD22_OLEA_β":tree.common_ancestor({'name':'Olea_europaea'},{'name':'Jasminum_sambac'}),"WGD23_ASTE":tree.common_ancestor({'name':'Carthamus_tinctorius'},{'name':'Helianthus_annuus'}),"WGD24_APIO_α":tree.common_ancestor({'name':'Daucus_carota'},{'name':'Coriandrum_sativum'}),"WGD25_APIO_β":tree.common_ancestor({'name':'Panax_ginseng'},{'name':'Daucus_carota'}),"WGD27_AMAR":tree.common_ancestor({'name':'Amaranthus_hypochondriacus'},{'name':'Amaranthus_tuberculatus'}),"WGD28_KALA":tree.common_ancestor({'name':'Kalanchoe_fedtschenkoi'},{'name':'Kalanchoe_fedtschenkoi'}),"WGD29_RANU":tree.common_ancestor({'name':'Aquilegia_oxysepala'},{'name':'Corydalis_tomentella'}),"WGD30_POAC":tree.common_ancestor({'name':'Triticum_aestivum'},{'name':'Streptochaeta_angustifolia'}),"WGD31_ZEAM":tree.common_ancestor({'name':'Zea_mays'},{'name':'Zea_mays'}),"WGD32_POAL":tree.common_ancestor({'name':'Triticum_aestivum'},{'name':'Ananas_comosus'}),"WGD33_MUSA_α":tree.common_ancestor({'name':'Ensete_ventricosum'},{'name':'Musa_acuminata'}),"WGD34_MUSA_β":tree.common_ancestor({'name':'Musa_acuminata'},{'name':'Zingiber_officinale'})}
    WGD_clade2 = {"WGD35_AREC":tree.common_ancestor({'name':'Elaeis_guineensis'},{'name':'Calamus_simplicifolius'}),"WGD38_ASPA":tree.common_ancestor({'name':'Asparagus_setaceus'},{'name':'Asparagus_officinalis'}),"WGD40_ARAC_α":tree.common_ancestor({'name':'Zantedeschia_elliottiana'},{'name':'Lemna_minuta'}),"WGD41_ARAC_β":tree.common_ancestor({'name':'Zantedeschia_elliottiana'},{'name':'Lemna_minuta'}),"WGD42_LAUR_α":tree.common_ancestor({'name':'Phoebe_bournei'},{'name':'Cinnamomum_kanehirae'}),"WGD43_LAUR_β":tree.common_ancestor({'name':'Phoebe_bournei'},{'name':'Magnolia_biondii'}),"WGD44_NYMP":tree.common_ancestor({'name':'Nymphaea_colorata'},{'name':'Brasenia_schreberi'}),"WGD44_BRSH_α":tree.common_ancestor({'name':'Brasenia_schreberi'},{'name':'Brasenia_schreberi'}),"WGD44_BRSH_β":tree.common_ancestor({'name':'Brasenia_schreberi'},{'name':'Brasenia_schreberi'}),"WGD45_ACOR":tree.common_ancestor({'name':'Acorus_americanus'},{'name':'Acorus_tatarinowii'}),"WGD46_DIOS":tree.common_ancestor({'name':'Dioscorea_alata'},{'name':'Trichopus_zeylanicus'}),"WGD47_TRIC":tree.common_ancestor({'name':'Trichopus_zeylanicus'},{'name':'Trichopus_zeylanicus'}),"WGD48_ACAN":tree.common_ancestor({'name':'Acanthochlamys_bracteata'},{'name':'Acanthochlamys_bracteata'}),"WGD50_CUCU":tree.common_ancestor({'name':'Cucumis_melo'},{'name':'Datisca_glomerata'}),"WGD51_CUBI":tree.common_ancestor({'name':'Cucurbita_maxima'},{'name':'Cucurbita_argyrosperma'}),"WGD52_BEGO":tree.common_ancestor({'name':'Begonia_fuchsioides'},{'name':'Begonia_darthvaderiana'}),"WGD52_BEFU":tree.common_ancestor({'name':'Begonia_fuchsioides'},{'name':'Begonia_fuchsioides'}),"WGD53_LUPI":tree.common_ancestor({'name':'Lupinus_angustifolius'},{'name':'Lupinus_angustifolius'}),"WGD54_CONV":tree.common_ancestor({'name':'Ipomoea_nil'},{'name':'Cuscuta_australis'}),"WGD57_ACTI_α":tree.common_ancestor({'name':'Actinidia_chinensis'},{'name':'Actinidia_eriantha'}),"WGD58_ACTI_β":tree.common_ancestor({'name':'Rhododendron_simsii'},{'name':'Gilia_yorkii'}),"WGD59_NELU":tree.common_ancestor({'name':'Nelumbo_nucifera'},{'name':'Nelumbo_nucifera'}),"WGD60_PASE_β":tree.common_ancestor({'name':'Papaver_somniferum'},{'name':'Papaver_setigerum'}),"WGD61_PASE_α":tree.common_ancestor({'name':'Papaver_setigerum'},{'name':'Papaver_setigerum'}),"WGD63_ORCH":tree.common_ancestor({'name':'Cremastra_appendiculata'},{'name':'Apostasia_shenzhenica'}),"WGD64_ZANT":tree.common_ancestor({'name':'Zanthoxylum_bungeanum'},{'name':'Zanthoxylum_armatum'}),"WGD65_EURY":tree.common_ancestor({'name':'Euryale_ferox'},{'name':'Euryale_ferox'}),"WGD66_PIPE":tree.common_ancestor({'name':'Piper_nigrum'},{'name':'Piper_nigrum'}),"WGD67_CALY":tree.common_ancestor({'name':'Chimonanthus_salicifolius'},{'name':'Chimonanthus_salicifolius'})}
    WGD_clade3 = {"WGD68_CHLO":tree.common_ancestor({'name':'Chloranthus_spicatus'},{'name':'Chloranthus_sessilifolius'}),"WGD68_LARD":tree.common_ancestor({'name':'Akebia_trifoliata'},{'name':'Akebia_trifoliata'}),"WGD69_BUXA":tree.common_ancestor({'name':'Buxus_sinica'},{'name':'Buxus_austroyunnanensis'}),"WGD70_TROC_α":tree.common_ancestor({'name':'Trochodendron_aralioides'},{'name':'Tetracentron_sinense'}),"WGD70_TROC_β":tree.common_ancestor({'name':'Trochodendron_aralioides'},{'name':'Tetracentron_sinense'}),"WGD71_ACHN":tree.common_ancestor({'name':'Achnatherum_splendens'},{'name':'Achnatherum_splendens'}),"WGD72_BONI":tree.common_ancestor({'name':'Bonia_amplexicaulis'},{'name':'Bonia_amplexicaulis'}),"WGD76_PHRA":tree.common_ancestor({'name':'Phragmites_australis'},{'name':'Phragmites_australis'}),"WGD77_ZIZA":tree.common_ancestor({'name':'Zizania_latifolia'},{'name':'Zizania_palustris'}),"WGD78_ZING":tree.common_ancestor({'name':'Zingiber_officinale'},{'name':'Lanxangia_tsaoko'}),"WGD79_POLG":tree.common_ancestor({'name':'Fagopyrum_tataricum'},{'name':'Rheum_tanguticum'}),"WGD80_CARY":tree.common_ancestor({'name':'Gypsophila_paniculata'},{'name':'Silene_latifolia'}),"WGD81_ACPT":tree.common_ancestor({'name':'Portulaca_amilis'},{'name':'Hylocereus_undatus'}),"WGD82_NYSS":tree.common_ancestor({'name':'Davidia_involucrata'},{'name':'Nyssa_yunnanensis'}),"WGD83_PRIM":tree.common_ancestor({'name':'Aegiceras_corniculatum'},{'name':'Primula_veris'}),"WGD84_UTRI":tree.common_ancestor({'name':'Utricularia_gibba'},{'name':'Utricularia_gibba'}),"WGD85_OROB":tree.common_ancestor({'name':'Origanum_majorana'},{'name':'Buddleja_alternifolia'}),"WGD86_STRI":tree.common_ancestor({'name':'Striga_hermonthica'},{'name':'Striga_asiatica'}),"WGD87_STRE":tree.common_ancestor({'name':'Streptocarpus_rexii'},{'name':'Streptocarpus_rexii'}),"WGD88_ANTI":tree.common_ancestor({'name':'Antirrhinum_majus'},{'name':'Antirrhinum_majus'}),"WGD89_SALV":tree.common_ancestor({'name':'Salvia_splendens'},{'name':'Salvia_splendens'}),"WGD91_BORA_α":tree.common_ancestor({'name':'Lithospermum_erythrorhizon'},{'name':'Lithospermum_erythrorhizon'}),"WGD92_BORA_β":tree.common_ancestor({'name':'Lithospermum_erythrorhizon'},{'name':'Lithospermum_erythrorhizon'}),"WGD93_SMAL":tree.common_ancestor({'name':'Smallanthus_sonchifolius'},{'name':'Smallanthus_sonchifolius'}),"WGD95_HELI":tree.common_ancestor({'name':'Mikania_micrantha'},{'name':'Helianthus_annuus'}),"WGD97_ARAL":tree.common_ancestor({'name':'Eleutherococcus_senticosus'},{'name':'Panax_ginseng'}),"WGD98_ELEU":tree.common_ancestor({'name':'Eleutherococcus_senticosus'},{'name':'Eleutherococcus_senticosus'}),"WGD99_DIPS":tree.common_ancestor({'name':'Lonicera_japonica'},{'name':'Lonicera_japonica'})}
    WGD_clade4 = {"WGD100_AQUI":tree.common_ancestor({'name':'Ilex_polyneura'},{'name':'Ilex_polyneura'}),"WGD101_BAUH":tree.common_ancestor({'name':'Bauhinia_variegata'},{'name':'Bauhinia_variegata'}),"WGD102_HIPP_α":tree.common_ancestor({'name':'Hippophae_rhamnoides'},{'name':'Elaeagnus_angustifolia'}),"WGD102_HIPP_β":tree.common_ancestor({'name':'Hippophae_rhamnoides'},{'name':'Elaeagnus_angustifolia'}),"WGD103_TRIP":tree.common_ancestor({'name':'Tripterygium_wilfordii'},{'name':'Tripterygium_wilfordii'}),"WGD104_PASS":tree.common_ancestor({'name':'Passiflora_edulis'},{'name':'Passiflora_organensis'}),"WGD105_RHIZ":tree.common_ancestor({'name':'Kandelia_candel'},{'name':'Bruguiera_sexangula'}),"WGD106_MELA_α":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Melastoma_dodecandrum'}),"WGD107_MELA_β":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Melastoma_dodecandrum'}),"WGD108_HIBI":tree.common_ancestor({'name':'Hibiscus_cannabinus'},{'name':'Abelmoschus_esculentus'}),"WGD109_DIPT":tree.common_ancestor({'name':'Shorea_leprosula'},{'name':'Shorea_leprosula'}),"WGD110_THYM":tree.common_ancestor({'name':'Aquilaria_sinensis'},{'name':'Stellera_chamaejasme'}),"WGD111_STAP":tree.common_ancestor({'name':'Euscaphis_japonica'},{'name':'Euscaphis_japonica'}),"WGD112_BRET":tree.common_ancestor({'name':'Bretschneidera_sinensis'},{'name':'Bretschneidera_sinensis'}),"WGD113_CAPP":tree.common_ancestor({'name':'Capparis_spinosa'},{'name':'Capparis_spinosa'}),"WGD115_LOME":tree.common_ancestor({'name':'Lobularia_maritima'},{'name':'Megacarpaea_delavayi'}),"WGD116_ORYC":tree.common_ancestor({'name':'Orychophragmus_violaceus'},{'name':'Orychophragmus_violaceus'}),"WGD117_TOON":tree.common_ancestor({'name':'Toona_sinensis'},{'name':'Toona_sinensis'}),"WGD119_MANG":tree.common_ancestor({'name':'Mangifera_indica'},{'name':'Mangifera_indica'}),"WGD120_PROT":tree.common_ancestor({'name':'Macadamia_jansenii'},{'name':'Protea_cynaroides'}),"WGD122_CYMO":tree.common_ancestor({'name':'Cymodocea_nodosa'},{'name':'Cymodocea_nodosa'}),"WGD123_POZO":tree.common_ancestor({'name':'Potamogeton_acutifolius'},{'name':'Zostera_muelleri'}),"WGD124_SEAG":tree.common_ancestor({'name':'Thalassia_testudinum'},{'name':'Zostera_muelleri'}),"WGD125_CERA_α":tree.common_ancestor({'name':'Ceratophyllum_demersum'},{'name':'Ceratophyllum_demersum'}),"WGD126_CERA_β":tree.common_ancestor({'name':'Ceratophyllum_demersum'},{'name':'Ceratophyllum_demersum'}),"WGD126b_TAMA":tree.common_ancestor({'name':'Tamarix_chinensis'},{'name':'Tamarix_chinensis'}),"WGD126c_SAPI":tree.common_ancestor({'name':'Saururus_chinensis'},{'name':'Piper_nigrum'}),"WGD126d_SANT":tree.common_ancestor({'name':'Santalum_yasi'},{'name':'Santalum_yasi'}),"WGD127_LEPI":tree.common_ancestor({'name':'Lepidium_sativum'},{'name':'Lepidium_sativum'}),"WGD128_ACON":tree.common_ancestor({'name':'Aconitum_vilmorinianum'},{'name':'Aconitum_vilmorinianum'}),"WGD131_POLY_α":tree.common_ancestor({'name':'Polygala_tenuifolia'},{'name':'Polygala_tenuifolia'}),"WGD132_POLY_β":tree.common_ancestor({'name':'Polygala_tenuifolia'},{'name':'Polygala_tenuifolia'})}
    WGD_clade5 = {"WGD133_MESU":tree.common_ancestor({'name':'Mesua_ferrea'},{'name':'Mesua_ferrea'}),"WGD135_CORN":tree.common_ancestor({'name':'Cornus_wilsoniana'},{'name':'Cornus_wilsoniana'}),"WGD136_AESC":tree.common_ancestor({'name':'Aesculus_chinensis'},{'name':'Aesculus_chinensis'}),"WGD137_ABEL":tree.common_ancestor({'name':'Abelmoschus_esculentus'},{'name':'Abelmoschus_esculentus'}),"WGD138_RHSE":tree.common_ancestor({'name':'Rhodiola_crenulata'},{'name':'Sedum_album'}),"WGD141_ROSM":tree.common_ancestor({'name':'Rosmarinus_officinalis'},{'name':'Rosmarinus_officinalis'}),"WGD142_OCIM":tree.common_ancestor({'name':'Ocimum_basilicum'},{'name':'Ocimum_basilicum'}),"WGD143_NEOL":tree.common_ancestor({'name':'Neolamarckia_cadamba'},{'name':'Neolamarckia_cadamba'}),"WGD144_COEU":tree.common_ancestor({'name':'Tamarix_chinensis'},{'name':'Vigna_radiata'}),"WGD145_COMO":tree.common_ancestor({'name':'Acanthochlamys_bracteata'},{'name':'Triticum_aestivum'}),"WGD146_BOTH":tree.common_ancestor({'name':'Bothriochloa_decipiens'},{'name':'Bothriochloa_decipiens'}),"WGD147_ZOYS":tree.common_ancestor({'name':'Zoysia_matrella'},{'name':'Zoysia_japonica'})}
    return {**WGD_clade1,**WGD_clade2,**WGD_clade3,**WGD_clade4,**WGD_clade5}

def getTree_nodes_WGDIDs(WGDIDs_treenodeids):
    Tree_nodes_WGDids = {}
    for wgdid,nodeid in WGDIDs_treenodeids.items():
        if nodeid in Tree_nodes_WGDids:
            Tree_nodes_WGDids[nodeid] += [wgdid]
        else:
            Tree_nodes_WGDids[nodeid] = [wgdid]
    return Tree_nodes_WGDids

if __name__ == '__main__':
        cli()
