# Usage
# python simulatewgd.py simulatewgdontree Maintree.txt WGDdate.tsv

from io import StringIO
from Bio import Phylo
import click
from rich.logging import RichHandler
import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy import stats
from scipy.stats import chi2
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import ks_2samp
from sklearn import mixture
from sklearn.metrics import silhouette_score
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_rgb
import colorsys
import math
import itertools

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    This script is to simulate WGD on a timetree
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

def addvvline(ax,xvalue,color,lstyle,labell):
    if labell == '': ax.axvline(xvalue,color=color, ls=lstyle, lw=1)
    else: ax.axvline(xvalue,color=color, ls=lstyle, lw=1, label='{}: {:.1f}'.format(labell,xvalue))
    return ax

def get_totalH(Hs):
    CHF = 0
    for i in Hs: CHF = CHF + i
    return CHF

def kde_mode(kde_x, kde_y):
    maxy_iloc = np.argmax(kde_y)
    mode = kde_x[maxy_iloc]
    return mode, max(kde_y)

def find_first_parent(tree,target):
    All_Parents = [(node,node.distance(target)) for node in tree.get_nonterminals() if node.is_parent_of(target) and node != target]
    return sorted(All_Parents,key = lambda x:x[1])[0][0]

def adjust_saturation(color_name, saturation_factor):
    rgb = to_rgb(color_name)
    r, g, b = rgb
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    s = max(0, min(1, s * saturation_factor))
    return colorsys.hls_to_rgb(h, l, s)

def poisson_pmf(k, lambda_):
    return (lambda_**k) * np.exp(-lambda_) / math.factorial(k)

def modifiedF_c():
    Family_clades2 = {}
    for key,value in Family_clades.items():
        if value in ["Amborellales","Nymphaeales"]:
            Family_clades2[key] = "ANA grade"
        else:
            Family_clades2[key] = value
    return Family_clades2

taxonomy_prset={'Torenia_fournieri':{'family': 'Linderniaceae', 'order': 'Lamiales'},'Platanus_x_acerifolia':{'family': 'Platanaceae', 'order': 'Proteales'},'Pontederia_crassipes':{'family': 'Pontederiaceae', 'order': 'Commelinales'},'Abelmoschus_esculentus': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Abrus_precatorius': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Acacia_pycnantha': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Acanthochlamys_bracteata': {'family': 'Velloziaceae', 'order': 'Pandanales'}, 'Acer_pseudosieboldianum': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Acer_yangbiense': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Achnatherum_splendens': {'family': 'Poaceae', 'order': 'Poales'}, 'Aconitum_vilmorinianum': {'family': 'Ranunculaceae', 'order': 'Ranunculales'}, 'Acorus_americanus': {'family': 'Acoraceae', 'order': 'Acorales'}, 'Acorus_tatarinowii': {'family': 'Acoraceae', 'order': 'Acorales'}, 'Actinidia_chinensis': {'family': 'Actinidiaceae', 'order': 'Ericales'}, 'Actinidia_eriantha': {'family': 'Actinidiaceae', 'order': 'Ericales'}, 'Adenosma_buchneroides': {'family': 'Plantaginaceae', 'order': 'Lamiales'}, 'Aegiceras_corniculatum': {'family': 'Primulaceae', 'order': 'Ericales'}, 'Aegilops_tauschii': {'family': 'Poaceae', 'order': 'Poales'}, 'Aeschynomene_evenia': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Aesculus_chinensis': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Aethionema_arabicum': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Akebia_trifoliata': {'family': 'Lardizabalaceae', 'order': 'Ranunculales'}, 'Allium_cepa': {'family': 'Amaryllidaceae', 'order': 'Asparagales'}, 'Allium_sativum': {'family': 'Amaryllidaceae', 'order': 'Asparagales'}, 'Amaranthus_cruentus': {'family': 'Amaranthaceae', 'order': 'Caryophyllales'}, 'Amaranthus_hybridus': {'family': 'Amaranthaceae', 'order': 'Caryophyllales'}, 'Amaranthus_hypochondriacus': {'family': 'Amaranthaceae', 'order': 'Caryophyllales'}, 'Amaranthus_palmeri': {'family': 'Amaranthaceae', 'order': 'Caryophyllales'}, 'Amaranthus_tuberculatus': {'family': 'Amaranthaceae', 'order': 'Caryophyllales'}, 'Amborella_trichopoda': {'family': 'Amborellaceae', 'order': 'Amborellales'}, 'Ammopiptanthus_nanus': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Amorphophallus_konjac': {'family': 'Araceae', 'order': 'Alismatales'}, 'Amphicarpaea_edgeworthii': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Ananas_comosus': {'family': 'Bromeliaceae', 'order': 'Poales'}, 'Andrographis_paniculata': {'family': 'Acanthaceae', 'order': 'Lamiales'}, 'Angelica_sinensis': {'family': 'Apiaceae', 'order': 'Apiales'}, 'Antirrhinum_majus': {'family': 'Plantaginaceae', 'order': 'Lamiales'}, 'Apium_graveolens': {'family': 'Apiaceae', 'order': 'Apiales'}, 'Apostasia_shenzhenica': {'family': 'Orchidaceae', 'order': 'Asparagales'}, 'Aquilaria_sinensis': {'family': 'Thymelaeaceae', 'order': 'Malvales'}, 'Aquilegia_coerulea': {'family': 'Ranunculaceae', 'order': 'Ranunculales'}, 'Aquilegia_oxysepala': {'family': 'Ranunculaceae', 'order': 'Ranunculales'}, 'Arabidopsis_lyrata': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Arabidopsis_thaliana': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Arabis_alpina': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Arachis_duranensis': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Arachis_hypogaea': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Aralia_elata': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Arctium_lappa': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Areca_catechu': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Argentina_anserina': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Aristolochia_contorta': {'family': 'Aristolochiaceae', 'order': 'Piperales'}, 'Aristolochia_fimbriata': {'family': 'Aristolochiaceae', 'order': 'Piperales'}, 'Artemisia_annua': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Artocarpus_nanchuanensis': {'family': 'Moraceae', 'order': 'Rosales'}, 'Asparagus_officinalis': {'family': 'Asparagaceae', 'order': 'Asparagales'}, 'Asparagus_setaceus': {'family': 'Asparagaceae', 'order': 'Asparagales'}, 'Atriplex_hortensis': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Avena_insularis': {'family': 'Poaceae', 'order': 'Poales'}, 'Avena_sativa': {'family': 'Poaceae', 'order': 'Poales'}, 'Averrhoa_carambola': {'family': 'Oxalidaceae', 'order': 'Oxalidales'}, 'Avicennia_marina': {'family': 'Acanthaceae', 'order': 'Lamiales'}, 'Bauhinia_variegata': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Begonia_darthvaderiana': {'family': 'Begoniaceae', 'order': 'Cucurbitales'}, 'Begonia_fuchsioides': {'family': 'Begoniaceae', 'order': 'Cucurbitales'}, 'Benincasa_hispida': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Beta_patula': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Beta_vulgaris': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Betula_pendula': {'family': 'Betulaceae', 'order': 'Fagales'}, 'Boechera_stricta': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Boehmeria_nivea': {'family': 'Urticaceae', 'order': 'Rosales'}, 'Bombax_ceiba': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Bonia_amplexicaulis': {'family': 'Poaceae', 'order': 'Poales'}, 'Bothriochloa_decipiens': {'family': 'Poaceae', 'order': 'Poales'}, 'Brachypodium_distachyon': {'family': 'Poaceae', 'order': 'Poales'}, 'Brachypodium_hybridum': {'family': 'Poaceae', 'order': 'Poales'}, 'Brassica_napus': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Brassica_oleracea': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Brassica_rapa': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Bretschneidera_sinensis': {'family': 'Akaniaceae', 'order': 'Brassicales'}, 'Bruguiera_sexangula': {'family': 'Rhizophoraceae', 'order': 'Malpighiales'}, 'Buddleja_alternifolia': {'family': 'Scrophulariaceae', 'order': 'Lamiales'}, 'Buxus_austroyunnanensis': {'family': 'Buxaceae', 'order': 'Buxales'}, 'Buxus_sinica': {'family': 'Buxaceae', 'order': 'Buxales'}, 'Cajanus_cajan': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Calamus_simplicifolius': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Callicarpa_americana': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Calotropis_gigantea': {'family': 'Apocynaceae', 'order': 'Gentianales'}, 'Camelina_sativa': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Camellia_chekiangoleosa': {'family': 'Theaceae', 'order': 'Ericales'}, 'Camellia_oleifera': {'family': 'Theaceae', 'order': 'Ericales'}, 'Camptotheca_acuminata': {'family': 'Nyssaceae', 'order': 'Cornales'}, 'Cannabis_sativa': {'family': 'Cannabaceae', 'order': 'Rosales'}, 'Capparis_spinosa': {'family': 'Capparaceae', 'order': 'Brassicales'}, 'Capsella_rubella': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Capsicum_annuum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Capsicum_baccatum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Capsicum_chinense': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Cardamine_hirsuta': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Carex_cristatella': {'family': 'Cyperaceae', 'order': 'Poales'}, 'Carex_scoparia': {'family': 'Cyperaceae', 'order': 'Poales'}, 'Carica_papaya': {'family': 'Caricaceae', 'order': 'Brassicales'}, 'Carpinus_fangiana': {'family': 'Betulaceae', 'order': 'Fagales'}, 'Carthamus_tinctorius': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Carya_illinoinensis': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Castanea_mollissima': {'family': 'Fagaceae', 'order': 'Fagales'}, 'Casuarina_glauca': {'family': 'Casuarinaceae', 'order': 'Fagales'}, 'Catharanthus_roseus': {'family': 'Apocynaceae', 'order': 'Gentianales'}, 'Cephalotus_follicularis': {'family': 'Cephalotaceae', 'order': 'Oxalidales'}, 'Ceratophyllum_demersum': {'family': 'Ceratophyllaceae', 'order': 'Ceratophyllales'}, 'Cercis_chinensis': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Chenopodium_quinoa': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Chimonanthus_salicifolius': {'family': 'Calycanthaceae', 'order': 'Laurales'}, 'Chiococca_alba': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Chloranthus_sessilifolius': {'family': 'Chloranthaceae', 'order': 'Chloranthales'}, 'Chloranthus_spicatus': {'family': 'Chloranthaceae', 'order': 'Chloranthales'}, 'Chosenia_arbutifolia': {'family': 'Salicaceae', 'order': 'Malpighiales'}, 'Chrysanthemum_seticuspe': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Cicer_arietinum': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Cichorium_endivia': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Cichorium_intybus': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Cinchona_pubescens': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Cinnamomum_kanehirae': {'family': 'Lauraceae', 'order': 'Laurales'}, 'Citrullus_lanatus': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Citrus_clementina': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Citrus_maxima': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Citrus_sinensis': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Cleome_violacea': {'family': 'Cleomaceae', 'order': 'Brassicales'}, 'Clethra_arborea': {'family': 'Clethraceae', 'order': 'Ericales'}, 'Cocos_nucifera': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Coffea_canephora': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Coffea_humblotiana': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Coptis_chinensis': {'family': 'Ranunculaceae', 'order': 'Ranunculales'}, 'Corchorus_capsularis': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Coriandrum_sativum': {'family': 'Apiaceae', 'order': 'Apiales'}, 'Cornus_wilsoniana': {'family': 'Cornaceae', 'order': 'Cornales'}, 'Corydalis_tomentella': {'family': 'Papaveraceae', 'order': 'Ranunculales'}, 'Corylus_heterophylla': {'family': 'Betulaceae', 'order': 'Fagales'}, 'Crataegus_pinnatifida': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Cremastra_appendiculata': {'family': 'Orchidaceae', 'order': 'Asparagales'}, 'Cucumis_melo': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Cucumis_sativus': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Cucurbita_argyrosperma': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Cucurbita_maxima': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Cuscuta_australis': {'family': 'Convolvulaceae', 'order': 'Solanales'}, 'Cuscuta_campestris': {'family': 'Convolvulaceae', 'order': 'Solanales'}, 'Cycas_panzhihuaensis': {'family': 'Cycadaceae', 'order': 'Cycadales'}, 'Cyclocarya_paliurus': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Cymodocea_nodosa': {'family': 'Cymodoceaceae', 'order': 'Alismatales'}, 'Cynara_cardunculus': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Daemonorops_jenkinsiana': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Dalbergia_odorifera': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Datisca_glomerata': {'family': 'Datiscaceae', 'order': 'Cucurbitales'}, 'Datura_stramonium': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Daucus_carota': {'family': 'Apiaceae', 'order': 'Apiales'}, 'Davidia_involucrata': {'family': 'Nyssaceae', 'order': 'Cornales'}, 'Dendrobium_nobile': {'family': 'Orchidaceae', 'order': 'Asparagales'}, 'Dianthus_caryophyllus': {'family': 'Caryophyllaceae', 'order': 'Caryophyllales'}, 'Dichanthelium_oligosanthes': {'family': 'Poaceae', 'order': 'Poales'}, 'Dimocarpus_longan': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Dioscorea_alata': {'family': 'Dioscoreaceae', 'order': 'Dioscoreales'}, 'Dioscorea_rotundata': {'family': 'Dioscoreaceae', 'order': 'Dioscoreales'}, 'Diospyros_lotus': {'family': 'Ebenaceae', 'order': 'Ericales'}, 'Draba_nivalis': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Dryas_drummondii': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Durio_zibethinus': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Echinochloa_crus-galli': {'family': 'Poaceae', 'order': 'Poales'}, 'Elaeagnus_angustifolia': {'family': 'Elaeagnaceae', 'order': 'Rosales'}, 'Elaeis_guineensis': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Eleutherococcus_senticosus': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Ensete_glaucum': {'family': 'Musaceae', 'order': 'Zingiberales'}, 'Ensete_ventricosum': {'family': 'Musaceae', 'order': 'Zingiberales'}, 'Entada_phaseoloides': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Epimedium_pubescens': {'family': 'Berberidaceae', 'order': 'Ranunculales'}, 'Eragrostis_nindensis': {'family': 'Poaceae', 'order': 'Poales'}, 'Eragrostis_tef': {'family': 'Poaceae', 'order': 'Poales'}, 'Erigeron_canadensis': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Eriobotrya_japonica': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Eruca_sativa': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Erysimum_cheiranthoides': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Eucalyptus_grandis': {'family': 'Myrtaceae', 'order': 'Myrtales'}, 'Euphorbia_lathyris': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Euryale_ferox': {'family': 'Nymphaeaceae', 'order': 'Nymphaeales'}, 'Euscaphis_japonica': {'family': 'Staphyleaceae', 'order': 'Crossosomatales'}, 'Eutrema_salsugineum': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Fagopyrum_tataricum': {'family': 'Polygonaceae', 'order': 'Caryophyllales'}, 'Fagus_sylvatica': {'family': 'Fagaceae', 'order': 'Fagales'}, 'Faidherbia_albida': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Ficus_hispida': {'family': 'Moraceae', 'order': 'Rosales'}, 'Ficus_microcarpa': {'family': 'Moraceae', 'order': 'Rosales'}, 'Fragaria_pentaphylla': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Fragaria_vesca': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Fraxinus_pennsylvanica': {'family': 'Oleaceae', 'order': 'Lamiales'}, 'Gardenia_jasminoides': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Genlisea_aurea': {'family': 'Lentibulariaceae', 'order': 'Lamiales'}, 'Gilia_yorkii': {'family': 'Polemoniaceae', 'order': 'Ericales'}, 'Glebionis_coronaria': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Glycine_max': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Glycine_soja': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Gossypium_barbadense': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Gossypium_hirsutum': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Gossypium_raimondii': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Guadua_angustifolia': {'family': 'Poaceae', 'order': 'Poales'}, 'Gypsophila_paniculata': {'family': 'Caryophyllaceae', 'order': 'Caryophyllales'}, 'Haloxylon_ammodendron': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Hamamelis_virginiana': {'family': 'Hamamelidaceae', 'order': 'Saxifragales'}, 'Handroanthus_impetiginosus': {'family': 'Bignoniaceae', 'order': 'Lamiales'}, 'Helianthus_annuus': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Heritiera_littoralis': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Hevea_brasiliensis': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Hibiscus_cannabinus': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Hippophae_rhamnoides': {'family': 'Elaeagnaceae', 'order': 'Rosales'}, 'Hirschfeldia_incana': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Hordeum_vulgare': {'family': 'Poaceae', 'order': 'Poales'}, 'Humulus_lupulus': {'family': 'Cannabaceae', 'order': 'Rosales'}, 'Hydrangea_macrophylla': {'family': 'Hydrangeaceae', 'order': 'Cornales'}, 'Hylocereus_undatus': {'family': 'Cactaceae', 'order': 'Caryophyllales'}, 'Ilex_polyneura': {'family': 'Aquifoliaceae', 'order': 'Aquifoliales'}, 'Ipomoea_nil': {'family': 'Convolvulaceae', 'order': 'Solanales'}, 'Isatis_indigotica': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Jacaranda_mimosifolia': {'family': 'Bignoniaceae', 'order': 'Lamiales'}, 'Jasminum_sambac': {'family': 'Oleaceae', 'order': 'Lamiales'}, 'Jatropha_curcas': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Juglans_mandshurica': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Juglans_regia': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Juglans_sigillata': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Juncus_effusus': {'family': 'Juncaceae', 'order': 'Poales'}, 'Juncus_inflexus': {'family': 'Juncaceae', 'order': 'Poales'}, 'Kalanchoe_fedtschenkoi': {'family': 'Crassulaceae', 'order': 'Saxifragales'}, 'Kandelia_candel': {'family': 'Rhizophoraceae', 'order': 'Malpighiales'}, 'Kandelia_obovata': {'family': 'Rhizophoraceae', 'order': 'Malpighiales'}, 'Kingdonia_uniflora': {'family': 'Circaeasteraceae', 'order': 'Ranunculales'}, 'Lablab_purpureus': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Lactuca_sativa': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Lagerstroemia_indica': {'family': 'Lythraceae', 'order': 'Myrtales'}, 'Lanxangia_tsaoko': {'family': 'Zingiberaceae', 'order': 'Zingiberales'}, 'Leersia_perrieri': {'family': 'Poaceae', 'order': 'Poales'}, 'Lemna_minor': {'family': 'Araceae', 'order': 'Alismatales'}, 'Lemna_minuta': {'family': 'Araceae', 'order': 'Alismatales'}, 'Lens_culinaris': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Lens_ervoides': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Lepidium_sativum': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Lindenbergia_luchunensis': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Linum_usitatissimum': {'family': 'Linaceae', 'order': 'Malpighiales'}, 'Liriodendron_chinense': {'family': 'Magnoliaceae', 'order': 'Magnoliales'}, 'Liriodendron_tulipifera': {'family': 'Magnoliaceae', 'order': 'Magnoliales'}, 'Litchi_chinensis': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Lithospermum_erythrorhizon': {'family': 'Boraginaceae', 'order': 'Boraginales'}, 'Lobularia_maritima': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Iochroma_cyaneum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Lonicera_japonica': {'family': 'Caprifoliaceae', 'order': 'Dipsacales'}, 'Lotus_japonicus': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Lupinus_angustifolius': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Macadamia_integrifolia': {'family': 'Proteaceae', 'order': 'Proteales'}, 'Macadamia_jansenii': {'family': 'Proteaceae', 'order': 'Proteales'}, 'Macleaya_cordata': {'family': 'Papaveraceae', 'order': 'Ranunculales'}, 'Macrotyloma_uniflorum': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Magnolia_biondii': {'family': 'Magnoliaceae', 'order': 'Magnoliales'}, 'Malania_oleifera': {'family': 'Ximeniaceae', 'order': 'Santalales'}, 'Malus_domestica': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Malus_sieversii': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Mangifera_indica': {'family': 'Anacardiaceae', 'order': 'Sapindales'}, 'Manihot_esculenta': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Medicago_sativa': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Medicago_truncatula': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Megacarpaea_delavayi': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Megadenia_pygmaea': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Melaleuca_alternifolia': {'family': 'Myrtaceae', 'order': 'Myrtales'}, 'Melastoma_dodecandrum': {'family': 'Melastomataceae', 'order': 'Myrtales'}, 'Melilotus_albus': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Mesua_ferrea': {'family': 'Calophyllaceae', 'order': 'Malpighiales'}, 'Microthlaspi_erraticum': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Mikania_micrantha': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Mimulus_guttatus': {'family': 'Phrymaceae', 'order': 'Lamiales'}, 'Miscanthus_lutarioriparius': {'family': 'Poaceae', 'order': 'Poales'}, 'Morella_rubra': {'family': 'Myricaceae', 'order': 'Fagales'}, 'Moringa_oleifera': {'family': 'Moringaceae', 'order': 'Brassicales'}, 'Morus_notabilis': {'family': 'Moraceae', 'order': 'Rosales'}, 'Mucuna_pruriens': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Musa_acuminata': {'family': 'Musaceae', 'order': 'Zingiberales'}, 'Musa_balbisiana': {'family': 'Musaceae', 'order': 'Zingiberales'}, 'Nelumbo_nucifera': {'family': 'Nelumbonaceae', 'order': 'Proteales'}, 'Neolamarckia_cadamba': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Nephelium_lappaceum': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Nicotiana_attenuata': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Nicotiana_tabacum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Nicotiana_tomentosiformis': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Nymphaea_colorata': {'family': 'Nymphaeaceae', 'order': 'Nymphaeales'}, 'Nymphaea_thermarum': {'family': 'Nymphaeaceae', 'order': 'Nymphaeales'}, 'Nymphoides_indica': {'family': 'Menyanthaceae', 'order': 'Asterales'}, 'Nyssa_sinensis': {'family': 'Nyssaceae', 'order': 'Cornales'}, 'Nyssa_yunnanensis': {'family': 'Nyssaceae', 'order': 'Cornales'}, 'Ochroma_pyramidale': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Ocimum_basilicum': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Olea_europaea': {'family': 'Oleaceae', 'order': 'Lamiales'}, 'Olyra_latifolia': {'family': 'Poaceae', 'order': 'Poales'}, 'Ophiorrhiza_pumila': {'family': 'Rubiaceae', 'order': 'Gentianales'}, 'Origanum_majorana': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Origanum_vulgare': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Orobanche_cumana': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Oropetium_thomaeum': {'family': 'Poaceae', 'order': 'Poales'}, 'Orychophragmus_violaceus': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Oryza_glaberrima': {'family': 'Poaceae', 'order': 'Poales'}, 'Oryza_punctata': {'family': 'Poaceae', 'order': 'Poales'}, 'Oryza_rufipogon': {'family': 'Poaceae', 'order': 'Poales'}, 'Oryza_sativa': {'family': 'Poaceae', 'order': 'Poales'}, 'Ostrya_rehderiana': {'family': 'Betulaceae', 'order': 'Fagales'}, 'Oxytropis_ochrocephala': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Paeonia_ostii': {'family': 'Paeoniaceae', 'order': 'Saxifragales'}, 'Panax_ginseng': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Panax_japonicus': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Panax_quinquefolius': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Panax_stipuleanatus': {'family': 'Araliaceae', 'order': 'Apiales'}, 'Panicum_virgatum': {'family': 'Poaceae', 'order': 'Poales'}, 'Papaver_rhoeas': {'family': 'Papaveraceae', 'order': 'Ranunculales'}, 'Papaver_setigerum': {'family': 'Papaveraceae', 'order': 'Ranunculales'}, 'Papaver_somniferum': {'family': 'Papaveraceae', 'order': 'Ranunculales'}, 'Parasponia_andersonii': {'family': 'Cannabaceae', 'order': 'Rosales'}, 'Passiflora_edulis': {'family': 'Passifloraceae', 'order': 'Malpighiales'}, 'Passiflora_organensis': {'family': 'Passifloraceae', 'order': 'Malpighiales'}, 'Paulownia_fortunei': {'family': 'Paulowniaceae', 'order': 'Lamiales'}, 'Pennisetum_glaucum': {'family': 'Poaceae', 'order': 'Poales'}, 'Perilla_frutescens': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Persea_americana': {'family': 'Lauraceae', 'order': 'Laurales'}, 'Petunia_axillaris': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Phalaenopsis_equestris': {'family': 'Orchidaceae', 'order': 'Asparagales'}, 'Pharus_latifolius': {'family': 'Poaceae', 'order': 'Poales'}, 'Phaseolus_vulgaris': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Phelipanche_aegyptiaca': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Phoebe_bournei': {'family': 'Lauraceae', 'order': 'Laurales'}, 'Phoenix_dactylifera': {'family': 'Arecaceae', 'order': 'Arecales'}, 'Phragmites_australis': {'family': 'Poaceae', 'order': 'Poales'}, 'Phtheirospermum_japonicum': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Physalis_grisea': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Physalis_pruinosa': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Picrorhiza_kurrooa': {'family': 'Plantaginaceae', 'order': 'Lamiales'}, 'Piper_nigrum': {'family': 'Piperaceae', 'order': 'Piperales'}, 'Pistacia_vera': {'family': 'Anacardiaceae', 'order': 'Sapindales'}, 'Pistia_stratiotes': {'family': 'Araceae', 'order': 'Alismatales'}, 'Pisum_sativm': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Platycarya_strobilacea': {'family': 'Juglandaceae', 'order': 'Fagales'}, 'Pogostemon_cablin': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Polygala_tenuifolia': {'family': 'Polygalaceae', 'order': 'Fabales'}, 'Poncirus_trifoliata': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Populus_trichocarpa': {'family': 'Salicaceae', 'order': 'Malpighiales'}, 'Populus_wilsonii': {'family': 'Salicaceae', 'order': 'Malpighiales'}, 'Portulaca_amilis': {'family': 'Portulacaceae', 'order': 'Caryophyllales'}, 'Posidonia_oceanica': {'family': 'Posidoniaceae', 'order': 'Alismatales'}, 'Potamogeton_acutifolius': {'family': 'Potamogetonaceae', 'order': 'Alismatales'}, 'Primula_veris': {'family': 'Primulaceae', 'order': 'Ericales'}, 'Protea_cynaroides': {'family': 'Proteaceae', 'order': 'Proteales'}, 'Prunus_mume': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Prunus_persica': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Punica_granatum': {'family': 'Lythraceae', 'order': 'Myrtales'}, 'Pyrus_betulifolia': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Pyrus_bretschneideri': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Quercus_acutissima': {'family': 'Fagaceae', 'order': 'Fagales'}, 'Quercus_robur': {'family': 'Fagaceae', 'order': 'Fagales'}, 'Quercus_rubra': {'family': 'Fagaceae', 'order': 'Fagales'}, 'Raddia_guianensis': {'family': 'Poaceae', 'order': 'Poales'}, 'Raphanus_sativus': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Rheum_tanguticum': {'family': 'Polygonaceae', 'order': 'Caryophyllales'}, 'Rhizophora_apiculata': {'family': 'Rhizophoraceae', 'order': 'Malpighiales'}, 'Rhodiola_crenulata': {'family': 'Crassulaceae', 'order': 'Saxifragales'}, 'Rhododendron_griersonianum': {'family': 'Ericaceae', 'order': 'Ericales'}, 'Rhododendron_simsii': {'family': 'Ericaceae', 'order': 'Ericales'}, 'Rhynchospora_breviuscula': {'family': 'Cyperaceae', 'order': 'Poales'}, 'Rhynchospora_pubera': {'family': 'Cyperaceae', 'order': 'Poales'}, 'Rhynchospora_tenuis': {'family': 'Cyperaceae', 'order': 'Poales'}, 'Ricinus_communis': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Roridula_gorgonias': {'family': 'Roridulaceae', 'order': 'Ericales'}, 'Rosa_chinensis': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Rosa_rugosa': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Rosmarinus_officinalis': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Rubus_corchorifolius': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Rubus_occidentalis': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Rumex_hastatulus': {'family': 'Polygonaceae', 'order': 'Caryophyllales'}, 'Saccharum_officinarum': {'family': 'Poaceae', 'order': 'Poales'}, 'Salix_arbutifolia': {'family': 'Salicaceae', 'order': 'Malpighiales'}, 'Salix_chaenomeloides': {'family': 'Salicaceae', 'order': 'Malpighiales'}, 'Santalum_yasi': {'family': 'Santalaceae', 'order': 'Santalales'}, 'Sapindus_mukorossi': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Saururus_chinensis': {'family': 'Saururaceae', 'order': 'Piperales'}, 'Scalesia_atractyloides': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Schrenkiella_parvula': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Sclerocarya_birrea': {'family': 'Anacardiaceae', 'order': 'Sapindales'}, 'Scutellaria_baicalensis': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Secale_cereale': {'family': 'Poaceae', 'order': 'Poales'}, 'Sedum_album': {'family': 'Crassulaceae', 'order': 'Saxifragales'}, 'Senna_tora': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Sesamum_indicum': {'family': 'Pedaliaceae', 'order': 'Lamiales'}, 'Setaria_italica': {'family': 'Poaceae', 'order': 'Poales'}, 'Setaria_viridis': {'family': 'Poaceae', 'order': 'Poales'}, 'Shorea_leprosula': {'family': 'Dipterocarpaceae', 'order': 'Malvales'}, 'Silene_latifolia': {'family': 'Caryophyllaceae', 'order': 'Caryophyllales'}, 'Simmondsia_chinensis': {'family': 'Simmondsiaceae', 'order': 'Caryophyllales'}, 'Sinapis_alba': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Siraitia_grosvenorii': {'family': 'Cucurbitaceae', 'order': 'Cucurbitales'}, 'Smallanthus_sonchifolius': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Solanum_lycopersicum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Solanum_pennellii': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Solanum_tuberosum': {'family': 'Solanaceae', 'order': 'Solanales'}, 'Sonneratia_alba': {'family': 'Lythraceae', 'order': 'Myrtales'}, 'Sonneratia_caseolaris': {'family': 'Lythraceae', 'order': 'Myrtales'}, 'Sorbus_pohuashanensis': {'family': 'Rosaceae', 'order': 'Rosales'}, 'Sorghum_bicolor': {'family': 'Poaceae', 'order': 'Poales'}, 'Spatholobus_suberectus': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Spinacia_oleracea': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Spirodela_intermedia': {'family': 'Araceae', 'order': 'Alismatales'}, 'Spirodela_polyrhiza': {'family': 'Araceae', 'order': 'Alismatales'}, 'Stellera_chamaejasme': {'family': 'Thymelaeaceae', 'order': 'Malvales'}, 'Stevia_rebaudiana': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Streptocarpus_rexii': {'family': 'Gesneriaceae', 'order': 'Lamiales'}, 'Streptochaeta_angustifolia': {'family': 'Poaceae', 'order': 'Poales'}, 'Striga_asiatica': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Striga_hermonthica': {'family': 'Orobanchaceae', 'order': 'Lamiales'}, 'Strobilanthes_cusia': {'family': 'Acanthaceae', 'order': 'Lamiales'}, 'Suaeda_aralocaspica': {'family': 'Chenopodiaceae', 'order': 'Caryophyllales'}, 'Syringa_oblata': {'family': 'Oleaceae', 'order': 'Lamiales'}, 'Syzygium_aromaticum': {'family': 'Myrtaceae', 'order': 'Myrtales'}, 'Tamarix_chinensis': {'family': 'Tamaricaceae', 'order': 'Caryophyllales'}, 'Taraxacum_kok-saghyz': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Taraxacum_mongolicum': {'family': 'Asteraceae', 'order': 'Asterales'}, 'Tarenaya_hassleriana': {'family': 'Cleomaceae', 'order': 'Brassicales'}, 'Tectona_grandis': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Telopea_speciosissima': {'family': 'Proteaceae', 'order': 'Proteales'}, 'Tetracentron_sinense': {'family': 'Trochodendraceae', 'order': 'Trochodendrales'}, 'Thalassia_testudinum': {'family': 'Hydrocharitaceae', 'order': 'Alismatales'}, 'Thalictrum_thalictroides': {'family': 'Ranunculaceae', 'order': 'Ranunculales'}, 'Thellungiella_parvula': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Theobroma_cacao': {'family': 'Malvaceae', 'order': 'Malvales'}, 'Thinopyrum_intermedium': {'family': 'Poaceae', 'order': 'Poales'}, 'Thlaspi_arvense': {'family': 'Brassicaceae', 'order': 'Brassicales'}, 'Thymus_quinquecostatus': {'family': 'Lamiaceae', 'order': 'Lamiales'}, 'Toona_sinensis': {'family': 'Meliaceae', 'order': 'Sapindales'}, 'Trapa_bicornis': {'family': 'Lythraceae', 'order': 'Myrtales'}, 'Trema_orientale': {'family': 'Cannabaceae', 'order': 'Rosales'}, 'Trichopus_zeylanicus': {'family': 'Dioscoreaceae', 'order': 'Dioscoreales'}, 'Trifolium_pratense': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Trifolium_subterraneum': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Tripterygium_wilfordii': {'family': 'Celastraceae', 'order': 'Celastrales'}, 'Triticum_aestivum': {'family': 'Poaceae', 'order': 'Poales'}, 'Triticum_turgidum': {'family': 'Poaceae', 'order': 'Poales'}, 'Trochodendron_aralioides': {'family': 'Trochodendraceae', 'order': 'Trochodendrales'}, 'Turnera_subulata': {'family': 'Passifloraceae', 'order': 'Malpighiales'}, 'Utricularia_gibba': {'family': 'Lentibulariaceae', 'order': 'Lamiales'}, 'Vaccinium_caesariense': {'family': 'Ericaceae', 'order': 'Ericales'}, 'Vaccinium_corymbosum': {'family': 'Ericaceae', 'order': 'Ericales'}, 'Vernicia_fordii': {'family': 'Euphorbiaceae', 'order': 'Malpighiales'}, 'Vicia_sativa': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Vigna_angularis': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Vigna_radiata': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Vigna_subterranea': {'family': 'Fabaceae', 'order': 'Fabales'}, 'Vitellaria_paradoxa': {'family': 'Sapotaceae', 'order': 'Ericales'}, 'Vitis_arizonica': {'family': 'Vitaceae', 'order': 'Vitales'}, 'Vitis_vinifera': {'family': 'Vitaceae', 'order': 'Vitales'}, 'Wolffia_australiana': {'family': 'Araceae', 'order': 'Alismatales'}, 'Wurfbainia_villosa': {'family': 'Zingiberaceae', 'order': 'Zingiberales'}, 'Xanthoceras_sorbifolium': {'family': 'Sapindaceae', 'order': 'Sapindales'}, 'Zantedeschia_elliottiana': {'family': 'Araceae', 'order': 'Alismatales'}, 'Zanthoxylum_armatum': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Zanthoxylum_bungeanum': {'family': 'Rutaceae', 'order': 'Sapindales'}, 'Zea_mays': {'family': 'Poaceae', 'order': 'Poales'}, 'Zingiber_officinale': {'family': 'Zingiberaceae', 'order': 'Zingiberales'}, 'Zizania_latifolia': {'family': 'Poaceae', 'order': 'Poales'}, 'Zizania_palustris': {'family': 'Poaceae', 'order': 'Poales'}, 'Ziziphus_jujuba': {'family': 'Rhamnaceae', 'order': 'Rosales'}, 'Zostera_marina': {'family': 'Zosteraceae', 'order': 'Alismatales'}, 'Zostera_muelleri': {'family': 'Zosteraceae', 'order': 'Alismatales'}, 'Zoysia_japonica': {'family': 'Poaceae', 'order': 'Poales'}, 'Zoysia_matrella': {'family': 'Poaceae', 'order': 'Poales'}, 'Brasenia_schreberi': {'family': 'Cabombaceae', 'order': 'Nymphaeales'}, 'Platycodon_grandifloras': {'family': 'Campanulaceae', 'order': 'Asterales'}, 'Salvia_splendens': {'family': 'Lamiaceae', 'order': 'Lamiales'},'Hopea_hainanensis':{'family': 'Dipterocarpaceae', 'order': 'Malvales'},'Dipterocarpus_turbinatus':{'family': 'Dipterocarpaceae', 'order': 'Malvales'},'Vatica_odorata':{'family': 'Dipterocarpaceae', 'order': 'Malvales'},'Crocus_sativus':{'family': 'Iridaceae', 'order': 'Asparagales'},'Prunella_vulgaris':{'family': 'Lamiaceae', 'order': 'Lamiales'}}

Family_clades = {'Linderniaceae':'Eudicots','Platanaceae':'Eudicots','Pontederiaceae':'Monocots','Malvaceae': 'Eudicots', 'Fabaceae': 'Eudicots', 'Velloziaceae': 'Monocots', 'Sapindaceae': 'Eudicots', 'Poaceae': 'Monocots', 'Ranunculaceae': 'Eudicots', 'Acoraceae': 'Monocots', 'Actinidiaceae': 'Eudicots', 'Plantaginaceae': 'Eudicots', 'Primulaceae': 'Eudicots', 'Brassicaceae': 'Eudicots', 'Lardizabalaceae': 'Eudicots', 'Amaryllidaceae': 'Monocots', 'Amaranthaceae': 'Eudicots', 'Amborellaceae': 'Amborellales', 'Araceae': 'Monocots', 'Bromeliaceae': 'Monocots', 'Acanthaceae': 'Eudicots', 'Apiaceae': 'Eudicots', 'Orchidaceae': 'Monocots', 'Thymelaeaceae': 'Eudicots', 'Araliaceae': 'Eudicots', 'Asteraceae': 'Eudicots', 'Arecaceae': 'Monocots', 'Rosaceae': 'Eudicots', 'Aristolochiaceae': 'Magnoliids', 'Moraceae': 'Eudicots', 'Asparagaceae': 'Monocots', 'Chenopodiaceae': 'Eudicots', 'Oxalidaceae': 'Eudicots', 'Begoniaceae': 'Eudicots', 'Cucurbitaceae': 'Eudicots', 'Betulaceae': 'Eudicots', 'Urticaceae': 'Eudicots', 'Akaniaceae': 'Eudicots', 'Rhizophoraceae': 'Eudicots', 'Scrophulariaceae': 'Eudicots', 'Buxaceae': 'Eudicots', 'Lamiaceae': 'Eudicots', 'Apocynaceae': 'Eudicots', 'Theaceae': 'Eudicots', 'Nyssaceae': 'Eudicots', 'Cannabaceae': 'Eudicots', 'Capparaceae': 'Eudicots', 'Solanaceae': 'Eudicots', 'Cyperaceae': 'Monocots', 'Caricaceae': 'Eudicots', 'Juglandaceae': 'Eudicots', 'Fagaceae': 'Eudicots', 'Casuarinaceae': 'Eudicots', 'Cephalotaceae': 'Eudicots', 'Ceratophyllaceae': 'Ceratophyllales', 'Calycanthaceae': 'Magnoliids', 'Rubiaceae': 'Eudicots', 'Chloranthaceae': 'Chloranthales', 'Salicaceae': 'Eudicots', 'Lauraceae': 'Magnoliids', 'Rutaceae': 'Eudicots', 'Cleomaceae': 'Eudicots', 'Clethraceae': 'Eudicots', 'Cornaceae': 'Eudicots', 'Papaveraceae': 'Eudicots', 'Convolvulaceae': 'Eudicots', 'Cycadaceae': 'Cycadales', 'Cymodoceaceae': 'Monocots', 'Datiscaceae': 'Eudicots', 'Caryophyllaceae': 'Eudicots', 'Dioscoreaceae': 'Monocots', 'Ebenaceae': 'Eudicots', 'Elaeagnaceae': 'Eudicots', 'Musaceae': 'Monocots', 'Berberidaceae': 'Eudicots', 'Myrtaceae': 'Eudicots', 'Euphorbiaceae': 'Eudicots', 'Nymphaeaceae': 'Nymphaeales', 'Staphyleaceae': 'Eudicots', 'Polygonaceae': 'Eudicots', 'Oleaceae': 'Eudicots', 'Lentibulariaceae': 'Eudicots', 'Polemoniaceae': 'Eudicots', 'Hamamelidaceae': 'Eudicots', 'Bignoniaceae': 'Eudicots', 'Hydrangeaceae': 'Eudicots', 'Cactaceae': 'Eudicots', 'Aquifoliaceae': 'Eudicots', 'Juncaceae': 'Monocots', 'Crassulaceae': 'Eudicots', 'Circaeasteraceae': 'Eudicots', 'Lythraceae': 'Eudicots', 'Zingiberaceae': 'Monocots', 'Orobanchaceae': 'Eudicots', 'Linaceae': 'Eudicots', 'Magnoliaceae': 'Magnoliids', 'Boraginaceae': 'Eudicots', 'Caprifoliaceae': 'Eudicots', 'Proteaceae': 'Eudicots', 'Ximeniaceae': 'Eudicots', 'Anacardiaceae': 'Eudicots', 'Melastomataceae': 'Eudicots', 'Calophyllaceae': 'Eudicots', 'Phrymaceae': 'Eudicots', 'Myricaceae': 'Eudicots', 'Moringaceae': 'Eudicots', 'Nelumbonaceae': 'Eudicots', 'Menyanthaceae': 'Eudicots', 'Paeoniaceae': 'Eudicots', 'Passifloraceae': 'Eudicots', 'Paulowniaceae': 'Eudicots', 'Piperaceae': 'Magnoliids', 'Polygalaceae': 'Eudicots', 'Portulacaceae': 'Eudicots', 'Posidoniaceae': 'Monocots', 'Potamogetonaceae': 'Monocots', 'Ericaceae': 'Eudicots', 'Roridulaceae': 'Eudicots', 'Santalaceae': 'Eudicots', 'Saururaceae': 'Magnoliids', 'Pedaliaceae': 'Eudicots', 'Dipterocarpaceae': 'Eudicots', 'Simmondsiaceae': 'Eudicots', 'Gesneriaceae': 'Eudicots', 'Tamaricaceae': 'Eudicots', 'Trochodendraceae': 'Eudicots', 'Hydrocharitaceae': 'Monocots', 'Meliaceae': 'Eudicots', 'Celastraceae': 'Eudicots', 'Sapotaceae': 'Eudicots', 'Vitaceae': 'Eudicots', 'Rhamnaceae': 'Eudicots', 'Zosteraceae': 'Monocots', 'Cabombaceae': 'Nymphaeales', 'Campanulaceae': 'Eudicots','Iridaceae':'Monocots'}

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('tree', type=click.Path(exists=True))
@click.argument('date', type=click.Path(exists=True))
def simulatewgdontree(tree,date):
    tree = Phylo.read(tree,'newick');allsp = [i.name for i in tree.get_terminals()];root_age = tree.distance(tree.get_terminals()[0])
    df = pd.read_csv(date,header=0,index_col=None,sep='\t')
    df["Species"] = df["Species"].apply(lambda x:x.replace(" ","_"))
    sporig = df.shape[0]
    df = df[df["Species"].isin(allsp)]
    df = df.drop_duplicates(subset=['WGD ID'])
    WGD_IDs_Ages = {wgdid:age for wgdid,age in zip(df["WGD ID"],df["Consensus Mean"])}
    WGD_clades_dics = getwgdlocations(tree)
    y = lambda x:" ".join(x.split("_")[1:])
    yy = lambda x:", ".join([tip.name for tip in x.get_terminals()])
    yyy = lambda x:{key.replace("α","a").replace("β","b"):value for key,value in x.items()}
    WGD_clades_dics = yyy(WGD_clades_dics)
    Num_WGD = len(WGD_clades_dics)
    dic_WGD_date = {y(wgdid):yy(mrca) for wgdid,mrca in WGD_clades_dics.items()}
    observed_mean,observed_median = np.mean(df['Consensus Mean']),np.median(df['Consensus Mean'])
    observed_WGD_dates = df['Consensus Mean'].to_list()
    observed_WGD_dates_peak = df['Consensus Peak'].to_list()
    for indice,i in enumerate(tree.get_nonterminals()): i.name = "Node_{}".format(indice)
    WGDIDs_treenodeids = {y(key):value.name for key,value in WGD_clades_dics.items()}
    Tree_nodes_Ages = getTree_nodes_Ages(WGD_IDs_Ages,WGDIDs_treenodeids)
    Tree_nodes_WGDids = getTree_nodes_WGDIDs(WGDIDs_treenodeids)
    Ages_nodes = get_nodes_ages(tree)
    WGD_rate = Num_WGD/sum([i[0] for i in Ages_nodes.values()])
    WGD_rates_categories,node_categories,rates_var_categories_ = get_rate_clades(Ages_nodes,Tree_nodes_WGDids,tree)
    observed_WGD_dates_clades = get_WGDdates_clades(node_categories,Tree_nodes_Ages)
    WGD_rate_var = get_wgd_rate_var(Ages_nodes,Tree_nodes_Ages)
    simulated_WGD_Dates,sim_WGDs,sim_means,sim_medians = [],{},[],[]
    Basic_Probability,K_PG_Probability = 0,0
    categories = sorted(list(WGD_rates_categories.keys()))
    x0 = [WGD_rates_categories[c] for c in categories]
    yx0 = lambda x:{c:v for c,v in zip(categories,x)}
    ML_WGD_rate_categories = minimize(lambda x: -get_prob_Background(Ages_nodes,Tree_nodes_Ages,Tree_nodes_WGDids,WGD_rate,lognormal=False,sample_var=WGD_rate_var,rate_catetory=yx0(x),node_catetory=node_categories,rates_var_categories=rates_var_categories_),x0,bounds = [[0,None] for i in range(len(WGD_rates_categories))],method='L-BFGS-B')
    WGD_rates_categories = ML_WGD_rates_categories = {c:r for c,r in zip(categories,ML_WGD_rate_categories.x)}
    Basic_Probability_background_fun = minimize(lambda x: -get_prob_Background(Ages_nodes,Tree_nodes_Ages,Tree_nodes_WGDids,x),[WGD_rate],bounds=[[0,None]],method='L-BFGS-B')
    Basic_Probability_background_ML_Rate = Basic_Probability_background_fun.x[0]
    Basic_Probability_background = -Basic_Probability_background_fun.fun[0]
    Basic_Probability_background_autocor = get_prob_Background(Ages_nodes,Tree_nodes_Ages,Tree_nodes_WGDids,WGD_rate,lognormal=False,sample_var=WGD_rate_var,rate_catetory=WGD_rates_categories,node_catetory=node_categories,rates_var_categories=rates_var_categories_)
    AIC_Basic_background = 2*1 - 2*Basic_Probability_background
    AIC_Variable = 2*2 - 2*Basic_Probability_background_autocor
    AIC_dics = {"Constant-rate":AIC_Basic_background,"Variable-rate":AIC_Variable}
    best_AIC = sorted(AIC_dics.items(), key = lambda x:x[1])[0]
    p_value = LRT_comparison(Basic_Probability_background_autocor,Basic_Probability_background,1)
    pvalues_log10_reversed = -np.log10(p_value)
    K_PG_Probabilitys,olds,pvalues,alpha_MLs = [],[],[],[]
    alpha_ = 1
    log_likelihood_ratios = []
    for old,young in zip(np.arange(0,205,5),np.arange(5,210,5)):
        old,young = 205-old,205-young
        Time_range_test = (old,young)
        alpha_ML_ = minimize(lambda x: -get_prob_K_pg(x,Ages_nodes,Tree_nodes_Ages,WGD_rate,Time_range_test,lognormal=False,sample_var=WGD_rate_var,rate_catetory=WGD_rates_categories,node_catetory=node_categories,rates_var_categories=rates_var_categories_), [1], bounds = [[1,None]],method='L-BFGS-B')
        prob_ML = -alpha_ML_.fun[0]
        alpha_ML = alpha_ML_.x[0]
        K_PG_Probability_control = get_prob_K_pg(1,Ages_nodes,Tree_nodes_Ages,WGD_rate,Time_range_test,lognormal=False,sample_var=WGD_rate_var,rate_catetory=WGD_rates_categories,node_catetory=node_categories,rates_var_categories=rates_var_categories_)
        if prob_ML <= K_PG_Probability_control: alpha_ML = 1
        if LRT_comparison(prob_ML,K_PG_Probability_control,1) >= 0.0001: alpha_ML = 1
        alpha_MLs.append(alpha_ML)
        p_value = LRT_comparison(prob_ML,K_PG_Probability_control,1)
        log_likelihood_ratios.append(-2*(K_PG_Probability_control-prob_ML))
        olds.append(old);K_PG_Probabilitys.append(K_PG_Probability);pvalues.append(p_value)
    writepvaluealpha(pvalues,alpha_MLs,olds,log_likelihood_ratios=log_likelihood_ratios)
    simulated_WGD_dates_ = []
    y2 = lambda x:", ".join([str(i) for i in x])
    Iterations = 100
    for time in range(Iterations):
        simulated_WGD_Dates = []
        for key,value in Ages_nodes.items():
            duration_time,start_time = value[0],value[2]
            simulated_WGD_dates = simulation_WGD(key,WGD_rate,duration_time,start_time,lognormal=False,sample_var=WGD_rate_var,rate_catetory=WGD_rates_categories,node_catetory=node_categories,rates_var_categories=rates_var_categories_)
            sim_WGDs[key] = [y2(simulated_WGD_dates)] if sim_WGDs.get(key) is None else sim_WGDs[key] + [y2(simulated_WGD_dates)]
            for i in simulated_WGD_dates: simulated_WGD_Dates.append(i)
        simulated_WGD_dates_.append(simulated_WGD_Dates)
        sim_means.append(np.mean(simulated_WGD_Dates))
        sim_medians.append(np.median(simulated_WGD_Dates))
    assert len(simulated_WGD_dates_) == Iterations
    plotcombined_complete_paper(simulated_WGD_dates_,observed_WGD_dates,alpha_MLs,pvalues,olds,outfile="Simulation_Results.pdf")
    df = pd.DataFrame.from_dict(sim_WGDs)
    df.index = ["Iteration_{}".format(i) for i in range(df.shape[0])]
    df.to_csv("Multi_Simulation_WGDates_{}.tsv".format(Iterations),header=True,index=True,sep='\t')
    writesimudates(simulated_WGD_dates_,observed_WGD_dates)

def get_WGDdates_clades(node_categories,Tree_nodes_Ages):
    clades = {clade:[] for clade in set(node_categories.values())}
    for nodeid,wgdates in Tree_nodes_Ages.items():
        clade = node_categories[nodeid]
        clades[clade] += wgdates
    return clades

def writesimudates(simulated_WGD_dates_,observed_WGD_dates):
    X = list(itertools.chain(*simulated_WGD_dates_))
    labels_X = ["simulation" for i in X]
    labels_O = ["observation" for i in observed_WGD_dates]
    combined = X + observed_WGD_dates
    labels_combined = labels_X + labels_O
    df_X = pd.DataFrame(data=combined, index=labels_combined, columns=['WGD date'])
    df_X.to_csv("Simulation_Iterations_Observation.tsv",header=True,index=True,sep='\t')

def processolds(olds):
    new_olds = []
    for old in olds:
        num = int(old)
        new_olds += ["-".join([str(num),str(num-5)])]
    return new_olds
    
def writepvaluealpha(pvalues,alpha_MLs,olds,log_likelihood_ratios=None):
    olds = processolds(olds)
    df_alpha = pd.DataFrame(data=alpha_MLs, index=olds, columns=['Pulled ratio'])
    df_alpha.index.name = "Time window"
    df_alpha.to_csv("Pulled_ratio_per_window.tsv",header=True,index=True,sep='\t')
    pvalues = [-np.log10(i) for i in pvalues]
    df_pvalue = pd.DataFrame(data=pvalues, index=olds, columns=['P-values'])
    df_pvalue.index.name = "Time window"
    df_pvalue.to_csv("Pvalue_per_window.tsv",header=True,index=True,sep='\t')
    if log_likelihood_ratios is not None:
        df_llr = pd.DataFrame(data=log_likelihood_ratios, index=olds, columns=['log-likelihood ratio'])
        df_llr.index.name = "Time window"
        df_llr.to_csv("log_likelihood_ratio_per_window.tsv",header=True,index=True,sep='\t')

def plotcombined_complete_paper(simulated_dates,observed_dates,alpha_MLs,pvalues,olds,outfile="Simulation_Results.pdf",presim=None):
    fig,ax1,ax2,ax3,ax4,ax5,ax6 = constructfigbasic6p()
    if presim is not None:
        df = pd.read_csv(presim,header=0,index_col=0,sep='\t',low_memory=False).fillna("")
        X = []
        y = lambda x:[float(i) for i in x]
        for i in df.index:
            dates_per_iter = []
            for dates in df.loc[i,:].dropna().astype(str):
                if dates == "": continue
                dates_per_iter += y(dates.split(", "))
            X.append(dates_per_iter)
        simulated_dates = X
    plot_sim_fulldis_paper_ax(ax1,simulated_dates,observed_dates,alpha=0.4)
    plot_sim_testMD_ax(ax2,simulated_dates,observed_dates,(40,0),alpha=0.4)
    plot_sim_testMD_ax(ax3,simulated_dates,observed_dates,(80,40),alpha=0.4)
    plot_sim_testMD_ax(ax4,simulated_dates,observed_dates,(120,80),alpha=0.4)
    plot_chipvalue_ax(ax5,pvalues,olds,noxlabel=False,lw=1)
    plot_MLalpha_ax(ax6,pvalues,alpha_MLs,olds,noxlabel=False,lw=1)
    plt.tight_layout()
    plt.savefig(outfile,format ='pdf', bbox_inches='tight')
    plt.close()

def getmedian(observed_dates):
    differences = [abs(d1-d2) for d1,d2 in itertools.combinations(observed_dates, 2)]
    return np.median(differences)

def plot_sim_testMD_ax(ax,simulated_WGD_dates_,observed_WGD_dates,range_,alpha=0.8):
    olds = []
    judges = {True:"Lower",False:"Higher"}
    old,young = range_
    simu_mds = []
    for simulated_WGD_dates in simulated_WGD_dates_:
        simulated_WGD_dates_bounded = [d for d in simulated_WGD_dates if young<=d<=old]
        simu_md = getmedian(simulated_WGD_dates_bounded)
        if simu_md>=0: simu_mds.append(simu_md)
    observed_WGD_dates_bounded = [d for d in observed_WGD_dates if young<=d<=old]
    observed_WGD_dates_md = getmedian(observed_WGD_dates_bounded)
    X = simu_mds
    observed_median = observed_WGD_dates_md
    if not observed_median >=0: observed_median=0
    lower,upper = min([min(X),observed_median])*0.9,max([observed_median,max(X)])*1.1
    assert lower<upper
    assert upper-lower > 0.1
    if upper <= observed_median+2: upper = observed_median+2
    Hs, Bins, patches = ax.hist(X, bins = np.arange(lower,upper,0.1), color='green', alpha=0, rwidth=0.8, density = True)
    kde_x = np.linspace(lower,upper,num=10000)
    kde_y=stats.gaussian_kde(X,bw_method='silverman').pdf(kde_x)
    mode, maxim = kde_mode(kde_x, kde_y)
    CHF = get_totalH(Hs)
    scaling = CHF*0.1
    ax.plot(kde_x, kde_y*scaling, color='k',alpha=0.8, ls = '-',lw = 1)
    ax.fill_between(kde_x, kde_y*scaling, color='gray',alpha=alpha,label="Simulation")
    ax.axvline(observed_median,color='g', ls='--', lw=1,label="Observation")
    ax.set_title("Test of cluster {}-{}mya".format(young,old),fontsize = 15)
    ax.legend(fontsize=15,loc=1,frameon=False)
    ax.set_xlim(lower,upper)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.set_ylabel("Density", fontsize = 15)
    ax.set_xlabel("Median WGD distance (my)", fontsize = 15)

def plot_sim_fulldis_paper_ax(ax,X,observed_dates,colors=['black','green'],lw=1,alpha=0.8):
    colors = [adjust_saturation(x,0.8) for x in colors]
    X = list(itertools.chain(*X))
    lower,upper = min(X)*0.9,max(X)*1.1
    Hs, Bins, patches = ax.hist(X, bins = np.arange(lower,upper,1), color=colors[0], alpha=0, rwidth=0.8, density=True)
    kde_x = np.linspace(lower,upper,num=500)
    kde_y=stats.gaussian_kde(X,bw_method='silverman').pdf(kde_x)
    mode, maxim = kde_mode(kde_x, kde_y)
    CHF = get_totalH(Hs)
    scaling = CHF*1
    max_1 = max(kde_y*scaling)
    kde_y_simu = kde_y*scaling
    ax.plot(kde_x, kde_y*scaling, color=colors[0],alpha=0.8, ls = '-',lw = lw)
    ax.fill_between(kde_x, kde_y*scaling, color='gray',alpha=alpha,label="Simulation")
    Hs, Bins, patches = ax.hist(observed_dates, bins = np.arange(lower,upper,1), color='green', alpha=0, rwidth=0.8, density=True)
    kde_y =stats.gaussian_kde(observed_dates,bw_method='silverman').pdf(kde_x)
    CHF = get_totalH(Hs)
    scaling = CHF*1
    ax.plot(kde_x, kde_y*scaling, color=colors[1],alpha=0.8, ls = '-',lw = lw)
    max_2 = max(kde_y*scaling)
    kde_y_obse = kde_y*scaling
    ax.fill_between(kde_x, kde_y*scaling, color=colors[1],alpha=0.2, label="Observation")
    statistic, p_value = ks_2samp(kde_y_simu, kde_y_obse)
    ax.plot([],[],alpha=0,label='P-value: {:.4f}'.format(p_value))
    handles, labels = ax.get_legend_handles_labels()
    order_list = ["Simulation","Observation",'P-value: {:.4f}'.format(p_value)]
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda x: order_list.index(x[0])))
    ax.legend(handles,labels,fontsize=15,loc=1,frameon=False)
    ax.set_xlim(lower,upper)
    ax.set_title("Nonrandom temporal distribution of WGDs",fontsize = 15)
    ax.set_ylim(0,0.02)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.set_ylabel("Density", fontsize = 15)
    ax.set_xlabel("WGD dates (mya)", fontsize = 15)

def constructfigbasic6p():
    fig = plt.figure(figsize=(16, 18))
    gs = GridSpec(3, 2, figure=fig, wspace=0.15, hspace=0.20, left=0.08, right=0.92, top=0.95, bottom=0.08)
    axes = []
    for i in range(3):
        for j in range(2):
            axes.append(fig.add_subplot(gs[i, j]))
    ax1,ax2,ax3,ax4,ax5,ax6 = axes
    return fig,ax1,ax2,ax3,ax4,ax5,ax6

def get_rate_clades(Ages_nodes,Tree_nodes_WGDids,tree):
    y = lambda x:taxonomy_prset[x.replace(" ","_")]['family']
    Family_clades2 = modifiedF_c()
    rates_categories = {}
    rates_categories_ = {}
    rates_var_categories_ = {}
    node_categories = {}
    for key,value in Ages_nodes.items():
        nodes = list(tree.get_nonterminals()) + list(tree.get_terminals())
        nodes_dic = {node.name:node for node in nodes}
        node = nodes_dic[key]
        spnames = [sp.name for sp in node.get_terminals()]
        clades = set([Family_clades2[y(spname)] for spname in spnames])
        if len(clades) == 1:
            clade = list(clades)[0]
            node_categories[key] = clade
            if clade not in rates_categories:
                rates_categories[clade] = {"Duration":[value[0]],"Num_WGDs":[len(Tree_nodes_WGDids.get(key,[]))]}
            else:
                rates_categories[clade]["Duration"].append(value[0])
                rates_categories[clade]["Num_WGDs"].append(len(Tree_nodes_WGDids.get(key,[])))
        else:
            clade = "Background"
            node_categories[key] = clade
            if clade not in rates_categories:
                rates_categories[clade] = {"Duration":[value[0]],"Num_WGDs":[len(Tree_nodes_WGDids.get(key,[]))]}
                rates_categories[clade]["Num_WGDs"].append(len(Tree_nodes_WGDids.get(key,[])))
    for clade,dic in rates_categories.items():
        rates_categories_[clade] = sum(dic["Num_WGDs"])/sum(dic["Duration"])
        rate_per_subclade = [occur/time for occur,time in zip(dic["Num_WGDs"],dic["Duration"])]
        rates_var_categories_[clade] = np.var(rate_per_subclade)
    return rates_categories_, node_categories, rates_var_categories_

def get_wgd_rate_var(Ages_nodes,Tree_nodes_Ages):
    rates = []
    for key,value in Ages_nodes.items():
        duration = value[0]
        occurance = len(Tree_nodes_Ages.get(key,[]))
        rate = occurance/duration
        rates.append(rate)
    var = np.var(rates)
    return var


def plottestbranch(Ages_nodes,Tree_nodes_WGDids,nodes_occurance):
    p_values = []
    olds = []
    for age in range(0,205,5):
        old_bound = age+5;young_bound = age
        olds.append(old_bound)
        observed_occurence = 0
        simu_occurences = []
        for key,value in Tree_nodes_WGDids.items():
            stem_age = Ages_nodes[key][-1]
            if stem_age > old_bound or stem_age < young_bound:
                continue
            observed_occurence += len(value)
            simu_occurences += nodes_occurance[key]
        t_statistic, p_value = ttest_1samp(simu_occurences, popmean=observed_occurence)
        p_values.append(p_value)


def addvvline(ax,xvalue,color,lstyle,labell):
    if labell == '': ax.axvline(xvalue,color=color, ls=lstyle, lw=1)
    elif labell == 'Observation': ax.axvline(xvalue,color=color, ls=lstyle, lw=1)
    else: ax.axvline(xvalue,color=color, ls=lstyle, lw=1, label='{}: {:.1f}'.format(labell,xvalue))
    return ax

def addhhline(ax,yvalue,color,lstyle,labell):
    if labell == '': ax.axhline(yvalue,color=color, ls=lstyle, lw=1)
    elif labell == 'Observation': ax.axhline(yvalue,color=color, ls=lstyle, lw=1)
    else: ax.axhline(yvalue,color=color, ls=lstyle, lw=1, label='{}: {:.1f}'.format(labell,xvalue))
    return ax

def plot_chipvalue_ax(ax,pvalues,olds,noxlabel=True,lw=5):
    pvalues_log10_reversed = -np.log10(pvalues)
    for v,o in zip(pvalues_log10_reversed,olds):
        if v < 4: ax.scatter(o,v, c='gray', ls='-', s=4, alpha=0.8)
        else: ax.scatter(o,v, c='red', ls='-', s=16, alpha=0.8)
    ax.plot(olds,pvalues_log10_reversed, c='k', ls='-', lw=lw, alpha=1)
    ax.axhline(4,color='r', ls='--', lw=1,label = "P-value: 0.0001")
    ax.set_ylabel(r'-$\log_{10}($p-value$)$',fontsize=15)
    if not noxlabel: ax.set_xlabel("Date (mya)",fontsize=15)
    max_,min_ = np.ceil(max(pvalues_log10_reversed)),np.floor(min(pvalues_log10_reversed))
    ax.set_ylim(min_-1,max_+1)
    ax.set_title("Model test of pulled rate",fontsize=15)
    ax.legend(loc=1,fontsize=15,frameon=False)

def plot_MLalpha_ax(ax,pvalues,alpha_MLs,olds,noxlabel=False,lw=5):
    pvalues_log10_reversed = -np.log10(pvalues)
    for v,o,a in zip(pvalues_log10_reversed,olds,alpha_MLs):
        if v < 4: ax.scatter(o,a, c='gray', ls='-', s=4, alpha=0.8)
        else: ax.scatter(o,a, c='red', ls='-', s=16, alpha=0.8)
    ax.plot(olds,np.array(alpha_MLs), c='k', ls='-', lw=lw, alpha=1)
    ax = addvvline(ax,66.0,'k','--','K-Pg event')
    ax = addvvline(ax,55.8,'k','--','PETM event')
    ax.legend(loc=0,fontsize=15,frameon=False)
    ax.set_title("Optimized pulled ratio",fontsize=15)
    ax.set_ylabel("Pulled ratio",fontsize=15)
    if not noxlabel: ax.set_xlabel("Date (mya)",fontsize=15)
    max_,min_ = np.ceil(max(alpha_MLs)),np.floor(min(alpha_MLs))
    ax.set_yticks(np.arange(min_, max_+1, 1))

def simulation_WGD(node,WGD_rate,duration_time,start_time,lognormal=False,sample_var=None,rate_catetory=None,node_catetory=None,rates_var_categories=None):
    if rate_catetory is not None and node_catetory is not None and rates_var_categories is not None:
        catetory = node_catetory[node]
        WGD_rate = clade_rate = rate_catetory[catetory]
        sample_var = rates_var_categories[catetory]
        if WGD_rate <=0: WGD_rate = 1e-11
    if lognormal:
        sigma = np.sqrt(sample_var)
        if WGD_rate <=0: WGD_rate = 1e-11
        if sigma <=0: sigma = 1e-11
        sampled_WGD_rate = lognorm.rvs(s=sigma, scale=WGD_rate, size=1)
        WGD_rate = sampled_WGD_rate[0]
    lambda_ = WGD_rate
    current_time = start_time
    event_dates = []
    while current_time >= (start_time - duration_time):
        time_to_next_event = np.random.exponential(1/lambda_)
        current_time = current_time - time_to_next_event
        if (start_time - duration_time) <= current_time <= start_time:
            event_dates.append(current_time)
    return event_dates

def filterWGDate(dic,younger_bound,older_bound):
    dic_filtered = {key:[] for key in dic.keys()}
    for key,value in dic.items():
        for date in value:
            if younger_bound <= date <= older_bound:
                dic_filtered[key].append(date)
    return dic_filtered

def LRT_comparison(alternative,null,df):
    LR = -2*(null-alternative)
    p_value = stats.chi2.sf(LR, df)
    return p_value

def get_prob_K_pg(alpha,Ages_nodes,Tree_nodes_Ages,WGD_rate,Time_range_test,lognormal=False,sample_var=None,rate_catetory=None,node_catetory=None,rates_var_categories=None):
    K_PG_Probability = 0
    for key,value in Ages_nodes.items():
        duration_time,start_time,end_time = value[0],value[2],value[1]
        K_PG_prob = Kpg_model_WGD(alpha,key,WGD_rate,Tree_nodes_Ages,duration_time,start_time,end_time,Time_range_test,lognormal=lognormal,sample_var=sample_var,rate_catetory=rate_catetory,node_catetory=node_catetory,rates_var_categories=rates_var_categories)
        K_PG_Probability += K_PG_prob
    return K_PG_Probability

def Kpg_model_WGD(alpha,node,WGD_rate,Tree_nodes_Ages,duration_time,start_time,end_time,Time_range_test,lognormal=False,sample_var=None,rate_catetory=None,node_catetory=None,rates_var_categories=None):
    Probability = 0
    Conditions = {'a':0,'b':0,'c':0,'d':0,'e':0}
    if rate_catetory is not None and node_catetory is not None and rates_var_categories is not None:
        catetory = node_catetory[node]
        WGD_rate = clade_rate = rate_catetory[catetory]
        if WGD_rate<=0: WGD_rate = 1e-11
        sample_var = rates_var_categories[catetory]
    if lognormal:
        sigma = np.sqrt(sample_var)
        if sigma <=0: sigma = 1e-11
        sampled_WGD_rate = lognorm.rvs(s=sigma, scale=WGD_rate, size=1)
        pdf_sample = lognorm.logpdf(sampled_WGD_rate[0], s=sigma, scale=WGD_rate)
        Probability += pdf_sample
        WGD_rate = sampled_WGD_rate[0]
    if start_time <= Time_range_test[1] or end_time >= Time_range_test[0]:
        Probability += model_WGD(node,WGD_rate,Tree_nodes_Ages,duration_time,start_time)
        Conditions['a'] +=1
    elif Time_range_test[1]<=start_time<=Time_range_test[0] and Time_range_test[1]<=end_time<=Time_range_test[0]:
        Probability += model_WGD(node,WGD_rate*alpha,Tree_nodes_Ages,duration_time,start_time)
        Conditions['b'] +=1
    elif end_time<= Time_range_test[1] <= start_time and end_time<= Time_range_test[0] <= start_time:
        duration_time_first_Background = start_time - Time_range_test[0]
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,Time_range_test[0],start_time)
        Probability += model_WGD(node,WGD_rate,Tree_nodes_Ages_filter,duration_time_first_Background,start_time)
        duration_time_second_Pulled = Time_range_test[0] - Time_range_test[1]
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,Time_range_test[1],Time_range_test[0])
        Probability += model_WGD(node,WGD_rate*alpha,Tree_nodes_Ages_filter,duration_time_second_Pulled,Time_range_test[0])
        duration_time_third_Background = Time_range_test[1] - end_time
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,end_time,Time_range_test[1])
        Probability += model_WGD(node,WGD_rate,Tree_nodes_Ages_filter,duration_time_third_Background,Time_range_test[1])
        Conditions['c'] +=1
    elif end_time<=Time_range_test[0]<=start_time and Time_range_test[1] <= end_time:
        duration_time_first_Background = start_time - Time_range_test[0]
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,Time_range_test[0],start_time)
        Probability += model_WGD(node,WGD_rate,Tree_nodes_Ages_filter,duration_time_first_Background,start_time)
        duration_time_second_Pulled = Time_range_test[0] - end_time
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,end_time,Time_range_test[0])
        Probability += model_WGD(node,WGD_rate*alpha,Tree_nodes_Ages_filter,duration_time_second_Pulled,Time_range_test[0])
        Conditions['d'] +=1
    elif end_time<=Time_range_test[1]<=start_time and start_time<=Time_range_test[0]:
        duration_time_first_Pulled = start_time - Time_range_test[1]
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,Time_range_test[1],start_time)
        Probability += model_WGD(node,WGD_rate*alpha,Tree_nodes_Ages_filter,duration_time_first_Pulled,start_time)
        duration_time_second_Background = Time_range_test[1] - end_time
        Tree_nodes_Ages_filter = filterWGDate(Tree_nodes_Ages,end_time,Time_range_test[1])
        Probability += model_WGD(node,WGD_rate,Tree_nodes_Ages_filter,duration_time_second_Background,Time_range_test[1])
        Conditions['e'] +=1
    return Probability

def getTree_nodes_WGDIDs(WGDIDs_treenodeids):
    Tree_nodes_WGDids = {}
    for wgdid,nodeid in WGDIDs_treenodeids.items():
        if nodeid in Tree_nodes_WGDids:
            Tree_nodes_WGDids[nodeid] += [wgdid]
        else:
            Tree_nodes_WGDids[nodeid] = [wgdid]
    return Tree_nodes_WGDids

def getTree_nodes_Ages(WGD_IDs_Ages,WGDIDs_treenodeids):
    Tree_nodes_Ages = {}
    for wgdid,age in WGD_IDs_Ages.items():
        treenodeid = WGDIDs_treenodeids[wgdid]
        if treenodeid in Tree_nodes_Ages:
            Tree_nodes_Ages[treenodeid].append(age)
        else: 
            Tree_nodes_Ages[treenodeid] = [age]
    return Tree_nodes_Ages

def getfirstparent(target,tree):
    All_parents = [(node,node.distance(target)) for node in tree.get_nonterminals() if node.is_parent_of(target) and node != target]
    return sorted(All_parents,key = lambda x:x[1])[0][0]

def calcomment(comment):
    comment = comment.replace("&95%HPD={","").replace("}","")
    younger, older = comment.split(", ")
    return float(younger), float(older)

def get_nodes_ages_withuncertainty(tree):
    nodes = [i for i in tree.get_nonterminals()]
    Ages_nodes = {}
    for clade in nodes:
        if type(clade.branch_length) is not float:
            continue
        first_parent_node = getfirstparent(clade,tree)
        older_bound_pre = max(clade.depths().values())
        younger_CI,older_CI = calcomment(first_parent_node.comment)
        older_bound = older_CI
        younger_CI,older_CI = calcomment(clade.comment)
        younger_bound_pre = older_bound_pre - clade.branch_length
        younger_bound = younger_CI
        duration = older_bound - younger_bound
        Ages_nodes[clade.name] = (duration*100,younger_bound*100,older_bound*100)
    tips = [i for i in tree.get_terminals()]
    for clade in tips:
        younger_bound = 0
        older_bound = clade.branch_length
        Ages_nodes[clade.name] = (clade.branch_length*100,younger_bound*100,older_bound*100)
    return Ages_nodes

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

def get_prob_Background(Ages_nodes,Tree_nodes_Ages,Tree_nodes_WGDids,WGD_rate,lognormal=False,sample_var=None,rate_catetory=None,node_catetory=None,rates_var_categories=None):
    Basic_Probability = 0
    for node,value in Ages_nodes.items():
        duration_time,start_time,end_time = value[0],value[2],value[1]
        assert len(Tree_nodes_WGDids.get(node,[])) == len(Tree_nodes_Ages.get(node,[]))
        if rate_catetory is not None and node_catetory is not None and rates_var_categories is not None:
            catetory = node_catetory[node]
            WGD_rate = clade_rate = rate_catetory[catetory]
            sample_var = rates_var_categories[catetory]
            if WGD_rate <=0: WGD_rate = 1e-11
        if lognormal:
            sigma = np.sqrt(sample_var)
            if sigma <=0: sigma = 1e-11
            sampled_WGD_rate = lognorm.rvs(s=sigma, scale=WGD_rate, size=1)
            pdf_sample = lognorm.logpdf(sampled_WGD_rate[0], s=sigma, scale=WGD_rate)
            Basic_Probability += pdf_sample
            prob = model_WGD(node,sampled_WGD_rate[0],Tree_nodes_Ages,duration_time,start_time)
        else: prob = model_WGD(node,WGD_rate,Tree_nodes_Ages,duration_time,start_time)
        Basic_Probability += prob
    return Basic_Probability

def model_WGD(node,WGD_rate,Tree_nodes_Ages,duration_time,start_time):
    lambda_ = WGD_rate * duration_time
    dates = Tree_nodes_Ages.get(node,[])
    Num_WGD_observed = len(dates)
    probability_poisson = np.log(poisson_pmf(Num_WGD_observed, lambda_))
    probability_exponential = 0
    for date in dates:
        occur_interval = start_time - date
        if occur_interval <0:
            continue
        probability = np.log(exponential_pdf(occur_interval, lambda_))
        probability_exponential += probability
    probability = probability_poisson + probability_exponential
    return probability

def exponential_pdf(x, lambda_):
    return lambda_ * np.exp(-lambda_ * x)

def getwgdlocations(tree):
    WGD_clade1 = {"WGD1_MALE":tree.common_ancestor({'name':'Pyrus_bretschneideri'},{'name':'Crataegus_pinnatifida'}),"WGD2_PAPI":tree.common_ancestor({'name':'Ammopiptanthus_nanus'},{'name':'Vigna_radiata'}),"WGD2_CAES":tree.common_ancestor({'name':'Acacia_pycnantha'},{'name':'Senna_tora'}),"WGD3_GLYC":tree.common_ancestor({'name':'Glycine_soja'},{'name':'Glycine_max'}),"WGD4_JUGL":tree.common_ancestor({'name':'Carya_illinoinensis'},{'name':'Juglans_regia'}),"WGD5_BRCE":tree.common_ancestor({'name':'Sinapis_alba'},{'name':'Brassica_oleracea'}),"WGD6_BRAS_α":tree.common_ancestor({'name':'Aethionema_arabicum'},{'name':'Brassica_oleracea'}),"WGD7_BRAS_β":tree.common_ancestor({'name':'Capparis_spinosa'},{'name':'Brassica_oleracea'}),"WGD9_TARE":tree.common_ancestor({'name':'Tarenaya_hassleriana'},{'name':'Tarenaya_hassleriana'}),"WGD12_MALV":tree.common_ancestor({'name':'Gossypium_barbadense'},{'name':'Durio_zibethinus'}),"WGD13_MAHE":tree.common_ancestor({'name':'Hevea_brasiliensis'},{'name':'Manihot_esculenta'}),"WGD14_POSA":tree.common_ancestor({'name':'Populus_wilsonii'},{'name':'Salix_arbutifolia'}),"WGD15_LINU_α":tree.common_ancestor({'name':'Linum_usitatissimum'},{'name':'Linum_usitatissimum'}),"WGD15_LINU_β":tree.common_ancestor({'name':'Linum_usitatissimum'},{'name':'Linum_usitatissimum'}),"WGD16_MYRT":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Sonneratia_caseolaris'}),"WGD17_LYTH":tree.common_ancestor({'name':'Sonneratia_caseolaris'},{'name':'Lagerstroemia_indica'}),"WGD19_TRAP":tree.common_ancestor({'name':'Trapa_bicornis'},{'name':'Trapa_bicornis'}),"WGD20_SOLA":tree.common_ancestor({'name':'Capsicum_chinense'},{'name':'Petunia_axillaris'}),"WGD21_OLEA_α":tree.common_ancestor({'name':'Syringa_oblata'},{'name':'Olea_europaea'}),"WGD22_OLEA_β":tree.common_ancestor({'name':'Olea_europaea'},{'name':'Jasminum_sambac'}),"WGD23_ASTE":tree.common_ancestor({'name':'Carthamus_tinctorius'},{'name':'Helianthus_annuus'}),"WGD24_APIO_α":tree.common_ancestor({'name':'Daucus_carota'},{'name':'Coriandrum_sativum'}),"WGD25_APIO_β":tree.common_ancestor({'name':'Panax_ginseng'},{'name':'Daucus_carota'}),"WGD27_AMAR":tree.common_ancestor({'name':'Amaranthus_hypochondriacus'},{'name':'Amaranthus_tuberculatus'}),"WGD28_KALA":tree.common_ancestor({'name':'Kalanchoe_fedtschenkoi'},{'name':'Kalanchoe_fedtschenkoi'}),"WGD29_RANU":tree.common_ancestor({'name':'Aquilegia_oxysepala'},{'name':'Corydalis_tomentella'}),"WGD30_POAC":tree.common_ancestor({'name':'Triticum_aestivum'},{'name':'Streptochaeta_angustifolia'}),"WGD31_ZEAM":tree.common_ancestor({'name':'Zea_mays'},{'name':'Zea_mays'}),"WGD32_POAL":tree.common_ancestor({'name':'Triticum_aestivum'},{'name':'Ananas_comosus'}),"WGD33_MUSA_α":tree.common_ancestor({'name':'Ensete_ventricosum'},{'name':'Musa_acuminata'}),"WGD34_MUSA_β":tree.common_ancestor({'name':'Musa_acuminata'},{'name':'Zingiber_officinale'})}
    WGD_clade2 = {"WGD35_AREC":tree.common_ancestor({'name':'Elaeis_guineensis'},{'name':'Calamus_simplicifolius'}),"WGD38_ASPA":tree.common_ancestor({'name':'Asparagus_setaceus'},{'name':'Asparagus_officinalis'}),"WGD40_ARAC_α":tree.common_ancestor({'name':'Zantedeschia_elliottiana'},{'name':'Lemna_minuta'}),"WGD41_ARAC_β":tree.common_ancestor({'name':'Zantedeschia_elliottiana'},{'name':'Lemna_minuta'}),"WGD42_LAUR_α":tree.common_ancestor({'name':'Phoebe_bournei'},{'name':'Cinnamomum_kanehirae'}),"WGD43_LAUR_β":tree.common_ancestor({'name':'Phoebe_bournei'},{'name':'Magnolia_biondii'}),"WGD44_NYMP":tree.common_ancestor({'name':'Nymphaea_colorata'},{'name':'Brasenia_schreberi'}),"WGD44_BRSH_α":tree.common_ancestor({'name':'Brasenia_schreberi'},{'name':'Brasenia_schreberi'}),"WGD44_BRSH_β":tree.common_ancestor({'name':'Brasenia_schreberi'},{'name':'Brasenia_schreberi'}),"WGD45_ACOR":tree.common_ancestor({'name':'Acorus_americanus'},{'name':'Acorus_tatarinowii'}),"WGD46_DIOS":tree.common_ancestor({'name':'Dioscorea_alata'},{'name':'Trichopus_zeylanicus'}),"WGD47_TRIC":tree.common_ancestor({'name':'Trichopus_zeylanicus'},{'name':'Trichopus_zeylanicus'}),"WGD48_ACAN":tree.common_ancestor({'name':'Acanthochlamys_bracteata'},{'name':'Acanthochlamys_bracteata'}),"WGD50_CUCU":tree.common_ancestor({'name':'Cucumis_melo'},{'name':'Datisca_glomerata'}),"WGD51_CUBI":tree.common_ancestor({'name':'Cucurbita_maxima'},{'name':'Cucurbita_argyrosperma'}),"WGD52_BEGO":tree.common_ancestor({'name':'Begonia_fuchsioides'},{'name':'Begonia_darthvaderiana'}),"WGD52_BEFU":tree.common_ancestor({'name':'Begonia_fuchsioides'},{'name':'Begonia_fuchsioides'}),"WGD53_LUPI":tree.common_ancestor({'name':'Lupinus_angustifolius'},{'name':'Lupinus_angustifolius'}),"WGD54_CONV":tree.common_ancestor({'name':'Ipomoea_nil'},{'name':'Cuscuta_australis'}),"WGD57_ACTI_α":tree.common_ancestor({'name':'Actinidia_chinensis'},{'name':'Actinidia_eriantha'}),"WGD58_ACTI_β":tree.common_ancestor({'name':'Rhododendron_simsii'},{'name':'Gilia_yorkii'}),"WGD59_NELU":tree.common_ancestor({'name':'Nelumbo_nucifera'},{'name':'Nelumbo_nucifera'}),"WGD60_PASE_β":tree.common_ancestor({'name':'Papaver_somniferum'},{'name':'Papaver_setigerum'}),"WGD61_PASE_α":tree.common_ancestor({'name':'Papaver_setigerum'},{'name':'Papaver_setigerum'}),"WGD63_ORCH":tree.common_ancestor({'name':'Cremastra_appendiculata'},{'name':'Apostasia_shenzhenica'}),"WGD64_ZANT":tree.common_ancestor({'name':'Zanthoxylum_bungeanum'},{'name':'Zanthoxylum_armatum'}),"WGD65_EURY":tree.common_ancestor({'name':'Euryale_ferox'},{'name':'Euryale_ferox'}),"WGD66_PIPE":tree.common_ancestor({'name':'Piper_nigrum'},{'name':'Piper_nigrum'}),"WGD67_CALY":tree.common_ancestor({'name':'Chimonanthus_salicifolius'},{'name':'Chimonanthus_salicifolius'})}
    WGD_clade3 = {"WGD68_CHLO":tree.common_ancestor({'name':'Chloranthus_spicatus'},{'name':'Chloranthus_sessilifolius'}),"WGD68_LARD":tree.common_ancestor({'name':'Akebia_trifoliata'},{'name':'Akebia_trifoliata'}),"WGD69_BUXA":tree.common_ancestor({'name':'Buxus_sinica'},{'name':'Buxus_austroyunnanensis'}),"WGD70_TROC_α":tree.common_ancestor({'name':'Trochodendron_aralioides'},{'name':'Tetracentron_sinense'}),"WGD70_TROC_β":tree.common_ancestor({'name':'Trochodendron_aralioides'},{'name':'Tetracentron_sinense'}),"WGD71_ACHN":tree.common_ancestor({'name':'Achnatherum_splendens'},{'name':'Achnatherum_splendens'}),"WGD72_BONI":tree.common_ancestor({'name':'Bonia_amplexicaulis'},{'name':'Bonia_amplexicaulis'}),"WGD76_PHRA":tree.common_ancestor({'name':'Phragmites_australis'},{'name':'Phragmites_australis'}),"WGD77_ZIZA":tree.common_ancestor({'name':'Zizania_latifolia'},{'name':'Zizania_palustris'}),"WGD78_ZING":tree.common_ancestor({'name':'Zingiber_officinale'},{'name':'Lanxangia_tsaoko'}),"WGD79_POLG":tree.common_ancestor({'name':'Fagopyrum_tataricum'},{'name':'Rheum_tanguticum'}),"WGD80_CARY":tree.common_ancestor({'name':'Gypsophila_paniculata'},{'name':'Silene_latifolia'}),"WGD81_ACPT":tree.common_ancestor({'name':'Portulaca_amilis'},{'name':'Hylocereus_undatus'}),"WGD82_NYSS":tree.common_ancestor({'name':'Davidia_involucrata'},{'name':'Nyssa_yunnanensis'}),"WGD83_PRIM":tree.common_ancestor({'name':'Aegiceras_corniculatum'},{'name':'Primula_veris'}),"WGD84_UTRI":tree.common_ancestor({'name':'Utricularia_gibba'},{'name':'Utricularia_gibba'}),"WGD85_OROB":tree.common_ancestor({'name':'Origanum_majorana'},{'name':'Buddleja_alternifolia'}),"WGD86_STRI":tree.common_ancestor({'name':'Striga_hermonthica'},{'name':'Striga_asiatica'}),"WGD87_STRE":tree.common_ancestor({'name':'Streptocarpus_rexii'},{'name':'Streptocarpus_rexii'}),"WGD88_ANTI":tree.common_ancestor({'name':'Antirrhinum_majus'},{'name':'Antirrhinum_majus'}),"WGD89_SALV":tree.common_ancestor({'name':'Salvia_splendens'},{'name':'Salvia_splendens'}),"WGD91_BORA_α":tree.common_ancestor({'name':'Lithospermum_erythrorhizon'},{'name':'Lithospermum_erythrorhizon'}),"WGD92_BORA_β":tree.common_ancestor({'name':'Lithospermum_erythrorhizon'},{'name':'Lithospermum_erythrorhizon'}),"WGD93_SMAL":tree.common_ancestor({'name':'Smallanthus_sonchifolius'},{'name':'Smallanthus_sonchifolius'}),"WGD95_HELI":tree.common_ancestor({'name':'Mikania_micrantha'},{'name':'Helianthus_annuus'}),"WGD97_ARAL":tree.common_ancestor({'name':'Eleutherococcus_senticosus'},{'name':'Panax_ginseng'}),"WGD98_ELEU":tree.common_ancestor({'name':'Eleutherococcus_senticosus'},{'name':'Eleutherococcus_senticosus'}),"WGD99_DIPS":tree.common_ancestor({'name':'Lonicera_japonica'},{'name':'Lonicera_japonica'})}
    WGD_clade4 = {"WGD100_AQUI":tree.common_ancestor({'name':'Ilex_polyneura'},{'name':'Ilex_polyneura'}),"WGD101_BAUH":tree.common_ancestor({'name':'Bauhinia_variegata'},{'name':'Bauhinia_variegata'}),"WGD102_HIPP_α":tree.common_ancestor({'name':'Hippophae_rhamnoides'},{'name':'Elaeagnus_angustifolia'}),"WGD102_HIPP_β":tree.common_ancestor({'name':'Hippophae_rhamnoides'},{'name':'Elaeagnus_angustifolia'}),"WGD103_TRIP":tree.common_ancestor({'name':'Tripterygium_wilfordii'},{'name':'Tripterygium_wilfordii'}),"WGD104_PASS":tree.common_ancestor({'name':'Passiflora_edulis'},{'name':'Passiflora_organensis'}),"WGD105_RHIZ":tree.common_ancestor({'name':'Kandelia_candel'},{'name':'Bruguiera_sexangula'}),"WGD106_MELA_α":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Melastoma_dodecandrum'}),"WGD107_MELA_β":tree.common_ancestor({'name':'Melastoma_dodecandrum'},{'name':'Melastoma_dodecandrum'}),"WGD108_HIBI":tree.common_ancestor({'name':'Hibiscus_cannabinus'},{'name':'Abelmoschus_esculentus'}),"WGD109_DIPT":tree.common_ancestor({'name':'Shorea_leprosula'},{'name':'Shorea_leprosula'}),"WGD110_THYM":tree.common_ancestor({'name':'Aquilaria_sinensis'},{'name':'Stellera_chamaejasme'}),"WGD111_STAP":tree.common_ancestor({'name':'Euscaphis_japonica'},{'name':'Euscaphis_japonica'}),"WGD112_BRET":tree.common_ancestor({'name':'Bretschneidera_sinensis'},{'name':'Bretschneidera_sinensis'}),"WGD113_CAPP":tree.common_ancestor({'name':'Capparis_spinosa'},{'name':'Capparis_spinosa'}),"WGD115_LOME":tree.common_ancestor({'name':'Lobularia_maritima'},{'name':'Megacarpaea_delavayi'}),"WGD116_ORYC":tree.common_ancestor({'name':'Orychophragmus_violaceus'},{'name':'Orychophragmus_violaceus'}),"WGD117_TOON":tree.common_ancestor({'name':'Toona_sinensis'},{'name':'Toona_sinensis'}),"WGD119_MANG":tree.common_ancestor({'name':'Mangifera_indica'},{'name':'Mangifera_indica'}),"WGD120_PROT":tree.common_ancestor({'name':'Macadamia_jansenii'},{'name':'Protea_cynaroides'}),"WGD122_CYMO":tree.common_ancestor({'name':'Cymodocea_nodosa'},{'name':'Cymodocea_nodosa'}),"WGD123_POZO":tree.common_ancestor({'name':'Potamogeton_acutifolius'},{'name':'Zostera_muelleri'}),"WGD124_SEAG":tree.common_ancestor({'name':'Thalassia_testudinum'},{'name':'Zostera_muelleri'}),"WGD125_CERA_α":tree.common_ancestor({'name':'Ceratophyllum_demersum'},{'name':'Ceratophyllum_demersum'}),"WGD126_CERA_β":tree.common_ancestor({'name':'Ceratophyllum_demersum'},{'name':'Ceratophyllum_demersum'}),"WGD126b_TAMA":tree.common_ancestor({'name':'Tamarix_chinensis'},{'name':'Tamarix_chinensis'}),"WGD126c_SAPI":tree.common_ancestor({'name':'Saururus_chinensis'},{'name':'Piper_nigrum'}),"WGD126d_SANT":tree.common_ancestor({'name':'Santalum_yasi'},{'name':'Santalum_yasi'}),"WGD127_LEPI":tree.common_ancestor({'name':'Lepidium_sativum'},{'name':'Lepidium_sativum'}),"WGD128_ACON":tree.common_ancestor({'name':'Aconitum_vilmorinianum'},{'name':'Aconitum_vilmorinianum'}),"WGD131_POLY_α":tree.common_ancestor({'name':'Polygala_tenuifolia'},{'name':'Polygala_tenuifolia'}),"WGD132_POLY_β":tree.common_ancestor({'name':'Polygala_tenuifolia'},{'name':'Polygala_tenuifolia'})}
    WGD_clade5 = {"WGD133_MESU":tree.common_ancestor({'name':'Mesua_ferrea'},{'name':'Mesua_ferrea'}),"WGD135_CORN":tree.common_ancestor({'name':'Cornus_wilsoniana'},{'name':'Cornus_wilsoniana'}),"WGD136_AESC":tree.common_ancestor({'name':'Aesculus_chinensis'},{'name':'Aesculus_chinensis'}),"WGD137_ABEL":tree.common_ancestor({'name':'Abelmoschus_esculentus'},{'name':'Abelmoschus_esculentus'}),"WGD138_RHSE":tree.common_ancestor({'name':'Rhodiola_crenulata'},{'name':'Sedum_album'}),"WGD141_ROSM":tree.common_ancestor({'name':'Rosmarinus_officinalis'},{'name':'Rosmarinus_officinalis'}),"WGD142_OCIM":tree.common_ancestor({'name':'Ocimum_basilicum'},{'name':'Ocimum_basilicum'}),"WGD143_NEOL":tree.common_ancestor({'name':'Neolamarckia_cadamba'},{'name':'Neolamarckia_cadamba'}),"WGD144_COEU":tree.common_ancestor({'name':'Tamarix_chinensis'},{'name':'Vigna_radiata'}),"WGD145_COMO":tree.common_ancestor({'name':'Acanthochlamys_bracteata'},{'name':'Triticum_aestivum'}),"WGD146_BOTH":tree.common_ancestor({'name':'Bothriochloa_decipiens'},{'name':'Bothriochloa_decipiens'}),"WGD147_ZOYS":tree.common_ancestor({'name':'Zoysia_matrella'},{'name':'Zoysia_japonica'})}
    return {**WGD_clade1,**WGD_clade2,**WGD_clade3,**WGD_clade4,**WGD_clade5}

if __name__ == '__main__':
        cli()
