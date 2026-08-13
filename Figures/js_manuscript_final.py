# This contains all the supporting functions utilized in the Scheib et. al., 2026 manuscript
# Changelog:
#   2026-08-07 (claude-opus-5): removed redundant get_allAnmParams3(); documented get_allAnmParams()
#                               vs get_allAnmParams2() (verbose/extended-phase vs quiet/standard-phase)
#   2026-08-10 (claude-opus-5): added reduce_summaryInfo_for_trial_structure() to prune summaryInfo down
#                               to the fields load_trial_structure_times() reads; it now lives in
#                               sled_summarize_anms.py
#   2026-08-11 (claude-opus-5): added load_all_grouped_aligned_traces_split() (plus its _align_open /
#                               _align_restore_floats helpers) to read back the compressed
#                               all_grouped_aligned_traces slices written by
#                               jsm_figs.save_all_grouped_aligned_traces_split,
#                               so the *_indiv_examples notebooks can load their example ROIs' trace data,
#                               cueInfo and summaryInfo straight out of the repo; its trace_types filter
#                               resolves the 'preShiftOnly_byTtype'/'postShiftOnly_byTtype' display modes to
#                               the stored 'byPrePost_byTtype' branch they are a view of, via
#                               resolve_align_trace_type(), so filtering by one no longer prunes the data away
###################################################
import os
import warnings
from datetime import datetime, date
import sys, numpy.core
sys.modules.setdefault('numpy._core', numpy.core)
sys.modules.setdefault('numpy._core.multiarray', numpy.core.multiarray)
sys.modules.setdefault('numpy._core.numeric', numpy.core.numeric)
import numpy as np
from tqdm.auto import tqdm
import copy
import statsmodels.api as sm
from scipy import stats
import time
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.stats import skew,mode,gmean,ranksums,ttest_ind, zscore, kstest, pearsonr
import random
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
###################################################
plt.rcParams['font.family']           = 'sans-serif'
plt.rcParams['font.sans-serif']       = ['Arial']
plt.rcParams['pdf.fonttype']          = 42       # embed fonts in PDF (Type 42)
plt.rcParams['axes.unicode_minus']    = False    # use ASCII minus sign
plt.rcParams['mathtext.default']      = 'regular'
plt.rcParams['animation.embed_limit'] = 2**128
##################################################################################################
def masterData_loading(DATA_DIR):
    masterData = {}
    masterData['B00002121749'] = {}
    masterData['B00002121749']['behavior'] = {}
    masterData['B00002121749']['shiftInfo'] = {'anm': 'B00002121749', 'shiftDay': 221015, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002121749']['somas'] = {}
    masterData['B00002121749']['somas']['daysRelativeToShift'] = {}
    masterData['B00002121749']['somas']['roiScores'] = {}
    masterData['B00002121749']['somas']['sessions'] = {}
    masterData['B00002121749']['somas']['sessions']['0'] = {}
    masterData['B00002121749']['somas']['sessions']['1'] = {}
    masterData['B00002121749']['somas']['sessions']['2'] = {}
    masterData['B00002121749']['somas']['sessions']['3'] = {}
    masterData['B00002121749']['somas']['sessions']['4'] = {}
    masterData['B00002121749']['dendrites'] = 'NOT RECORDED'
    masterData['B00002121774'] = {}
    masterData['B00002121774']['behavior'] = {}
    masterData['B00002121774']['shiftInfo'] = {'anm': 'B00002121774', 'shiftDay': 221012, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002121774']['somas'] = {}
    masterData['B00002121774']['somas']['daysRelativeToShift'] = {}
    masterData['B00002121774']['somas']['roiScores'] = {}
    masterData['B00002121774']['somas']['sessions'] = {}
    masterData['B00002121774']['somas']['sessions']['0'] = {}
    masterData['B00002121774']['somas']['sessions']['1'] = {}
    masterData['B00002121774']['somas']['sessions']['2'] = {}
    masterData['B00002121774']['somas']['sessions']['3'] = {}
    masterData['B00002121774']['somas']['sessions']['4'] = {}
    masterData['B00002121774']['dendrites'] = {}
    masterData['B00002121774']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002121774']['dendrites']['roiScores'] = {}
    masterData['B00002121774']['dendrites']['sessions'] = {}
    masterData['B00002121774']['dendrites']['sessions']['0'] = {}
    masterData['B00002121774']['dendrites']['sessions']['1'] = {}
    masterData['B00002121774']['dendrites']['sessions']['2'] = {}
    masterData['B00002121774']['dendrites']['sessions']['3'] = {}
    masterData['B00002121777'] = {}
    masterData['B00002121777']['behavior'] = {}
    masterData['B00002121777']['shiftInfo'] = {'anm': 'B00002121777', 'shiftDay': 221012, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002121777']['somas'] = {}
    masterData['B00002121777']['somas']['daysRelativeToShift'] = {}
    masterData['B00002121777']['somas']['roiScores'] = {}
    masterData['B00002121777']['somas']['sessions'] = {}
    masterData['B00002121777']['somas']['sessions']['0'] = {}
    masterData['B00002121777']['somas']['sessions']['1'] = {}
    masterData['B00002121777']['somas']['sessions']['2'] = {}
    masterData['B00002121777']['somas']['sessions']['3'] = {}
    masterData['B00002121777']['somas']['sessions']['4'] = {}
    masterData['B00002121777']['dendrites'] = {}
    masterData['B00002121777']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002121777']['dendrites']['roiScores'] = {}
    masterData['B00002121777']['dendrites']['sessions'] = {}
    masterData['B00002121777']['dendrites']['sessions']['0'] = {}
    masterData['B00002121777']['dendrites']['sessions']['1'] = {}
    masterData['B00002121777']['dendrites']['sessions']['2'] = {}
    masterData['B00002121777']['dendrites']['sessions']['3'] = {}
    masterData['B00002121777']['dendrites']['sessions']['4'] = {}
    masterData['B00002213772'] = {}
    masterData['B00002213772']['behavior'] = {}
    masterData['B00002213772']['shiftInfo'] = {'anm': 'B00002213772', 'shiftDay': 221223, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213772']['somas'] = {}
    masterData['B00002213772']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213772']['somas']['roiScores'] = {}
    masterData['B00002213772']['somas']['sessions'] = {}
    masterData['B00002213772']['somas']['sessions']['0'] = {}
    masterData['B00002213772']['somas']['sessions']['1'] = {}
    masterData['B00002213772']['somas']['sessions']['2'] = {}
    masterData['B00002213772']['somas']['sessions']['3'] = {}
    masterData['B00002213772']['somas']['sessions']['4'] = {}
    masterData['B00002213772']['dendrites'] = {}
    masterData['B00002213772']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213772']['dendrites']['roiScores'] = {}
    masterData['B00002213772']['dendrites']['sessions'] = {}
    masterData['B00002213772']['dendrites']['sessions']['0'] = {}
    masterData['B00002213772']['dendrites']['sessions']['1'] = {}
    masterData['B00002213772']['dendrites']['sessions']['2'] = {}
    masterData['B00002213772']['dendrites']['sessions']['3'] = {}
    masterData['B00002213772']['dendrites']['sessions']['4'] = {}
    masterData['B00002213772']['dendrites']['sessions']['5'] = {}
    masterData['B00002213773'] = {}
    masterData['B00002213773']['behavior'] = {}
    masterData['B00002213773']['shiftInfo'] = {'anm': 'B00002213773', 'shiftDay': 'na', 'shiftTrial': 'na', 'shiftDir': 'na'}
    masterData['B00002213773']['somas'] = {}
    masterData['B00002213773']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213773']['somas']['roiScores'] = {}
    masterData['B00002213773']['somas']['sessions'] = {}
    masterData['B00002213773']['somas']['sessions']['0'] = {}
    masterData['B00002213773']['dendrites'] = {}
    masterData['B00002213773']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213773']['dendrites']['roiScores'] = {}
    masterData['B00002213773']['dendrites']['sessions'] = {}
    masterData['B00002213773']['dendrites']['sessions']['0'] = {}
    masterData['B00002213773']['dendrites']['sessions']['1'] = {}
    masterData['B00002213784'] = {}
    masterData['B00002213784']['behavior'] = {}
    masterData['B00002213784']['shiftInfo'] = {'anm': 'B00002213784', 'shiftDay': 230129, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213784']['somas'] = {}
    masterData['B00002213784']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213784']['somas']['roiScores'] = {}
    masterData['B00002213784']['somas']['sessions'] = {}
    masterData['B00002213784']['somas']['sessions']['0'] = {}
    masterData['B00002213784']['somas']['sessions']['1'] = {}
    masterData['B00002213784']['somas']['sessions']['2'] = {}
    masterData['B00002213784']['somas']['sessions']['3'] = {}
    masterData['B00002213784']['dendrites'] = {}
    masterData['B00002213784']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213784']['dendrites']['roiScores'] = {}
    masterData['B00002213784']['dendrites']['sessions'] = {}
    masterData['B00002213784']['dendrites']['sessions']['0'] = {}
    masterData['B00002213784']['dendrites']['sessions']['1'] = {}
    masterData['B00002213784']['dendrites']['sessions']['2'] = {}
    masterData['B00002213784']['dendrites']['sessions']['3'] = {}
    masterData['B00002213784']['dendrites']['sessions']['4'] = {}
    masterData['B00002213785'] = {}
    masterData['B00002213785']['behavior'] = {}
    masterData['B00002213785']['shiftInfo'] = {'anm': 'B00002213785', 'shiftDay': 230131, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213785']['somas'] = {}
    masterData['B00002213785']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213785']['somas']['roiScores'] = {}
    masterData['B00002213785']['somas']['sessions'] = {}
    masterData['B00002213785']['somas']['sessions']['0'] = {}
    masterData['B00002213785']['somas']['sessions']['1'] = {}
    masterData['B00002213785']['somas']['sessions']['2'] = {}
    masterData['B00002213785']['somas']['sessions']['3'] = {}
    masterData['B00002213785']['dendrites'] = {}
    masterData['B00002213785']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213785']['dendrites']['roiScores'] = {}
    masterData['B00002213785']['dendrites']['sessions'] = {}
    masterData['B00002213785']['dendrites']['sessions']['0'] = {}
    masterData['B00002213785']['dendrites']['sessions']['1'] = {}
    masterData['B00002213785']['dendrites']['sessions']['2'] = {}
    masterData['B00002213785']['dendrites']['sessions']['3'] = {}
    masterData['B00002213785']['dendrites']['sessions']['4'] = {}
    masterData['B00002213889'] = {}
    masterData['B00002213889']['behavior'] = {}
    masterData['B00002213889']['shiftInfo'] = {'anm': 'B00002213889', 'shiftDay': 230920, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213889']['somas'] = {}
    masterData['B00002213889']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213889']['somas']['roiScores'] = {}
    masterData['B00002213889']['somas']['sessions'] = {}
    masterData['B00002213889']['somas']['sessions']['0'] = {}
    masterData['B00002213889']['somas']['sessions']['1'] = {}
    masterData['B00002213889']['somas']['sessions']['2'] = {}
    masterData['B00002213889']['somas']['sessions']['3'] = {}
    masterData['B00002213889']['dendrites'] = {}
    masterData['B00002213889']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213889']['dendrites']['roiScores'] = {}
    masterData['B00002213889']['dendrites']['sessions'] = {}
    masterData['B00002213889']['dendrites']['sessions']['0'] = {}
    masterData['B00002213889']['dendrites']['sessions']['1'] = {}
    masterData['B00002213889']['dendrites']['sessions']['2'] = {}
    masterData['B00002213908'] = {}
    masterData['B00002213908']['behavior'] = {}
    masterData['B00002213908']['shiftInfo'] = {'anm': 'B00002213908', 'shiftDay': 230920, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213908']['somas'] = {}
    masterData['B00002213908']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213908']['somas']['roiScores'] = {}
    masterData['B00002213908']['somas']['sessions'] = {}
    masterData['B00002213908']['somas']['sessions']['0'] = {}
    masterData['B00002213908']['somas']['sessions']['1'] = {}
    masterData['B00002213908']['somas']['sessions']['2'] = {}
    masterData['B00002213908']['somas']['sessions']['3'] = {}
    masterData['B00002213908']['somas']['sessions']['4'] = {}
    masterData['B00002213908']['dendrites'] = {}
    masterData['B00002213908']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213908']['dendrites']['roiScores'] = {}
    masterData['B00002213908']['dendrites']['sessions'] = {}
    masterData['B00002213908']['dendrites']['sessions']['0'] = {}
    masterData['B00002213908']['dendrites']['sessions']['1'] = {}
    masterData['B00002213908']['dendrites']['sessions']['2'] = {}
    masterData['B00002213908']['dendrites']['sessions']['3'] = {}
    masterData['B00002213909'] = {}
    masterData['B00002213909']['behavior'] = {}
    masterData['B00002213909']['shiftInfo'] = {'anm': 'B00002213909', 'shiftDay': 'na', 'shiftTrial': 'na', 'shiftDir': 'na'}
    masterData['B00002213909']['somas'] = {}
    masterData['B00002213909']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213909']['somas']['roiScores'] = {}
    masterData['B00002213909']['somas']['sessions'] = {}
    masterData['B00002213909']['somas']['sessions']['0'] = {}
    masterData['B00002213909']['somas']['sessions']['1'] = {}
    masterData['B00002213909']['dendrites'] = {}
    masterData['B00002213909']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213909']['dendrites']['roiScores'] = {}
    masterData['B00002213909']['dendrites']['sessions'] = {}
    masterData['B00002213909']['dendrites']['sessions']['0'] = {}
    masterData['B00002213909']['dendrites']['sessions']['1'] = {}
    masterData['B00002213909']['dendrites']['sessions']['2'] = {}
    masterData['B00002213920'] = {}
    masterData['B00002213920']['behavior'] = {}
    masterData['B00002213920']['shiftInfo'] = {'anm': 'B00002213920', 'shiftDay': 231024, 'shiftTrial': 102, 'shiftDir': 0}
    masterData['B00002213920']['somas'] = {}
    masterData['B00002213920']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213920']['somas']['roiScores'] = {}
    masterData['B00002213920']['somas']['sessions'] = {}
    masterData['B00002213920']['somas']['sessions']['0'] = {}
    masterData['B00002213920']['somas']['sessions']['1'] = {}
    masterData['B00002213920']['somas']['sessions']['2'] = {}
    masterData['B00002213920']['somas']['sessions']['3'] = {}
    masterData['B00002213920']['dendrites'] = {}
    masterData['B00002213920']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213920']['dendrites']['roiScores'] = {}
    masterData['B00002213920']['dendrites']['sessions'] = {}
    masterData['B00002213920']['dendrites']['sessions']['0'] = {}
    masterData['B00002213920']['dendrites']['sessions']['1'] = {}
    masterData['B00002213920']['dendrites']['sessions']['2'] = {}
    masterData['B00002213920']['dendrites']['sessions']['3'] = {}
    masterData['B00002213921'] = {}
    masterData['B00002213921']['behavior'] = {}
    masterData['B00002213921']['shiftInfo'] = {'anm': 'B00002213921', 'shiftDay': 231114, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213921']['somas'] = {}
    masterData['B00002213921']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213921']['somas']['roiScores'] = {}
    masterData['B00002213921']['somas']['sessions'] = {}
    masterData['B00002213921']['somas']['sessions']['0'] = {}
    masterData['B00002213921']['somas']['sessions']['1'] = {}
    masterData['B00002213921']['somas']['sessions']['2'] = {}
    masterData['B00002213921']['dendrites'] = {}
    masterData['B00002213921']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213921']['dendrites']['roiScores'] = {}
    masterData['B00002213921']['dendrites']['sessions'] = {}
    masterData['B00002213921']['dendrites']['sessions']['0'] = {}
    masterData['B00002213921']['dendrites']['sessions']['1'] = {}
    masterData['B00002213921']['dendrites']['sessions']['2'] = {}
    masterData['B00002213932'] = {}
    masterData['B00002213932']['behavior'] = {}
    masterData['B00002213932']['shiftInfo'] = {'anm': 'B00002213932', 'shiftDay': 231116, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213932']['somas'] = {}
    masterData['B00002213932']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213932']['somas']['roiScores'] = {}
    masterData['B00002213932']['somas']['sessions'] = {}
    masterData['B00002213932']['somas']['sessions']['0'] = {}
    masterData['B00002213932']['somas']['sessions']['1'] = {}
    masterData['B00002213932']['somas']['sessions']['2'] = {}
    masterData['B00002213932']['dendrites'] = {}
    masterData['B00002213932']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213932']['dendrites']['roiScores'] = {}
    masterData['B00002213932']['dendrites']['sessions'] = {}
    masterData['B00002213932']['dendrites']['sessions']['0'] = {}
    masterData['B00002213932']['dendrites']['sessions']['1'] = {}
    masterData['B00002213932']['dendrites']['sessions']['2'] = {}
    masterData['B00002213932']['dendrites']['sessions']['3'] = {}
    masterData['B00002213932']['dendrites']['sessions']['4'] = {}
    masterData['B00002213943'] = {}
    masterData['B00002213943']['behavior'] = {}
    masterData['B00002213943']['shiftInfo'] = {'anm': 'B00002213943', 'shiftDay': 231203, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213943']['somas'] = {}
    masterData['B00002213943']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213943']['somas']['roiScores'] = {}
    masterData['B00002213943']['somas']['sessions'] = {}
    masterData['B00002213943']['somas']['sessions']['0'] = {}
    masterData['B00002213943']['somas']['sessions']['1'] = {}
    masterData['B00002213943']['somas']['sessions']['2'] = {}
    masterData['B00002213943']['dendrites'] = {}
    masterData['B00002213943']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213943']['dendrites']['roiScores'] = {}
    masterData['B00002213943']['dendrites']['sessions'] = {}
    masterData['B00002213943']['dendrites']['sessions']['0'] = {}
    masterData['B00002213943']['dendrites']['sessions']['1'] = {}
    masterData['B00002213943']['dendrites']['sessions']['2'] = {}
    masterData['B00002213943']['dendrites']['sessions']['3'] = {}
    masterData['B00002213944'] = {}
    masterData['B00002213944']['behavior'] = {}
    masterData['B00002213944']['shiftInfo'] = {'anm': 'B00002213944', 'shiftDay': 231225, 'shiftTrial': 105, 'shiftDir': 0}
    masterData['B00002213944']['somas'] = 'NOT RECORDED'
    masterData['B00002213944']['dendrites'] = {}
    masterData['B00002213944']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213944']['dendrites']['roiScores'] = {}
    masterData['B00002213944']['dendrites']['sessions'] = {}
    masterData['B00002213944']['dendrites']['sessions']['0'] = {}
    masterData['B00002213944']['dendrites']['sessions']['1'] = {}
    masterData['B00002213970'] = {}
    masterData['B00002213970']['behavior'] = {}
    masterData['B00002213970']['shiftInfo'] = {'anm': 'B00002213970', 'shiftDay': 240218, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213970']['somas'] = {}
    masterData['B00002213970']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213970']['somas']['roiScores'] = {}
    masterData['B00002213970']['somas']['sessions'] = {}
    masterData['B00002213970']['somas']['sessions']['0'] = {}
    masterData['B00002213970']['somas']['sessions']['1'] = {}
    masterData['B00002213970']['dendrites'] = {}
    masterData['B00002213970']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213970']['dendrites']['roiScores'] = {}
    masterData['B00002213970']['dendrites']['sessions'] = {}
    masterData['B00002213970']['dendrites']['sessions']['0'] = {}
    masterData['B00002213970']['dendrites']['sessions']['1'] = {}
    masterData['B00002213970']['dendrites']['sessions']['2'] = {}
    masterData['B00002213970']['dendrites']['sessions']['3'] = {}
    masterData['B00002213970']['dendrites']['sessions']['4'] = {}
    masterData['B00002213970']['dendrites']['sessions']['5'] = {}
    masterData['B00002213970']['dendrites']['sessions']['6'] = {}
    masterData['B00002213970']['dendrites']['sessions']['7'] = {}
    masterData['B00002213982'] = {}
    masterData['B00002213982']['behavior'] = {}
    masterData['B00002213982']['shiftInfo'] = {'anm': 'B00002213982', 'shiftDay': 240218, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213982']['somas'] = {}
    masterData['B00002213982']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213982']['somas']['roiScores'] = {}
    masterData['B00002213982']['somas']['sessions'] = {}
    masterData['B00002213982']['somas']['sessions']['0'] = {}
    masterData['B00002213982']['somas']['sessions']['1'] = {}
    masterData['B00002213982']['somas']['sessions']['2'] = {}
    masterData['B00002213982']['somas']['sessions']['3'] = {}
    masterData['B00002213982']['dendrites'] = {}
    masterData['B00002213982']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213982']['dendrites']['roiScores'] = {}
    masterData['B00002213982']['dendrites']['sessions'] = {}
    masterData['B00002213982']['dendrites']['sessions']['0'] = {}
    masterData['B00002213982']['dendrites']['sessions']['1'] = {}
    masterData['B00002213982']['dendrites']['sessions']['2'] = {}
    masterData['B00002213982']['dendrites']['sessions']['3'] = {}
    masterData['B00002213985'] = {}
    masterData['B00002213985']['behavior'] = {}
    masterData['B00002213985']['shiftInfo'] = {'anm': 'B00002213985', 'shiftDay': 240307, 'shiftTrial': 103, 'shiftDir': 0}
    masterData['B00002213985']['somas'] = {}
    masterData['B00002213985']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213985']['somas']['roiScores'] = {}
    masterData['B00002213985']['somas']['sessions'] = {}
    masterData['B00002213985']['somas']['sessions']['0'] = {}
    masterData['B00002213985']['somas']['sessions']['1'] = {}
    masterData['B00002213985']['somas']['sessions']['2'] = {}
    masterData['B00002213985']['dendrites'] = {}
    masterData['B00002213985']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213985']['dendrites']['roiScores'] = {}
    masterData['B00002213985']['dendrites']['sessions'] = {}
    masterData['B00002213985']['dendrites']['sessions']['0'] = {}
    masterData['B00002213985']['dendrites']['sessions']['1'] = {}
    masterData['B00002213985']['dendrites']['sessions']['2'] = {}
    masterData['B00002213985']['dendrites']['sessions']['3'] = {}
    masterData['B00002213997'] = {}
    masterData['B00002213997']['behavior'] = {}
    masterData['B00002213997']['shiftInfo'] = {'anm': 'B00002213997', 'shiftDay': 240325, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213997']['somas'] = {}
    masterData['B00002213997']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213997']['somas']['roiScores'] = {}
    masterData['B00002213997']['somas']['sessions'] = {}
    masterData['B00002213997']['somas']['sessions']['0'] = {}
    masterData['B00002213997']['somas']['sessions']['1'] = {}
    masterData['B00002213997']['somas']['sessions']['2'] = {}
    masterData['B00002213997']['somas']['sessions']['3'] = {}
    masterData['B00002213997']['dendrites'] = {}
    masterData['B00002213997']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213997']['dendrites']['roiScores'] = {}
    masterData['B00002213997']['dendrites']['sessions'] = {}
    masterData['B00002213997']['dendrites']['sessions']['0'] = {}
    masterData['B00002213997']['dendrites']['sessions']['1'] = {}
    masterData['B00002213997']['dendrites']['sessions']['2'] = {}
    masterData['B00002213997']['dendrites']['sessions']['3'] = {}
    masterData['B00002213997']['dendrites']['sessions']['4'] = {}
    masterData['B00002213998'] = {}
    masterData['B00002213998']['behavior'] = {}
    masterData['B00002213998']['shiftInfo'] = {'anm': 'B00002213998', 'shiftDay': 240323, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213998']['somas'] = 'NOT RECORDED'
    masterData['B00002213998']['dendrites'] = {}
    masterData['B00002213998']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213998']['dendrites']['roiScores'] = {}
    masterData['B00002213998']['dendrites']['sessions'] = {}
    masterData['B00002213998']['dendrites']['sessions']['0'] = {}
    masterData['B00002213998']['dendrites']['sessions']['1'] = {}
    masterData['B00002213998']['dendrites']['sessions']['2'] = {}
    masterData['B00002213998']['dendrites']['sessions']['3'] = {}
    masterData['B00002213998']['dendrites']['sessions']['4'] = {}
    masterData['B00002213998']['dendrites']['sessions']['5'] = {}
    masterData['B00002213998']['dendrites']['sessions']['6'] = {}
    masterData['B00002213999'] = {}
    masterData['B00002213999']['behavior'] = {}
    masterData['B00002213999']['shiftInfo'] = {'anm': 'B00002213999', 'shiftDay': 240324, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002213999']['somas'] = {}
    masterData['B00002213999']['somas']['daysRelativeToShift'] = {}
    masterData['B00002213999']['somas']['roiScores'] = {}
    masterData['B00002213999']['somas']['sessions'] = {}
    masterData['B00002213999']['somas']['sessions']['0'] = {}
    masterData['B00002213999']['somas']['sessions']['1'] = {}
    masterData['B00002213999']['somas']['sessions']['2'] = {}
    masterData['B00002213999']['somas']['sessions']['3'] = {}
    masterData['B00002213999']['dendrites'] = {}
    masterData['B00002213999']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002213999']['dendrites']['roiScores'] = {}
    masterData['B00002213999']['dendrites']['sessions'] = {}
    masterData['B00002213999']['dendrites']['sessions']['0'] = {}
    masterData['B00002213999']['dendrites']['sessions']['1'] = {}
    masterData['B00002213999']['dendrites']['sessions']['2'] = {}
    masterData['B00002213999']['dendrites']['sessions']['3'] = {}
    masterData['B00002213999']['dendrites']['sessions']['4'] = {}
    masterData['B00002214001'] = {}
    masterData['B00002214001']['behavior'] = {}
    masterData['B00002214001']['shiftInfo'] = {'anm': 'B00002214001', 'shiftDay': 240413, 'shiftTrial': 100, 'shiftDir': 0}
    masterData['B00002214001']['somas'] = {}
    masterData['B00002214001']['somas']['daysRelativeToShift'] = {}
    masterData['B00002214001']['somas']['roiScores'] = {}
    masterData['B00002214001']['somas']['sessions'] = {}
    masterData['B00002214001']['somas']['sessions']['0'] = {}
    masterData['B00002214001']['somas']['sessions']['1'] = {}
    masterData['B00002214001']['somas']['sessions']['2'] = {}
    masterData['B00002214001']['somas']['sessions']['3'] = {}
    masterData['B00002214001']['somas']['sessions']['4'] = {}
    masterData['B00002214001']['dendrites'] = {}
    masterData['B00002214001']['dendrites']['daysRelativeToShift'] = {}
    masterData['B00002214001']['dendrites']['roiScores'] = {}
    masterData['B00002214001']['dendrites']['sessions'] = {}
    masterData['B00002214001']['dendrites']['sessions']['0'] = {}
    masterData['B00002214001']['dendrites']['sessions']['1'] = {}
    masterData['B00002214001']['dendrites']['sessions']['2'] = {}
    masterData['B00002214001']['dendrites']['sessions']['3'] = {}


    for anm in masterData.keys():
        print('--------------------------')
        print(f"Loading data for {anm}...")
        for field in masterData[anm].keys():
            print(f"Loading {field} data for {anm}...")
            if not field == 'shiftInfo':
                if field == 'behavior':
                    masterData[anm][field] = unpickler(os.path.join(DATA_DIR,anm,field),anm+"_"+field,extension = '.pkl')
                elif field == 'somas' or field == 'dendrites':
                    if isinstance(masterData[anm][field], dict):
                        temp = unpickler(os.path.join(DATA_DIR,anm,field),anm+"_"+field,extension = '.pkl')
                        masterData[anm][field]['daysRelativeToShift'] = copy.deepcopy(temp['daysRelativeToShift'])
                        masterData[anm][field]['roiScores'] = copy.deepcopy(temp['roiScores'])
                        del temp
                        for sess in masterData[anm][field]['sessions'].keys():
                            print(f"Loading session {sess} data for {anm}...")
                            masterData[anm][field]['sessions'][sess] = unpickler(os.path.join(DATA_DIR,anm,field),anm+"_"+field+"_session_"+str(sess),extension = '.pkl')


    return masterData
##################################################################################################
def load_ROI_clustering_reduced_split(savePath,saveName,ROI_type_keys=None,ROI_scores=None,big_keys=None,clusterList=None,trialGroupings=None,extension='.pkl',verbose=True):
    """
    Rebuild a ROI_clustering_reduced dict written by save_ROI_clustering_reduced_split, along with the
    cueInfo dict stored beside it. The metadata file '<saveName>_meta<extension>' is always read (it
    carries every cluster's 'nROIs'/'ROI_info'/'perROIs', the '_split_index' listing every part file,
    plus '_cueInfo' if one was saved); only the part files matching ROI_type_keys / ROI_scores /
    big_keys / clusterList / trialGroupings are then loaded and merged back in, so you can pull just
    one ROI type, a couple of clusters, or one trialGrouping without reading the whole multi-GB
    dataset. When the data was saved with split_trialGroupings=True, a trialGroupings filter skips the
    other trialGroupings' files entirely; otherwise the file is read once and then pruned in memory.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    savePath = os.path.join(DATA_DIR,'traceAlignments')  # directory holding the '_meta' and part files
    saveName = 'fig3_clustered_data'   # same base name passed to save_ROI_clustering_reduced_split
    ROI_type_keys = None               # None -> all ROI types; or list, e.g. ['dendrites','somas']
    ROI_scores = None                  # None -> all scores; or list, e.g. [1]
    big_keys = None                    # None -> every big key that was saved; or list, e.g. ['clusters']; [] -> metadata only
    clusterList = None                 # None/[] -> all clusters; or list of cluster indices to keep
    trialGroupings = ['byPrePost_byTtype','byPrePost_deltaTtypes']  # None/[] -> all saved trialGroupings; or str/list of trialGroupings to keep
    extension = '.pkl'                 # '.pkl' | '.p'  (extension used when the files were written)
    verbose = True                     # True | False  (prints each part file loaded)
    importlib.reload(jsm)
    ROI_clustering_reduced, cueInfo = load_ROI_clustering_reduced_split(
        savePath, saveName, ROI_type_keys=ROI_type_keys, ROI_scores=ROI_scores, big_keys=big_keys, clusterList=clusterList, trialGroupings=trialGroupings, extension=extension, verbose=verbose)
    # cueInfo comes back as None if the data was saved without one; then load it the old way:
    # cueInfo = safe_unpickler(summaryDataDir,"cueInfo",False)

    Parameters
    ----------
    savePath : str
        Directory holding the metadata and part files.
    saveName : str
        Base file name used when saving. Any '.pkl'/'.p' extension on it is stripped.
    ROI_type_keys : list of str or None
        ROI types to load. None keeps every ROI type present in the metadata file.
    ROI_scores : list or None
        ROI scores to load. None keeps every score present in the metadata file.
    big_keys : list of str or None
        Which of the saved big keys to load back in. None loads all of them; pass [] to get the
        metadata-only dict with no trace data.
    clusterList : list of int or None
        Cluster indices to keep. None or [] keeps all clusters. When the data was saved with
        split_clusters=True only these clusters' files are read; otherwise the whole 'clusters'
        file is read and then pruned to these indices.
    trialGroupings : str, list of str, or None
        Trial grouping(s) to keep, e.g. 'byPrePost_byTtype' or ['byPrePost_byTtype',
        'byPrePost_deltaTtypes']. None or [] keeps every trialGrouping in the saved data. When the
        data was saved with split_trialGroupings=True only the matching files are read; otherwise
        each file is read and then pruned to these trialGroupings.
    extension : str
        File extension used when the files were written, e.g. '.pkl'.
    verbose : bool
        If True, print each file loaded and a summary of what was reassembled.

    Returns
    -------
    ROI_clustering_reduced : dict
        Dict keyed [ROI_type_key][ROI_score][...] with the requested metadata and big keys merged
        back in, in the same layout reduce_ROI_clustering_for_raster_fig produced.
    cueInfo : dict or None
        The cueInfo dict stored in the metadata file, or None if the data was saved without one.
        Returned separately (rather than left in ROI_clustering_reduced) so the top level of
        ROI_clustering_reduced stays purely ROI_type_key entries.
    """
    import pickle

    base = saveName
    if base.lower().endswith('.pkl') or base.lower().endswith('.p'):
        base = os.path.splitext(base)[0]

    def _read(fileName):
        full_path = os.path.join(savePath,fileName)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"load_ROI_clustering_reduced_split: missing file {full_path}")
        with open(full_path,'rb') as f:
            obj = pickle.load(f)
        if verbose:
            print(f"  + loaded {fileName}  ({np.round(os.path.getsize(full_path)/1e9,decimals=3)} GB)")
        return obj

    if isinstance(trialGroupings,str):
        trialGroupings = [trialGroupings]
    if not trialGroupings:
        trialGroupings = None

    if verbose:
        print(f"load_ROI_clustering_reduced_split: savePath={savePath!r} saveName={base!r} "
              f"ROI_type_keys={ROI_type_keys} ROI_scores={ROI_scores} big_keys={big_keys} "
              f"clusterList={clusterList} trialGroupings={trialGroupings}")

    ROI_clustering_reduced = _read(f"{base}_meta{extension}")
    split_index = ROI_clustering_reduced.pop('_split_index',{'files':[]})
    cueInfo = ROI_clustering_reduced.pop('_cueInfo',None)
    if verbose:
        print(f"  = cueInfo: {'loaded from metadata file' if cueInfo is not None else 'not stored with this dataset'}")

    for ROI_type_key in list(ROI_clustering_reduced.keys()):
        if ROI_type_keys is not None and ROI_type_key not in ROI_type_keys:
            del ROI_clustering_reduced[ROI_type_key]
            continue
        for ROI_score in list(ROI_clustering_reduced[ROI_type_key].keys()):
            if ROI_scores is not None and ROI_score not in ROI_scores:
                del ROI_clustering_reduced[ROI_type_key][ROI_score]
        if not ROI_clustering_reduced[ROI_type_key]:
            del ROI_clustering_reduced[ROI_type_key]

    # prune the per-cluster descriptors that came from the metadata file down to clusterList
    if clusterList:
        for ROI_type_key in ROI_clustering_reduced:
            for ROI_score in ROI_clustering_reduced[ROI_type_key]:
                clusters_meta = ROI_clustering_reduced[ROI_type_key][ROI_score].get('clusters')
                if clusters_meta is not None:
                    for c in list(clusters_meta.keys()):
                        if c not in clusterList:
                            del clusters_meta[c]

    for entry in split_index['files']:
        ROI_type_key = entry['ROI_type_key']
        ROI_score = entry['ROI_score']
        if ROI_type_key not in ROI_clustering_reduced or ROI_score not in ROI_clustering_reduced[ROI_type_key]:
            continue
        if big_keys is not None and entry['key'] not in big_keys:
            continue
        if entry['cluster'] is not None and clusterList and entry['cluster'] not in clusterList:
            continue
        entry_tg = entry.get('trialGrouping')
        if entry_tg is not None and trialGroupings is not None and entry_tg not in trialGroupings:
            continue
        dst = ROI_clustering_reduced[ROI_type_key][ROI_score]
        data = _read(entry['file'])

        if entry['key'] != 'clusters':
            dst[entry['key']] = data
            continue

        # part files hold trace branches only; merge them onto the cluster descriptors from the
        # metadata file, accumulating across per-cluster / per-trialGrouping files
        prune_tg = trialGroupings if entry_tg is None else None
        clusters_dst = dst.setdefault('clusters',{})
        if entry['cluster'] is None:
            for c,cluster_part in data.items():
                if clusterList and c not in clusterList:
                    continue
                if prune_tg is not None:
                    cluster_part = _prune_cluster_trialGroupings(cluster_part,trialGroupings=prune_tg)
                if cluster_part:
                    _deep_merge(clusters_dst.setdefault(c,{}),cluster_part)
        else:
            cluster_part = data
            if prune_tg is not None:
                cluster_part = _prune_cluster_trialGroupings(cluster_part,trialGroupings=prune_tg)
            if cluster_part:
                _deep_merge(clusters_dst.setdefault(entry['cluster'],{}),cluster_part)

    # the raster figs read ['clusters'][c]['nROIs'] / ['ROI_info'] straight away, so flag a cluster
    # that came back without them instead of letting it fail deep inside the figure code
    for ROI_type_key in ROI_clustering_reduced:
        for ROI_score in ROI_clustering_reduced[ROI_type_key]:
            clusters_dst = ROI_clustering_reduced[ROI_type_key][ROI_score].get('clusters',{})
            missing = [c for c in clusters_dst if 'nROIs' not in clusters_dst[c]]
            if missing:
                print(f"<<WARNING>> [{ROI_type_key!r}][{ROI_score!r}]: clusters {missing} loaded without "
                      f"'nROIs'/'ROI_info' — the metadata file and the part files were written by "
                      f"different versions of save_ROI_clustering_reduced_split; re-save from the "
                      f"in-memory ROI_clustering_reduced to fix")

    if verbose:
        for ROI_type_key in ROI_clustering_reduced:
            for ROI_score in ROI_clustering_reduced[ROI_type_key]:
                clusters_dst = ROI_clustering_reduced[ROI_type_key][ROI_score].get('clusters',{})
                print(f"  = [{ROI_type_key!r}][{ROI_score!r}]: "
                      f"{len(ROI_clustering_reduced[ROI_type_key][ROI_score])} keys, {len(clusters_dst)} clusters, "
                      f"trialGroupings {_list_trialGroupings(clusters_dst)}")

    return ROI_clustering_reduced,cueInfo
# small per-cluster descriptor keys — everything else in a cluster dict is a [group][align_data][trialGrouping]
# trace branch. Used by save_ROI_clustering_reduced_split / load_ROI_clustering_reduced_split to decide what
# belongs in the metadata file vs the bulky part files.
_CLUSTER_META_KEYS = ['nROIs','ROI_info','perROIs']

def _cluster_trace_only(cluster_src,trialGroupings=None):
    """Return the [group][align_data][trialGrouping] trace branches of one cluster dict (i.e. everything
    outside _CLUSTER_META_KEYS), optionally pruned to trialGroupings. Returns {} if nothing matches."""
    trace_part = {}
    for group in cluster_src:
        if group in _CLUSTER_META_KEYS:
            continue
        for align_data in cluster_src[group]:
            for tg in cluster_src[group][align_data]:
                if trialGroupings is not None and tg not in trialGroupings:
                    continue
                trace_part.setdefault(group,{}).setdefault(align_data,{})[tg] = \
                    copy.deepcopy(cluster_src[group][align_data][tg])
    return trace_part

def _prune_cluster_trialGroupings(cluster_src,trialGroupings):
    """Prune one cluster dict's trace branches to trialGroupings while KEEPING its _CLUSTER_META_KEYS
    descriptors. Used on load, where a part file may still carry 'nROIs'/'ROI_info'/'perROIs' (part
    files written before those moved into the metadata file) — dropping them there would leave the
    raster figs with a cluster that has traces but no nROIs."""
    pruned = {k:v for k,v in cluster_src.items() if k in _CLUSTER_META_KEYS}
    pruned.update(_cluster_trace_only(cluster_src,trialGroupings=trialGroupings))
    return pruned

def _list_trialGroupings(clusters_src):
    """Every trialGrouping key present anywhere in a clusters dict, in first-seen order."""
    trialGroupings = []
    for c in clusters_src:
        for group in clusters_src[c]:
            if group in _CLUSTER_META_KEYS:
                continue
            for align_data in clusters_src[c][group]:
                for tg in clusters_src[c][group][align_data]:
                    if tg not in trialGroupings:
                        trialGroupings.append(tg)
    return trialGroupings

def _deep_merge(dst,src):
    """Recursively merge src into dst, so trace branches loaded from separate part files
    (different clusters, different trialGroupings) accumulate instead of overwriting."""
    for key,value in src.items():
        if isinstance(value,dict) and isinstance(dst.get(key),dict):
            _deep_merge(dst[key],value)
        else:
            dst[key] = value
    return dst

#########################################################################################################
# Loader for the compressed all_grouped_aligned_traces slices written by
# jsm_figs.save_all_grouped_aligned_traces_split — the individual-example
# trace data for the *_indiv_examples / *_indiv_pre_post_examples notebooks, keyed
#   [group]['grouped_aligned_traces'][align_data]['aligned_traces'][anm][ROI_type_key][ROI_score][ROI_idx]
#       [trace_type][<behavior_epoch | sess | bsess>][ttype]
# The by-ROI figs read the session list (and [0][ttype]['align_shift_fr']) out of the 'bySess_byTtype'
# branch for EVERY trace_type, so that branch is always kept whatever trace_types filter is given.
_ALIGN_SESS_TRACE_TYPE = 'bySess_byTtype'
# trace_types that are DISPLAY MODES rather than stored branches: the aligner never writes a
# 'preShiftOnly_byTtype' / 'postShiftOnly_byTtype' key, and trial_alignments_byROI_pages immediately
# rebinds trace_type to 'byPrePost_byTtype' and narrows behavior_epochs to one side. Anything that filters
# aligned_traces by trace_type has to resolve them the same way or it prunes away the branch the figs read.
_ALIGN_TRACE_TYPE_MODES = {'preShiftOnly_byTtype':('byPrePost_byTtype',['preShift']),
                           'postShiftOnly_byTtype':('byPrePost_byTtype',['postShift'])}

def resolve_align_trace_type(trace_type):
    """Map a figure trace_type onto the branch aligned_traces stores it under, plus the behavior epochs
    that trace_type restricts itself to. Returns (stored_trace_type, mode_epochs or None); see
    jsm_figs.resolve_align_trace_type."""
    return _ALIGN_TRACE_TYPE_MODES.get(trace_type,(trace_type,None))

def _align_restore_floats(node,float_dtype):
    """In-place upcast of every reduced-precision float ndarray (e.g. the float16 arrays written by
    save_all_grouped_aligned_traces_split) back to float_dtype. Arrays already at or above that precision
    are left alone, so a float64 array is never silently downcast on load."""
    if float_dtype is None:
        return node
    target = np.dtype(float_dtype)
    if isinstance(node,dict):
        for key,value in node.items():
            if isinstance(value,np.ndarray):
                if np.issubdtype(value.dtype,np.floating) and value.dtype.itemsize < target.itemsize:
                    node[key] = value.astype(target)
            elif isinstance(value,dict):
                _align_restore_floats(value,float_dtype)
    return node

def _align_open(full_path,mode,extension,compresslevel=9):
    """Open a part file with the compressor implied by extension: '.gz' -> gzip, '.xz'/'.lzma' -> lzma,
    '.bz2' -> bz2, anything else -> plain binary. compresslevel/preset is only passed when writing."""
    low = extension.lower()
    writing = 'w' in mode
    if low.endswith('.gz'):
        import gzip
        return gzip.open(full_path,mode,compresslevel=compresslevel) if writing else gzip.open(full_path,mode)
    if low.endswith('.xz') or low.endswith('.lzma'):
        import lzma
        return lzma.open(full_path,mode,preset=compresslevel) if writing else lzma.open(full_path,mode)
    if low.endswith('.bz2'):
        import bz2
        return bz2.open(full_path,mode,compresslevel=compresslevel) if writing else bz2.open(full_path,mode)
    return open(full_path,mode)

def load_all_grouped_aligned_traces_split(savePath,saveName,groups=None,align_datas=None,anms=None,ROI_type_keys=None,ROI_scores=None,ROI_idxs=None,trace_types=None,extension='.pkl.gz',restore_float_dtype='float32',verbose=True):
    """
    Rebuild an all_grouped_aligned_traces dict written by save_all_grouped_aligned_traces_split, along with
    the cueInfo and summaryInfo stored beside it. The metadata file '<saveName>_meta<extension>' is always
    read (it carries the per-group descriptors, each align_data's traceAlignParams, the '_split_index'
    listing every part file, plus '_cueInfo' / '_summaryInfo' if they were saved); only the part files
    matching groups / align_datas / anms / ROI_type_keys / ROI_scores / ROI_idxs are then read and merged
    back in, so a notebook can pull just the one example ROI its panel needs. The result is keyed exactly
    like the full all_grouped_aligned_traces, so
    all_grouped_aligned_traces[group]['grouped_aligned_traces'][align_data]['aligned_traces'] can be handed
    straight to trial_alignments_byROI_pages.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    savePath = os.path.join(DATA_DIR,'traceAlignments','fig6')  # directory holding the '_meta' and part files
    saveName = 'fig6_indiv_examples'   # same base name passed to save_all_grouped_aligned_traces_split
    groups = None                      # None -> every group saved; or list, e.g. ['goCue']
    align_datas = None                 # None -> every align_data saved; or list, e.g. ['allTC_dffbc_decon_sp_events_MatchDx2Sx1']
    anms = None                        # None -> every animal saved; or list, e.g. ['B00002213982']
    ROI_type_keys = None               # None -> all ROI types; or list, e.g. ['dendrites','somas']
    ROI_scores = None                  # None -> all scores; or list, e.g. [1]
    ROI_idxs = None                    # None -> all example ROIs; or list of ROI indices to keep
    trace_types = None                 # None -> every trace_type saved; or list, e.g. ['byPrePost_byTtype'] ('preShiftOnly_byTtype' resolves to it)
    extension = '.pkl.gz'              # '.pkl.gz' | '.pkl.xz' | '.pkl.bz2' | '.pkl'  (extension used when the files were written)
    restore_float_dtype = 'float32'    # None -> keep the stored dtype; or 'float32' to cast float16 traces back up
    verbose = True                     # True | False  (prints each part file loaded)
    all_grouped_aligned_traces, cueInfo, summaryInfo = load_all_grouped_aligned_traces_split(
        savePath, saveName, groups=groups, align_datas=align_datas, anms=anms, ROI_type_keys=ROI_type_keys, ROI_scores=ROI_scores, ROI_idxs=ROI_idxs, trace_types=trace_types, extension=extension, restore_float_dtype=restore_float_dtype, verbose=verbose)
    # cueInfo / summaryInfo come back as None if the data was saved without them; then load them the old way:
    # cueInfo = safe_unpickler(summaryDataDir,"cueInfo",False)

    Parameters
    ----------
    savePath : str
        Directory holding the metadata and part files.
    saveName : str
        Base file name used when saving. Any '.pkl'/'.p' extension (with or without a compression suffix)
        is stripped.
    groups : list of str or None
        Alignment groups to load, e.g. ['goCue']. None keeps every group in the metadata file.
    align_datas : list of str or None
        Align-data keys to load. None keeps every align_data in the metadata file.
    anms : list of str or None
        Animals to load. None keeps every animal that was saved.
    ROI_type_keys : list of str or None
        ROI types to load. None keeps every ROI type that was saved.
    ROI_scores : list or None
        ROI scores to load. None keeps every score that was saved.
    ROI_idxs : list of int or None
        Example ROI indices to keep. None keeps them all. When the data was saved with split_ROIs=True
        only these ROIs' files are read; otherwise the file is read and then pruned to these indices.
    trace_types : list of str or None
        Trial grouping(s) to keep, e.g. ['byPrePost_byTtype']. None keeps every trace_type in the saved
        data. Display modes ('preShiftOnly_byTtype'/'postShiftOnly_byTtype') are resolved to the stored
        'byPrePost_byTtype' branch they are a view of, so filtering by one does not prune away the data.
        The 'bySess_byTtype' support branch is always kept, since the figs read the session list from it
        for every trace_type.
    extension : str
        File extension used when the files were written; the decompressor is inferred from it.
    restore_float_dtype : str or None
        dtype to cast every float array back to after loading, e.g. 'float32' to undo a float16 save.
        None leaves the arrays at their stored precision.
    verbose : bool
        If True, print each file loaded and a summary of what was reassembled.

    Returns
    -------
    all_grouped_aligned_traces : dict
        Dict keyed [group]['grouped_aligned_traces'][align_data]['traceAlignParams' | 'aligned_traces'],
        in the same layout as the full dict the figures normally run on.
    cueInfo : dict or None
        The cueInfo dict stored in the metadata file, or None if the data was saved without one.
    summaryInfo : dict or None
        The summaryInfo dict stored in the metadata file, or None if the data was saved without one.
    """
    import pickle

    base = saveName
    for suffix in ['.gz','.xz','.lzma','.bz2']:
        if base.lower().endswith(suffix):
            base = base[:-len(suffix)]
    if base.lower().endswith('.pkl') or base.lower().endswith('.p'):
        base = os.path.splitext(base)[0]

    def _read(fileName):
        full_path = os.path.join(savePath,fileName)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"load_all_grouped_aligned_traces_split: missing file {full_path}")
        with _align_open(full_path,'rb',extension) as f:
            obj = pickle.load(f)
        if verbose:
            print(f"  + loaded {fileName}  ({np.round(os.path.getsize(full_path)/1e6,decimals=2)} MB)")
        return obj

    if isinstance(groups,str):
        groups = [groups]
    if isinstance(align_datas,str):
        align_datas = [align_datas]
    if isinstance(anms,str):
        anms = [anms]
    if isinstance(ROI_type_keys,str):
        ROI_type_keys = [ROI_type_keys]
    if isinstance(trace_types,str):
        trace_types = [trace_types]
    if isinstance(ROI_idxs,(int,np.integer)):
        ROI_idxs = [ROI_idxs]
    if not trace_types:
        trace_types = None
    if trace_types is not None:
        # a filter of 'preShiftOnly_byTtype' has to keep the 'byPrePost_byTtype' branch it is a view of,
        # otherwise it prunes away the only data the figs actually read
        trace_types = list(dict.fromkeys(resolve_align_trace_type(trace_type)[0] for trace_type in trace_types))

    if verbose:
        print(f"load_all_grouped_aligned_traces_split: savePath={savePath!r} saveName={base!r} "
              f"groups={groups} align_datas={align_datas} anms={anms} ROI_type_keys={ROI_type_keys} "
              f"ROI_scores={ROI_scores} ROI_idxs={ROI_idxs} trace_types={trace_types}")

    all_grouped_aligned_traces = _read(f"{base}_meta{extension}")
    split_index = all_grouped_aligned_traces.pop('_split_index',{'files':[]})
    cueInfo = all_grouped_aligned_traces.pop('_cueInfo',None)
    summaryInfo = all_grouped_aligned_traces.pop('_summaryInfo',None)
    if verbose:
        print(f"  = cueInfo: {'loaded from metadata file' if cueInfo is not None else 'not stored with this dataset'}"
              f" | summaryInfo: {'loaded from metadata file' if summaryInfo is not None else 'not stored with this dataset'}")

    for group in list(all_grouped_aligned_traces.keys()):
        if groups is not None and group not in groups:
            del all_grouped_aligned_traces[group]
            continue
        for align_data in list(all_grouped_aligned_traces[group].get('grouped_aligned_traces',{}).keys()):
            if align_datas is not None and align_data not in align_datas:
                del all_grouped_aligned_traces[group]['grouped_aligned_traces'][align_data]
        if not all_grouped_aligned_traces[group].get('grouped_aligned_traces',{}):
            del all_grouped_aligned_traces[group]

    for entry in split_index['files']:
        group = entry['group']
        align_data = entry['align_data']
        if group not in all_grouped_aligned_traces:
            continue
        if align_data not in all_grouped_aligned_traces[group].get('grouped_aligned_traces',{}):
            continue
        if anms is not None and entry['anm'] not in anms:
            continue
        if ROI_type_keys is not None and entry['ROI_type_key'] not in ROI_type_keys:
            continue
        if ROI_scores is not None and entry['ROI_score'] not in ROI_scores:
            continue
        if entry['ROI_idx'] is not None and ROI_idxs and entry['ROI_idx'] not in ROI_idxs:
            continue
        payload = _read(entry['file'])

        # part files always hold {ROI_idx: {trace_type: branch}}, whether or not they were split per ROI
        dst_traces = all_grouped_aligned_traces[group]['grouped_aligned_traces'][align_data].setdefault('aligned_traces',{})
        dst_ROIs = dst_traces.setdefault(entry['anm'],{}).setdefault(entry['ROI_type_key'],{}).setdefault(entry['ROI_score'],{})
        for ROI_idx,ROI_src in payload.items():
            if ROI_idxs and ROI_idx not in ROI_idxs:
                continue
            if trace_types is not None:
                ROI_src = {trace_type:branch for trace_type,branch in ROI_src.items()
                           if trace_type in trace_types or trace_type == _ALIGN_SESS_TRACE_TYPE}
            if restore_float_dtype is not None:
                for branch in ROI_src.values():
                    _align_restore_floats(branch,restore_float_dtype)
            _deep_merge(dst_ROIs.setdefault(ROI_idx,{}),ROI_src)

    if verbose:
        for group in all_grouped_aligned_traces:
            for align_data in all_grouped_aligned_traces[group]['grouped_aligned_traces']:
                dst_traces = all_grouped_aligned_traces[group]['grouped_aligned_traces'][align_data].get('aligned_traces',{})
                for anm in dst_traces:
                    for ROI_type_key in dst_traces[anm]:
                        for ROI_score in dst_traces[anm][ROI_type_key]:
                            for ROI_idx in dst_traces[anm][ROI_type_key][ROI_score]:
                                print(f"  = [{group!r}][{align_data!r}][{anm!r}][{ROI_type_key!r}][{ROI_score!r}]"
                                      f"[{ROI_idx!r}]: trace_types "
                                      f"{list(dst_traces[anm][ROI_type_key][ROI_score][ROI_idx].keys())}")

    return all_grouped_aligned_traces,cueInfo,summaryInfo

##################################################################################################    
def summary_stats(tempData,nBins=100,binWidth=0,temp_minVal=np.nan,temp_maxVal=np.nan,\
    calcHist=True,verbose=False, allow_overages = False, splitZeroBin = False, \
    calc_hist_log10 = True, calc_cumulative_hist = True, calc_norm_hist = True, calc_cumulative_freq = True, calc_pdf = True, 
    save_input = False, force_dType=False, dType='float32', exclude_NaNs = True, ddof = 0, warningsOn = True):
    """
    Calculate summary statistics for the provided data.

    Parameters:
    tempData: input data to analyze
    nBins: number of bins for histogram
    binWidth: width of each bin for histogram
    temp_minVal: minimum value for histogram
    temp_maxVal: maximum value for histogram
    calcHist: if True, calculate histogram and related statistics
    verbose: if True, print additional information
    allow_overages: if True, allow values greater than the histogram limits
    splitZeroBin: if True, split the zero bin in the histogram
    calc_hist_log10: if True, calculate logarithmic histogram
    calc_cumulative_hist: if True, calculate cumulative histogram
    calc_norm_hist: if True, calculate normalized histogram
    calc_cumulative_freq: if True, calculate cumulative frequency
    calc_pdf: if True, calculate probability density function
    save_input: if True, save the input data in the stats dictionary
    force_dType: if True, force the data to a specific data type
    dType: data type to force the data to
    exclude_NaNs: if True, exclude NaN values from the calculations
    warningsOn: if True, print warnings

    Returns:
    stats: a dictionary containing the calculated statistics
    """

    stats={}
    tempData = np.array(tempData)
    tempData=tempData.flatten()
    if len(tempData) == 0:
        if warningsOn:
            print("<<WARNING>> No data to analyze")
        stats['success'] = False
        stats['n'] = np.nan
        stats['mean'] = np.nan
        stats['gmean'] = np.nan
        stats['max'] = np.nan
        stats['min'] = np.nan
        stats['med'] = np.nan
        stats['std'] = np.nan
        stats['sem'] = np.nan
        stats['cv'] = np.nan
        stats['sum'] = np.nan
        stats['mode']=np.nan
        stats['q1']=np.nan
        stats['q3']=np.nan
        stats['iqr']=np.nan
        stats['np_skew']=np.nan
        stats['fmc_skew']=np.nan
    else:
        if force_dType:
            try:
                tempData = np.array(tempData).astype(dType)
            except:
                pass
        ##################
        #Basic Stats
        if save_input:
            stats['data'] = copy.deepcopy(tempData)
        stats['success'] = False
        if exclude_NaNs:
            stats['n']=sum(~np.isnan(tempData))
            stats['mean']=np.nanmean(tempData)
            stats['max']=np.nanmax(tempData)
            stats['min']=np.nanmin(tempData)
            try:
                stats['minPos']=np.nanmin(tempData[tempData>0])
            except:
                stats['minPos']=np.nan
            stats['med']=np.nanmedian(tempData,axis=0)
            stats['std']=np.nanstd(tempData,ddof=ddof)
            stats['sum']=np.nansum(tempData)
            try:
                stats['q1']=np.nanpercentile(tempData,25)
                stats['q3']=np.nanpercentile(tempData,75)
            except:
                stats['q1']=np.nan
                stats['q3']=np.nan
            # stats['mode']=mode(tempData, nan_policy='omit')
            stats['gmean']=gmean(tempData[np.isfinite(tempData)])

        else:
            stats['n']=sum(tempData)
            stats['mean']=np.mean(tempData)
            stats['max']=np.max(tempData)
            stats['min']=np.min(tempData)
            try:
                stats['minPos']=np.min(tempData[tempData>0])
            except:
                stats['minPos']=np.nan
            stats['med']=np.median(tempData,axis=0)
            stats['std']=np.std(tempData,ddof=ddof)
            stats['sum']=np.sum(tempData)
            stats['mode']=mode(tempData)
            try:
                stats['q1']=np.percentile(tempData,25)
                stats['q3']=np.percentile(tempData,75)
            except:
                stats['q1']=np.nan
                stats['q3']=np.nan
            # stats['mode']=mode(tempData, nan_policy='propogate')
            stats['gmean']=gmean(tempData)
        try:
            stats['iqr']=stats['q3']-stats['q1']
        except:
            stats['iqr']=np.nan
        try:
            stats['sem']=stats['std']/np.sqrt(stats['n'])
        except:
            stats['sem']=np.nan
        try:
            stats['cv']=stats['std']/stats['mean']
        except:
            stats['cv']=np.nan
        try:
            stats['np_skew']=(stats['mean']-stats['med'])/stats['std'] # nonparametric skew
        except:
            stats['np_skew']=np.nan
        try:
            if exclude_NaNs:
                stats['fmc_skew']=float(skew(tempData, axis=0, bias=True, nan_policy='omit')) #Fisher's moment coefficient of skewness
            else:
                stats['fmc_skew']=float(skew(tempData, axis=0, bias=True, nan_policy='propogate')) #Fisher's moment coefficient of skewness

        except:
            stats['fmc_skew'] = np.nan
        ##################
        if calcHist:
            stats['binWidth']=binWidth
            stats['minVal']=temp_minVal
            stats['maxVal']=temp_maxVal
            stats['nBins']=nBins
            try:
                # try:
                if nBins and ('float' in str(type(stats['minVal'])) or 'int' in str(type(stats['minVal']))) and \
                    ('float' in str(type(stats['maxVal'])) or 'int' in str(type(stats['maxVal']))):
                    binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                elif binWidth>0 and ('float' in str(type(stats['minVal'])) or 'int' in str(type(stats['minVal']))) and \
                    ('float' in str(type(stats['maxVal'])) or 'int' in str(type(stats['maxVal']))):
                    binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                elif nBins:
                    if exclude_NaNs:
                        stats['minVal'] = np.floor(np.nanmin(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.nanmax(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                    else:
                        stats['minVal'] = np.floor(np.min(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.max(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                elif binWidth>0:
                    if exclude_NaNs:
                        stats['minVal'] = np.floor(np.nanmin(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.nanmax(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                    else:
                        stats['minVal'] = np.floor(np.min(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.max(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                binCenters=binEdges[0:len(binEdges)-1]+(binEdges[1]-binEdges[0])/2
                zeroMin = False
                if exclude_NaNs:
                    if np.nanmin(tempData) == 0 and splitZeroBin:
                        zeroMin = True
                else:
                    if np.min(tempData) == 0 and splitZeroBin:
                        zeroMin = True
                if zeroMin:
                    if verbose:
                        print("NOTE: adding true zero bin")
                    binEdges[0] = 0.0000001
                    binEdges = np.insert(binEdges,0,0)
                    binCenters = np.insert(binCenters,0,0)

                hist,BinCenters=np.histogram(tempData,bins=binEdges)
                # except:
                #     hist,BinCenters=np.histogram(tempData)
                if allow_overages and np.any(tempData>binEdges[-1]):
                    stats['overage'] = True
                    if verbose:
                        print("NOTE: adding an extra bin for values greater than the limits")
                    overage = np.sum(tempData>binEdges[-1])
                    hist = np.append(hist,overage)
                    binCenters = np.append(binCenters,binCenters[-1]+(binCenters[-1]-binCenters[-3]))
                    if exclude_NaNs:
                        binEdges = np.append(binEdges,np.nanmax(tempData))
                    else:
                        binEdges = np.append(binEdges,np.max(tempData))

                hist_log10 = copy.deepcopy(hist)
                # hist_log10[hist_log10<1] = 0.1
                hist_log10 = np.log10(hist_log10)
                hist_log10[np.isinf(hist_log10)] = 0
                norm_hist=hist/np.nanmax(hist)
                cumulative_hist=np.cumsum(hist)
                pdf = hist/stats['n']
                if np.nanmax(cumulative_hist) == np.nan:
                    cumulative_freq=cumulative_hist
                else:
                    cumulative_freq=cumulative_hist/np.nanmax(cumulative_hist)
                stats['binEdges']=binEdges.astype('float32')
                stats['binCenters']=binCenters.astype('float32')
                stats['hist']=hist.astype('int32')
                if calc_hist_log10:
                    stats['hist_log10']=hist_log10.astype('float32')
                if calc_cumulative_hist:
                    stats['cumulative_hist']=cumulative_hist.astype('float32')
                if calc_norm_hist:
                    stats['norm_hist']=norm_hist.astype('float32')
                if calc_cumulative_freq:
                    stats['cumulative_freq']=cumulative_freq.astype('float32')
                if calc_pdf:
                    stats['pdf']=pdf.astype('float32')
                stats['success'] = True
                if not stats['nBins']:
                    stats['nBins'] = len(hist)
            except:
                if warningsOn:
                    print("<<WARNING>> Unable to calculate histograms...")
                # stats['hist'] = [np.nan]
                # stats['cumulative_hist'] = [np.nan]
                # stats['norm_hist'] = [np.nan]
                # stats['cumulative_freq'] = [np.nan]
        else:
            if verbose:
                print("Skipping Histograms")
    ##################
    if verbose:
        print("n    = "+str(stats['n']))
        print("mean = "+str(stats['mean']))
        print("max  = "+str(stats['max']))
        print("min  = "+str(stats['min']))
        print("med  = "+str(stats['med']))
        print("std  = "+str(stats['std']))
        print("sem  = "+str(stats['sem']))
        print("cv   = "+str(stats['cv']))
        print("sum  = "+str(stats['sum']))
    return stats
##################################################################################################

def generate_default_flagParams(maskKey = 'consensus_NMFtc_mask'):
    flagParams = {}
    flagParams['traceKey_mask'] = maskKey               # This is the source of frame-specific missing data
    ###########################
    flagParams['byTrial'] = {}
    flagParams['byTrial']['min_trial_frs'] = 10                        # Just to filter out any trials that are super short regardless of why
    flagParams['byTrial']['use_missing_percent_trial_frs'] = True      # Toggle to filter trials missing > % Frs
    flagParams['byTrial']['max_missing_percent_trial_frs'] = 10        # Flags any trials missing > % Frs
    ###########################
    flagParams['bySess'] = {}
    flagParams['bySess']['use_missing_percent_frs'] = True             # Toggle to filter sessions missing > % Frs
    flagParams['bySess']['max_missing_percent_frs'] = 10               # Threshold > % whole session frames missing
    ###########
    flagParams['bySess']['use_missing_percent_flagged_trials'] = True  # Toggle to filter sessions with > % Trials being flagged using max_missing_percent_trial_frs and min_trial_frs
    flagParams['bySess']['max_missing_percent_flagged_trials'] = 10    # Threshold > % Trials being flagged using max_missing_percent_trial_frs and min_trial_frs
    ###########
    flagParams['bySess']['use_missing_num_stretches'] = True           # Toggle to filter sessions with > # contig-missing-frame-stretches
    flagParams['bySess']['max_missing_num_stretches'] = 120            # Threshold  > # contig-missing-frame-stretches
    ###########
    flagParams['bySess']['use_missing_stretch_len_s_stats'] = True     # Toggle to check the contig-missing-frame-stretches length stats
    flagParams['bySess']['missing_stretch_len_s_hist_binWidth'] = 1    # Histogram config in seconds for contig-missing-frame-stretches bin width
    flagParams['bySess']['missing_stretch_len_s_hist_binMax'] = 10     # Histogram config with binWidth of 1s this will make 10, 1s bins, 0-1s, 1-2s, .... , 8-9s,9-inf
    flagParams['bySess']['max_missing_num_short_stretches'] = 100      # Threshold  > # contig-missing-frame-stretches in the smallest bin above
    flagParams['bySess']['max_missing_num_med_stretches'] = 30         # Threshold  > # contig-missing-frame-stretches in all middle bins
    flagParams['bySess']['max_missing_num_long_stretches'] = 10        # Threshold  > # contig-missing-frame-stretches in just largest bin including anything greater than binMax
    ###########################
    return flagParams

def get_tsepMatrix(tseps):
    minIdxs = []
    badIdxs = [[np.nan, np.nan],[np.nan, np.nan]]
    for d in range(len(tseps)):
        dayTseps = tseps[d]
        for t in range(len(dayTseps)):
            tTsep = dayTseps[t]
            if tTsep.shape[1]!=0:
                minIdxs.append(tTsep.shape[1])
            else:
                badIdxs.append([d,t])
    minIdx = np.min(minIdxs)
    badIdxs = np.array(badIdxs)

    nRois = tseps[0][0].shape[0]
    fillerTrial = np.ones([nRois,minIdx])*np.nan
    tsep = []
    for d in range(len(tseps)):
        dFlag = False
        dayTseps = tseps[d]
        if d in badIdxs[:,0]:
            dFlag = True
        for t in range(len(dayTseps)):
            if dFlag and t in badIdxs[:,1]:
                tsep.append(fillerTrial)
            else:
                tsep.append(tseps[d][t][:,:minIdx])
    tsep = np.array(tsep)
    return tsep
      
def getGoodRois(anatData, keepTrsh = 1, switch = False, dffKey = 'dff', useScores = False, goodROI_scores = [1], verbose = True):
    if switch:
        keepTrsh = -keepTrsh
    if 'NMF' in dffKey:
        try:
            roiScores = anatData['roiScores']['NMF']
        except:
            print('<<< dffKey = ' +str(dffKey)+' error in getGoodRois()>>>')
            print('either you asked for NMF somas which arent available yet or something else... bug jackson abt it')
    elif 'sp' in dffKey or 'impulse' in dffKey:
        roiScores = anatData['roiScores']['NMF']
        # if 'dendConfidence' in anatData:
        #     roiScores = anatData['roiScores']['NMF']
        # else:
        #     try:
        #         roiScores = anatData['roiScores']['manual']
        #     except:
        #         roiScores = anatData['roiScores']['NMF']
    else:
        try:
            roiScores = anatData['roiScores']['manual']
        except:
            print('<<< dffKey = ' +str(dffKey)+' error in getGoodRois()>>>')

    nRois = len(roiScores.keys())

    anatConf = []
    for roi in range(nRois):
        conf = roiScores[roi]
        if useScores:
            if conf in goodROI_scores:
                anatConf.append(True)
            else:
                anatConf.append(False)

        else:
            if not switch:
                if conf >= keepTrsh:
                    anatConf.append(True)
                else:
                    anatConf.append(False)
            else:
                if conf <= keepTrsh:
                    anatConf.append(True)
                else:
                    anatConf.append(False)
                    
    return np.array(anatConf)
       
def readFrames(cap, start, end):
    imgs = [];
    cap.set(cv2.CAP_PROP_POS_FRAMES, start);
    for a in range(start, end):
        _, frame = cap.retrieve();
        imgs.append(cv2.cvtColor(frame[:,512:,:], cv2.COLOR_BGR2GRAY));
        cap.grab()
    return np.dstack(imgs)  

def maskEdge(mask):
    dims = mask.shape
    edges = np.zeros([dims[0],dims[1]])
    edges = edges + np.vstack([np.zeros([1,dims[1]]), np.diff(mask, axis = 0)])
    edges = edges + np.hstack([np.zeros([dims[0],1]), np.diff(mask, axis = 1)])
    
    #edgeWidth = 100
    #x = np.where(np.argmax(edges,axis=0).astype(bool) == True)[0][0]-edgeWidth
    #y = np.where(np.argmax(edges,axis=1).astype(bool) == True)[0][-1]+edgeWidth
    #bottomLeft = np.array([x,y])

    #x = np.where(np.argmax(edges,axis=0).astype(bool) == True)[0][-1]+edgeWidth
    #y = np.where(np.argmax(edges,axis=1).astype(bool) == True)[0][0]-edgeWidth
    #topRight = np.array([x,y])

    #points = [bottomLeft,topRight]
    #for i in range(2):
    #    for j in range(2):
    #        if points[i][j]>511:
    #            points[i][j] = 511
    #        elif points[i][j]<1:
    #            points[i][j] = 1

    edges = edges.astype(bool)
    return edges#, points

def tsep_traces(data, dtype = 'consensus_NMFtc_dff_bc_decon_sp', maskKey = 'consensus_NMFtc_mask', Z=False, thresh=-np.inf, 
                binary=False, nanEdgeTrials=False, gatherPreTrial = False, nanMaskFrs = True, nanFractionTraces = True, maxNaNFraction = 0.5, verbose = True, ddof = 0):

    nTrials = len(data['trial_frames'])
    fr = data['fr']
    dff = copy.deepcopy(data[dtype])
    mask = copy.deepcopy(data[maskKey])
    if nanEdgeTrials:
        nTrialSI= len(data['trial_frames'])
        startFrame = data['trial_frames'][nanEdgeTrials].astype(int)[-1]
        endFrame = data['trial_frames'][nTrialSI - nanEdgeTrials].astype(int)[-1]
        dff[:,:startFrame]=np.nan
        dff[:,endFrame:]=np.nan
    if Z:
        if Z=='standard_noise':
            for m in range(dff.shape[0]):
                dff[m,:] = zscore(dff[m,:], nan_policy='omit')

        elif Z=='baseline_noise':
            for m in range(dff.shape[0]):
                mask = np.isfinite(dff[m,:])
                dff2 = dff[m,mask]
                sn = caiman_baseline_noise(dff2,range_ff = [0.25, 0.5])
                # ##copied from caiman deconvolution module
                # range_ff = [0.25, 0.5]
                # ff, Pxx = scipy.signal.welch(dff2)
                # ind1 = ff > range_ff[0]
                # ind2 = ff < range_ff[1]
                # ind = np.logical_and(ind1, ind2)
                # Pxx_ind = Pxx[ind]
                # sn =  np.sqrt(np.exp(np.mean(np.log(old_div(Pxx_ind, 2)))))
                # ##copied from caiman deconvolution module
                dff[m, :]=dff[m,:]/sn
        elif Z=='std_only':
            with warnings.catch_warnings():
                # all-NaN rows give 0/0 -> nan; expected, not a bug
                warnings.filterwarnings('ignore', message='invalid value encountered in divide')
                for m in range(dff.shape[0]):
                    dff[m,:] = dff[m,:] / np.nanstd(dff[m,:], ddof = ddof)
        else:
            if verbose:
                print('Z-score method not found: skipping')
    if nanFractionTraces:
        for m in range(dff.shape[0]):
            if np.mean(np.isfinite(dff[m,:]))<maxNaNFraction:
                dff[m,:]=np.nan
    if nanMaskFrs:
        dff[mask]=np.nan

    nRois = dff.shape[0]
    tsep_dff = []
    for t in range(nTrials):
        badTrialFlag = False
        try:
            if not np.isnan(data['trial_frames'][t]).all():
                frames = data['trial_frames'][t].astype(int) - 1
                for m in range(dff.shape[0]):
                    stFrame = frames[0]
                    edFrame = frames[-1]
                    if edFrame == -1:
                        badTrialFlag = True
                        edFrame = int(fr*25)
                    #try:
                    if not badTrialFlag:
                        if gatherPreTrial:
                            if t != 0:
                                stFrame = int(frames[0]-fr*4)
                                trace = copy.deepcopy(dff[m, stFrame:edFrame])
                            else:
                                firstTrialNans = np.ones(int(fr*4))
                                trace = np.hstack([firstTrialNans,dff[m, stFrame:edFrame]])
                        else:
                            trace = copy.deepcopy(dff[m,stFrame:edFrame])
                    else:
                        firstTrialNans = np.ones(int(fr*4))
                        trace = np.ones(edFrame)*np.nan
                    mask = np.isfinite(trace)
                    temp = trace[mask]
                    if np.isfinite(thresh):
                        temp = temp - thresh
                        temp[temp<0]=0
                    if binary:
                        temp[temp>binary] = binary
                    trace[mask]=temp
                    #except:
                    #trace = np.ones([nRois,1000])*np.nan
                    if gatherPreTrial:
                        if t == 0:
                            if not badTrialFlag:
                                dff[m, stFrame:edFrame]=trace[len(firstTrialNans):]
                        else:
                            dff[m, stFrame:edFrame]=trace
                    else:
                        dff[m, stFrame:edFrame]=trace
                try:
                    tsep_dff.append(dff[:,stFrame:edFrame])
                except:
                    tsep_dff.append(np.ones([nRois,1000])*np.nan)
            else:
                tsep_dff.append(np.ones([nRois,1000])*np.nan)
        except IndexError:
            tsep_dff.append(np.ones([nRois,1000])*np.nan)


    try:
        knownIssues = list(data['knownIssues'].keys())
        if 'nanTrialInsertion' in knownIssues:
            if verbose:
                print("<<WARNING>> [tsep_traces()] "+str(data['anmID'])+" "+str(data['session'])+" INSERTING A NANED TRIAL")
                for k in knownIssues:
                    print("            KNOWN ISSUE "+str(k)+" "+str(data['knownIssues'][k]))
            nanedTrial = np.ones([nRois,1000])*np.nan
            trial = data['knownIssues']['nanTrialInsertion']
            desc = data['knownIssues']['description']
            if desc == 'file_corruption':
                tsep_dff.insert(trial, nanedTrial)
            elif desc == 'early_stoppage':
                tsep_dff.append(nanedTrial)
        if 'ignore_session' in knownIssues:
            if verbose:
                print("<<WARNING>> [tsep_traces()] "+str(data['anmID'])+" "+str(data['session'])+" IGNORING SESSION")
                for k in knownIssues:
                    print("            KNOWN ISSUE "+str(k)+" "+str(data['knownIssues'][k]))
            nTrials = data['knownIssues']['ignore_session']
            nanedSession = np.ones([nTrials,nRois,1000])*np.nan
            tsep_dff = nanedSession
    except:
        print("<<ERROR>>   [tsep_traces()] "+str(data['anmID'])+" "+str(data['session'])+" out of date files or some other issue")
    return tsep_dff

def get_tsep_wFlag(anatData, dffKey, maskKey, Z, thresh, binary, nanEdgeTrials, gatherPreTrial, nanMaskFrs, \
                   nanFractionTraces, maxNaNFraction, flag_session, flag_trials, flagParams, verbose = True):
    tsep = list()
    flags = {}
    for sess in list(anatData['sessions'].keys()):
        ##################################################################
        #Collect all ROI traces for a given session
        tsep_sess = tsep_traces(anatData['sessions'][sess], dtype = dffKey, maskKey = maskKey, \
                                Z=Z, thresh=thresh, binary=binary, \
                                nanEdgeTrials=nanEdgeTrials, gatherPreTrial = gatherPreTrial, \
                                nanMaskFrs = nanMaskFrs, nanFractionTraces = nanFractionTraces, \
                                maxNaNFraction = maxNaNFraction, verbose = verbose)
        ##################################################################
        #Just warnings if there are trial or ROI mismatches
        trial_frames = copy.deepcopy(anatData['sessions'][sess]['trial_frames'])     
        if verbose:
            if not len(trial_frames) == len(tsep_sess):
                print("<<WARNING>> [get_tsep_wFlag()] "+str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+\
                      str(anatData['sessions'][sess]['session'])+" len(trial_frames) = "+str(len(trial_frames))+\
                        " len(tsep_sess) = "+str(len(tsep_sess)))
        tc_masks = copy.deepcopy(anatData['sessions'][sess][flagParams['traceKey_mask']])
        if verbose:
            if not tsep_sess[0].shape[0] == tc_masks.shape[0]:
                print("<<WARNING>> [get_tsep_wFlag()] "+str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+\
                      str(anatData['sessions'][sess]['session'])+" tsep_sess[0].shape[0] = "+str(tsep_sess[0].shape[0])+\
                        " tc_masks.shape[0] = "+str(tc_masks.shape[0]))
        ##################################################################
        # This is where we nan frames that are flagged either with the flag_session or flag_trials
        if flag_session or flag_trials:
            flags[sess] = {}
            flags[sess]['session'] = copy.deepcopy(anatData['sessions'][sess]['session'])
            flags[sess]['byROI_sess_flags'] = np.zeros((tc_masks.shape[0]),dtype=bool)
            flags[sess]['byROI_trial_flags'] = np.zeros((tc_masks.shape[0],len(tsep_sess)),dtype=bool)
            for ROI in range(tc_masks.shape[0]):
                traceFlags = trace_flagging(tc_masks[ROI,:],float(anatData['sessions'][sess]['fr']),flagParams,trial_frames,\
                                            str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+\
                                                str(anatData['sessions'][sess]['session'])+" ROI"+str(ROI), verbose)
                flags[sess]['byROI_sess_flags'][ROI] = traceFlags['flagged_sess']
                for t in range(traceFlags['ntrials']):
                    flags[sess]['byROI_trial_flags'][ROI,t] = traceFlags['trials'][t]['flagged_trial']
                    if (traceFlags['flagged_sess'] and flag_session) or (traceFlags['trials'][t]['flagged_trial'] and flag_trials):
                        tsep_sess[t][ROI,:] = np.nan
        ##################################################################
        tsep.append(tsep_sess)
    return tsep, flags

def get_tsep(anatData, dffKey, maskKey = 'consensus_NMFtc_mask', onlyGoodRois = True, goodROI_scores = [1], \
            Z = False, thresh = -np.inf, binary = False, nanEdgeTrials = False, gatherPreTrial = False, \
            nanMaskFrs = True, nanFractionTraces = True, maxNaNFraction = 0.5, \
            full_tsep = False, flag_session = False, flag_trials = False, \
            flagParams = {}, returnFlags = False, verbose = True):

    # [get_tsep()] NEW USAGE NOTES:
    # 1) full_tsep toggle generates tsepFull and tsepFullMatch (in addition to tsep) where tsepFull wont clip the trials to the min length
    #    but things have to be adjusted a little in the ROI cleanup and it will output tsep, tsepFull, tsepFullMatch
    #    tsep is the normal version, tsepFull is organized by tsepFull[session][trial][ROI] 
    #    while the tsepFullMatch should match the tsep but is derived from tsepFull as a check
    # 2) flag_session will use the flagParams and just NaN the whole ROI/Session if it hits the various thresholds
    # 3) flag_trial will only NaN trials that cross the trial % frame missing threshold 
    # 4) nanMaskFrs uses any TRUE maskKey frames to NaN individual trace frames. The consensus_NMFtc_mask default maskKey was derived from NMF drop out frames and missing deep interp frames
    # 5) nanFractionTraces will NaN any full, cross-session ROI traces that have maxNaNFraction of NaNs across the entire ROI trace/all trials
    # 6) get_tsep_wFlag() is just compartmentalizing tsep_traces() along with the flagging function trace_flagging()
    
    if not flagParams:
        if verbose:
            print("<<WARNING>> [get_tsep()] USING DEFAULT FLAGPARAMS")
        flagParams = generate_default_flagParams(maskKey = maskKey)

    tsep, tsepFlags = get_tsep_wFlag(anatData, dffKey, maskKey, Z, thresh, binary, \
                                     nanEdgeTrials, gatherPreTrial, nanMaskFrs, nanFractionTraces, maxNaNFraction, \
                                        flag_session, flag_trials, flagParams, verbose = verbose)
    if full_tsep:
        sess = list(anatData['sessions'].keys())[0]
        tsepFull = copy.deepcopy(tsep)
        if onlyGoodRois:
            include = getGoodRois(anatData, dffKey = dffKey, useScores = True, goodROI_scores = goodROI_scores, verbose = verbose)
            if not include.shape[0] == tsep[0][0].shape[0]:
                if verbose:
                    print("<<WARNING>> [get_tsep()] "+str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+str(anatData['sessions'][sess]['session'])+" ROI COUNT MISMATCH")
                    print("            #ROIs in getGoodRois = "+str(include.shape[0]))
                    print("            #ROIs in tsep =        "+str(tsep.shape[1]))
                try2Fix=False
                if include.shape[0] > tsep[0][0].shape[0]:
                    try2Fix = True
                    for i in range(tsep[0][0].shape[0]-1,include.shape[0]):
                        if verbose:
                            print("            ROI "+str(i)+" SCORE: "+str(anatData['roiScores']['NMF'][i]))
                        if anatData['roiScores']['NMF'][i] in goodROI_scores:
                            try2Fix=False
                if try2Fix:
                    if verbose:
                        print("            ATTEMPTING TO FIX")
                    include = include[0:tsep[0][0].shape[0]]
                    if verbose:
                        print("            #ROIs in getGoodRois = "+str(include.shape[0]))
            for s in range(len(tsep)):
                for t in range(len(tsep[s])):
                    tsepFull[s][t] = np.ones((include.shape[0],tsep[s][t].shape[1]))*np.nan
                    tsepFull[s][t] = copy.deepcopy(tsep[s][t][include,:])
    tsep = get_tsepMatrix(tsep)
    if onlyGoodRois:
        include = getGoodRois(anatData, dffKey = dffKey, useScores = True, goodROI_scores = goodROI_scores, verbose = verbose)
        sess = list(anatData['sessions'].keys())[0]
        if not include.shape[0] == tsep.shape[1]:
            if verbose:
                print("<<WARNING>> [get_tsep()] "+str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+str(anatData['sessions'][sess]['session'])+" ROI COUNT MISMATCH")
                print("            #ROIs in getGoodRois = "+str(include.shape[0]))
                print("            #ROIs in tsep =        "+str(tsep.shape[1]))
            try2Fix=False
            if include.shape[0] > tsep.shape[1]:
                try2Fix = True
                for i in range(tsep.shape[1]-1,include.shape[0]):
                    if verbose:
                        print("            ROI "+str(i)+" SCORE: "+str(anatData['roiScores']['NMF'][i]))
                    if anatData['roiScores']['NMF'][i] in goodROI_scores:
                        try2Fix=False
            if try2Fix:
                if verbose:
                    print("            ATTEMPTING TO FIX")
                include = include[0:tsep.shape[1]]
                if verbose:
                    print("            #ROIs in getGoodRois = "+str(include.shape[0]))
        tsep = tsep[:,include,:]

    if full_tsep:
        tsepFullMatch=np.ones_like(tsep)*np.nan
        for r in range(tsep.shape[1]):
            trial=-1
            for s in range(len(tsepFull)):
                for t in range(len(tsep[s])):
                    if not tsepFull[s][t].shape[0] == tsep.shape[1]:
                        print("<<ERROR>>   [get_tsep()] "+str(anatData['sessions'][sess]['anmID'])+" sess"+str(sess)+" "+str(anatData['sessions'][sess]['session'])+" tsepFull ROI MISMATCH")
                        print(tsepFull[s][t].shape[0])
                        print(tsep.shape[1])
                        raise Exception("ROI Mismatch")
                    trial+=1
                    tsepFullMatch[trial,r,:] = copy.deepcopy(tsepFull[s][t][r,0:tsep.shape[2]])

    if returnFlags:
        if full_tsep:
            return tsep, tsepFull, tsepFullMatch, tsepFlags
        else:
            return tsep, tsepFlags
    else:
        if full_tsep:
            return tsep, tsepFull, tsepFullMatch
        else:
            return tsep

def get_tracker(behavior):
    mats = behavior['mats']
    recorded_days = behavior['recorded_days']
    ttracker = []
    for d in range(len(recorded_days)):
        mat = mats[d]
        ttypes = copy.deepcopy(mat['TrialTypes'])
        nTrials = mat['nTrials']
        earlyLicks = [np.intersect1d(states,[3,5]).any() for states in mat['RawData']['OriginalStateData']]
        wrong = [13 in states for states in mat['RawData']['OriginalStateData']]
        noAttempt = [12 in states for states in mat['RawData']['OriginalStateData']]
        wrong = np.logical_xor(wrong, noAttempt)
        ttypes[earlyLicks] = ttypes[earlyLicks] + 10
        ttypes[wrong] = ttypes[wrong] + 5
        ttypes[noAttempt] = ttypes[noAttempt] + 20
        
        #cheat sheet
        #0,1 - correct right,left
        #5,6 - incorrect right,left
        #10,11 - correct but early lick right,left
        #15,16 - incorrect and early lick right,left
        #20,21 - no attempt right, left
        #30,31 - no attampt and early lick right, left
        
        for t in range(nTrials):
            ttracker.append([recorded_days[d],t,ttypes[t]])
    ttracker = np.array(ttracker)
    return ttracker
    
def get_dendDays(anmData):
    dendData = anmData['dendrites']
    somaData = anmData['somas']
    behavior = anmData['behavior']
    tracker = get_tracker(behavior)
    try:
        dendRecDayIdxs = dendData['daysRelativeToShift']

        rec_days = behavior['recorded_days']
        dendTracker = np.array([tracker[t,0] in dendRecDayIdxs for t in range(tracker.shape[0])])
        dendTracker = tracker[dendTracker,:]

        dendDays = []
        for d, day in enumerate(rec_days):
            if day in dendRecDayIdxs:
                nDays = np.sum(tracker[:,0]==day)
                for d2 in range(nDays):
                    dendDays.append(True)
            else:
                nDays = np.sum(tracker[:,0]==day)
                for d2 in range(nDays):
                    dendDays.append(False)
        dendDays = np.array(dendDays)
    except:
        dendDays = np.zeros(tracker.shape[0]).astype(bool)
    
    return dendDays

def identify_rad_3d(anmData, fraction = 0.9):

    mouth = anmData['mouth']
    #########################################################
    if mouth[2]>0:
        mouth[2] = -mouth[2]
    #########################################################
    rp = anmData['port_locs']['pre_shift']['right']
    lp = anmData['port_locs']['pre_shift']['left']
    rpPost = anmData['port_locs']['post_shift']['right']
    lpPost = anmData['port_locs']['post_shift']['left']

    rpDist = np.sqrt((mouth[0]-rp[0])**2 +(mouth[1]-rp[1])**2+(mouth[2]-rp[2])**2)
    lpDist = np.sqrt((mouth[0]-lp[0])**2 +(mouth[1]-lp[1])**2+(mouth[2]-lp[2])**2)
    lpPostDist = np.sqrt((mouth[0]-lpPost[0])**2 +(mouth[1]-lpPost[1])**2+(mouth[2]-lpPost[2])**2)
    rpPostDist = np.sqrt((mouth[0]-rpPost[0])**2 +(mouth[1]-rpPost[1])**2+(mouth[2]-rpPost[2])**2)
    closestPort = np.min([rpDist, lpDist, lpPostDist, rpPostDist])

    rad = closestPort*fraction
    return rad

def get_dists(lick_type): #wants lists of a certain lick type
    dists = []
    for d in range(len(lick_type)):
        _dists_ = []
        wData = lick_type[d]
        for t in range(len(wData)):
            test = wData[t]
            origin = test[0,:3]
            dists_ = []
            for i in range(test.shape[0]):
                dist = np.sqrt((origin[0] - test[i,0])**2 + (origin[1] - test[i,1])**2 + (origin[2] - test[i,2])**2)
                dists_.append(dist)
            dists_ = butter_lowpass_filter(dists_, 200, 500, order = 5)
            _dists_.append(dists_)
        dists.append(_dists_)
    return dists

def find_first(array, distance):
    for f in range(array.shape[0]):
        if array[0]<= distance:
            if array[f] > distance:
                first = int(f)
                break
    else:
        first = np.nan
    return first

def time_sort(array):
    argmaxs = np.array([np.nanargmax(array[i,:]) for i in range(array.shape[0])])
    test = copy.deepcopy(array)
    straces = np.zeros(test.shape)
    place = 0
    for idx in range(np.max(argmaxs)):
        if idx in argmaxs:
            idxs = np.where(idx == argmaxs)[0]
            for i in idxs:
                straces[place,:] = test[i,:]
                place = place+1
    return straces

def get_radIdxs_3d(anmData, lickArray, rad, count_if_noCross = True):
    maxIdx = lickArray.shape[1]
    mouth = anmData['mouth']
    #########################################################
    if mouth[2]>0:
        mouth[2] = -mouth[2]
    #########################################################
    mdists = []
    for t in range(lickArray.shape[0]):
        tdata = lickArray[t,:,:]
        if np.isnan(tdata).all():
            mdists.append([np.nan]*maxIdx)
        else:
            dists_ = []
            for f in range(tdata.shape[0]):
                try:
                    dists_.append(np.sqrt((tdata[f,0] - mouth[0])**2 + (tdata[f,1] - mouth[1])**2 + (tdata[f,2] - mouth[2])**2))
                except:
                    dists_.append(np.nan)
            mdists.append(dists_)
    mdists = np.array(mdists)

    rads = []
    for t in range(mdists.shape[0]):
        r = find_first(mdists[t,:], rad)
        if np.isnan(r):
            if count_if_noCross:
                try:
                    r = np.nanargmax(mdists[t,:])
                    rads.append(int(r))
                except ValueError:
                    rads.append(np.nan)
            else:
                rads.append(np.nan)
        else:
            rads.append(int(r))    
    return np.array(rads)

def get_compassCoors(anmData, lickArray, radI):
    el = []
    az = []
    mouth = anmData['mouth']
    ###################################################
    if mouth[2]>0:
        mouth[2] = -mouth[2]
    ###################################################
    for t in range(lickArray.shape[0]):
        if np.isfinite(radI[t]):
            i = int(radI[t])
            #el.append(lickArray[t,i-1,2]-mouth[2])
            el.append(np.arctan2((mouth[2]-lickArray[t,i-1,2]),(mouth[1]-lickArray[t,i-1,1])))
            az.append(np.arctan2((mouth[1]-lickArray[t,i-1,1]),(mouth[0]-lickArray[t,i-1,0])))
        else:
            el.append(np.nan)
            az.append(np.nan)
    el = -np.rad2deg(np.array(el))
    az = -np.rad2deg(np.array(az))
    return az, el
    
def get_lickAngle(azi, ele, lda):
    """
    1. Finds the projection of (azi, ele) onto the LDA line in 2D.
    2. Converts both points to 3D unit vectors.
    3. Calculates the 3D angle (Cosine Similarity) between them.
    """
    # --- Step 1: Find the 2D Projection in Angle Space ---
    # Line: w0*az + w1*el + b = 0
    w = lda.coef_[0]
    b = lda.intercept_[0]
    
    # Perpendicular projection of (azi, ele) onto the LDA line
    # Using the standard point-to-line projection formula
    denom = w[0]**2 + w[1]**2
    proj_azi = (w[1]*(w[1]*azi - w[0]*ele) - w[0]*b) / denom
    proj_ele = (w[0]*(w[0]*ele - w[1]*azi) - w[1]*b) / denom

    # --- Step 2: Convert to 3D Unit Vectors ---
    def to_unit_vec(az, el):
        az_r, el_r = np.radians(az), np.radians(el)
        return np.array([
            np.cos(el_r) * np.cos(az_r),
            np.cos(el_r) * np.sin(az_r),
            np.sin(el_r)
        ]).T

    u_lick = to_unit_vec(azi, ele)
    u_proj = to_unit_vec(proj_azi, proj_ele)

    # --- Step 3: Cosine Similarity / Angular Distance ---
    # dot(u1, u2) = cos(theta)
    # We use the dot product to find the actual 3D angle (theta)
    if u_lick.ndim == 1:
        dot_prod = np.dot(u_lick, u_proj)
    else:
        # Vectorized dot product across licks
        dot_prod = np.sum(u_lick * u_proj, axis=1)
        
    # Angular distance in degrees
    angle_deg = np.degrees(np.arccos(np.clip(dot_prod, -1, 1)))

    # --- Step 4: Add the Sign (L vs R) ---
    # Determine which side of the LDA boundary the lick is on
    # (Numerator of the distance formula gives the sign)
    sign = np.sign(w[0]*azi + w[1]*ele + b)
    
    return angle_deg * sign

def get_aziEle(anmData, trajs, radFraction = 1, count_of_noCross = False):
    rad = identify_rad_3d(anmData, fraction = radFraction)
    radI = get_radIdxs_3d(anmData, trajs, rad, count_if_noCross = count_of_noCross)
    az, el = get_compassCoors(anmData, trajs, radI)
    return az, el
    
def get_contactCoors(anmData):
    licks = get_contacts(anmData['trajs'])#, 1)
    tracker = get_tracker(anmData)

    rad = identify_rad_3d(anmData, fraction = 1)
    radI = get_radIdxs_3d(anmData, licks, rad, count_if_noCross = False)
    cTimes = np.hstack([get_cTimes(mat, sub_goCue = True) for mat in anmData['mats']])

    coors = []
    for t in range(tracker.shape[0]):
        cTime = cTimes[t]
        if np.sum(np.isfinite(licks[t,:,3]))>5: # if there is tracking there
            i = np.nanargmin(np.abs(licks[t,:,3]-cTime))
            coors.append(np.hstack([licks[t,i-1,0], licks[t,i-1,1], licks[t,i-1,2]]))
        else:
            coors.append(np.ones(3)*np.nan)
    
    return np.vstack(coors)
    
def get_contactPortDists(anmData, returnAssignedPort = False):
    coors = get_contactCoors(anmData)

    portLocs = anmData['port_locs']
    tracker = get_tracker(anmData)
    sLabels = ['right', 'left']
    dists = []
    for t in range(tracker.shape[0]):
        day = tracker[t,0]
        ttype = tracker[t,2]%5

        if day>0 or np.logical_and(day==0, tracker[t,1]>99):
            shiftPhase = 'post_shift'
        else:
            shiftPhase = 'pre_shift'

        if returnAssignedPort:
            port = portLocs[shiftPhase][sLabels[ttype]]
        else:
            if tracker[t,2]==5 or tracker[t,2]==15:
                port = portLocs[shiftPhase]['left']
            elif tracker[t,2]==6 or tracker[t,2]==16:
                port = portLocs[shiftPhase]['right']
            else:
                port = portLocs[shiftPhase][sLabels[ttype]]
                
        if np.sum(np.isfinite(coors[t,:])):
            dists.append(np.sqrt(np.sum((coors[t,:]-port)**2)))
        else:
            dists.append(np.nan)
        
    return np.hstack(dists)
    
def get_minDistToPort(anmData, lickType, portPhase, portSide):
    if lickType == 'contacts' or lickType == 'contact' or lickType == 'Contacts' or lickType == 'Contact':
        licks = get_contacts(anmData['trajs'])
    elif lickType == 'firsts' or lickType == 'first' or lickType == 'First' or lickType == 'Firsts':
        licks == get_firsts(anmData['trajs'])
    else:
        print('<<< error in get_minDistToPort() >>>')
        print('lickType must == firsts or contacts')
        
    try:
        port = anmData['port_locs'][portPhase][portSide]
    except:
        print('<<< error in get_minDistToPort() >>>')
        print('portPhase must == pre_shift or post_shift')
        print('portSide must == right or left')
        
    closestDist = []
    for t in range(licks.shape[0]):
        ttraj = licks[t,:,:3]
        nFinite = np.sum(np.isfinite(ttraj[:,0]))
        if nFinite>2:
            tDist = np.hstack([np.sqrt(np.sum((ttraj[f,:]-port)**2)) for f in range(nFinite)])
            closestDist.append(np.nanmin(tDist))
        else:
            closestDist.append(np.nan)
    
    return np.hstack(closestDist)
    
def get_lickDir(anmData, tracker, useErrorData = False):
    lickDir = []
    if useErrorData:
        licks = get_contacts(anmData['trajs'])
        azi, ele = get_aziEle(anmData, licks)
        me, de = get_errors_LDA(azi, ele, tracker)
    for t in range(tracker.shape[0]):
        ttype = tracker[t,2]
        if ttype == 0 or ttype == 10:
            lickDir.append(0)
        elif ttype == 1 or ttype == 11:
            lickDir.append(1)
        elif ttype == 5 or ttype == 15:
            if useErrorData:
                if me[t]:
                    lickDir.append(0)
                else:
                    lickDir.append(1)
            else:
                lickDir.append(1)
        elif ttype == 6 or ttype == 16:
            if useErrorData:
                if me[t]:
                    lickDir.append(1)
                else:
                    lickDir.append(0)
            else:
                lickDir.append(0)
        else:
            lickDir.append(np.nan)
    return np.hstack(lickDir)

def get_relDistances(azi, ele, tracker, absDist = False, returnLDA = False, returnProjs = False):
    
    def proj2d(points, origSlope, origB):
        xrange = np.linspace(np.nanpercentile(points[:,0], 5), np.nanpercentile(points[:,0], 95),100)
        projs = []
        for p in range(points.shape[0]):
            pt = points[p,:]
            if np.sum(np.isnan(pt))>0:
                projs.append([np.nan, np.nan])
                continue
            recM = -(1/origSlope)
            b2 = pt[1]-(pt[0]*recM)
            y2 = recM*xrange+b2

            xProj = (origB-b2)/(recM-origSlope)
            yProj = xProj*recM+b2
            projs.append(np.hstack([xProj, yProj]))
        return np.vstack(projs)

    coors = np.vstack([azi,ele]).T
    nanMask = np.any(np.isfinite(coors),axis = 1)
    fCoors = coors[nanMask,:]
    
    preMask = get_behaviorShiftMask(tracker, 'pre')
    preT = tracker[np.logical_and(nanMask, preMask),:]
    preCoors = coors[np.logical_and(nanMask, preMask),:]

    crclMask = np.logical_or(preT[:,2]==0, preT[:,2]==1)

    preppedCoors = preCoors[crclMask,:]
    Y = preT[crclMask,2]

    lda = LDA(solver = 'svd').fit(preppedCoors, Y)
    
    intercept = lda.intercept_
    coefs = lda.coef_[0]

    origSlope = -coefs[0]/coefs[1]
    origB = -intercept/coefs[1]

    projs = proj2d(coors, origSlope, origB)
    
    if absDist:
        dists = np.hstack([np.sqrt(np.sum((projs[t,:]-coors[t,:])**2)) for t in range(coors.shape[0])])
    else:
        dists = []
        for t in range(coors.shape[0]):
            direction = projs[t,:][0]-coors[t,:][0]
            if direction<0:
                dists.append(-np.sqrt(np.sum((projs[t,:]-coors[t,:])**2)))
            else:
                dists.append(np.sqrt(np.sum((projs[t,:]-coors[t,:])**2)))
        dists = np.hstack(dists)

    
    if returnLDA:
        if returnProjs:
            return dists, lda, projs
        else:
            return dists, lda
    elif returnProjs and not returnLDA:
        return dists, projs
    else:
        return dists

def get_firsts(trajs):
    maxIdx = 200
    firsts = []
    for d in range(len(trajs)):
        traj = trajs[d][0]
        for t in range(len(traj)):
            tdata = traj[t][0]
            data = np.ones((maxIdx, 4))*np.nan
            tMaxIdx = tdata.shape[0]
            if tMaxIdx < maxIdx:
                #data[:tMaxIdx, :] = tdata[:,:4]
                data[:tMaxIdx, :] = tdata[:,:]
            elif tMaxIdx >= maxIdx:
                #data = tdata[:maxIdx,:4]
                data = tdata[:maxIdx,:]
            firsts.append(data)
    firsts = np.array(firsts)
    return firsts

def get_contacts(trajs):
    maxIdx = 200
    contacts = []
    for d in range(len(trajs)):
        traj = trajs[d][0]
        for t in range(len(traj)):
            tdata = traj[t][-1]
            data = np.ones((maxIdx, 4))*np.nan
            tMaxIdx = tdata.shape[0]
            if tMaxIdx < maxIdx:
                data[:tMaxIdx, :] = tdata[:,:4]
            elif tMaxIdx >= maxIdx:
                data = tdata[:maxIdx,:4]
            contacts.append(data)     
    contacts = np.array(contacts)
    return contacts
       
def get_consumption_licks(trajs, lickIdx):
    maxIdx = 200
    licks = []
    for d in range(len(trajs)):
        traj = trajs[d][1]
        for t in range(len(traj)):
            data = np.ones((maxIdx, 4))*np.nan
            try:
                tdata = traj[t][lickIdx]
                tMaxIdx = tdata.shape[0]
                if tMaxIdx < maxIdx:
                    data[:tMaxIdx, :] = tdata[:,:4]
                elif tMaxIdx >= maxIdx:
                    data = tdata[:maxIdx,:4]
            except:
                pass
            licks.append(data)     
    licks = np.array(licks)
    return licks

def plot_lick_boundary_3d(azi, ele, lda, tracker):
    """
    Plots licks and the LDA decision boundary on a 3D unit sphere.
    """
    nanMask = np.all(np.isfinite(np.vstack([azi,ele])), axis = 0) 

    # 1. Convert Licks to 3D
    az_rad, el_rad = np.radians(azi), np.radians(ele)
    lx = np.cos(el_rad) * np.cos(az_rad)
    lx = lx[nanMask]
    ly = np.cos(el_rad) * np.sin(az_rad)
    ly = ly[nanMask]
    lz = np.sin(el_rad)
    lz = lz[nanMask]
    
    # Predict classes for coloring
    classes = lda.predict(np.column_stack([azi[nanMask], ele[nanMask]]))

    # 2. Generate the Boundary Line (Great Circle)
    # We find the normal to the LDA plane
    w, b = lda.coef_[0], lda.intercept_[0]
    m_az, m_el = np.nanmean(azi), np.nanmean(ele)
    
    def to_vec(a, e):
        ra, re = np.radians(a), np.radians(e)
        return np.array([np.cos(re)*np.cos(ra), np.cos(re)*np.sin(ra), np.sin(re)])

    # Get two points on the LDA line to define the plane
    az1, az2 = m_az - 20, m_az + 20
    el1, el2 = -(w[0]*az1 + b)/w[1], -(w[0]*az2 + b)/w[1]
    v1, v2 = to_vec(az1, el1), to_vec(az2, el2)
    
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)

    # Generate points along the great circle
    # We create a basis for the plane perpendicular to the normal
    if abs(normal[2]) < 0.9:
        v_basis = np.array([0, 0, 1])
    else:
        v_basis = np.array([0, 1, 0])
    
    b1 = np.cross(normal, v_basis)
    b1 /= np.linalg.norm(b1)
    b2 = np.cross(normal, b1)
    
    angles = np.linspace(0, 2*np.pi, 200)
    # This generates the full circle; we'll filter it to the lick region
    bx = b1[0] * np.cos(angles) + b2[0] * np.sin(angles)
    by = b1[1] * np.cos(angles) + b2[1] * np.sin(angles)
    bz = b1[2] * np.cos(angles) + b2[2] * np.sin(angles)

    # 3. Create the Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Draw a faint wireframe sphere for context
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    sx = np.cos(u)*np.sin(v)
    sy = np.sin(u)*np.sin(v)
    sz = np.cos(v)
    ax.plot_wireframe(sx, sy, sz, color="gray", alpha=0.1, linewidth=0.5)

    # Plot the licks
    labels = ['right', 'left']
    colors = ['b','r']
    for t in range(2):
        ax.scatter(lx[tracker[nanMask,2]==t], ly[tracker[nanMask,2]==t], lz[tracker[nanMask,2]==t], c=colors[t], label=labels[t], s=10, alpha=0.6)

    # Plot the boundary line (only the part near the licks)
    ax.plot(bx, by, bz, color='black', linewidth=2, label='LDA Boundary')

    # Formatting
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Spherical Lick Distribution & LDA Boundary')
    ax.scatter(0,0,0,c='k',s = 100)
    
    # Set view to focus on the data
    ax.view_init(elev=np.nanmean(ele), azim=np.nanmean(azi))
    ax.legend()
    
    plt.show()

def symmetric_orthogonalization(A):
    """
    Generated by Google Gemini
    Performs symmetric orthogoanalization on a set of vectors.

    Args:
        A: numpy.ndarray, matrix where each column is a vector to orthogonalize

    Returns:
        numpy.ndarray: symmetric orthogonalization matrix
    """
    S = A.T @ A  # Overlap matrix
    lam_s, l_s = np.linalg.eigh(S)  # Eigendecomposition
    lam_s_mat = np.diag(lam_s)  # Diagonal matrix of eigenvalues
    lam_sqrt_inv = np.sqrt(np.linalg.inv(lam_s_mat)) # Inverse square root of eigenvalues
    symm_orthog = l_s @ lam_sqrt_inv @ l_s.T # Symmetric orthogonalization matrix
    return A @ symm_orthog

def downSampleTsep(psTsep, factor = 2):
    fullTime = psTsep.shape[2]
    downSTime = int(fullTime/factor)
    
    new = []
    with warnings.catch_warnings():
        # start:stop bins can be entirely NaN near trial edges; that's expected, not a bug
        warnings.filterwarnings('ignore', message='Mean of empty slice')
        for i in range(0,psTsep.shape[2],factor):
            start = i
            stop = start+factor
            mean = np.nanmean(psTsep[:,:,start:stop], axis = 2)
            new.append(mean)
    
    return np.dstack(new)

def psTsep_func_old(tsep, times, fr, windowSize = [-1, 1], downSample = False, downSampleFactor = 2):
    aligned = []
    if type(windowSize) == list:
        windowSize = np.hstack(windowSize)
    winS = np.hstack(windowSize*fr).astype(int)
    nFr = np.sum(np.abs(winS))
    for t in range(tsep.shape[0]):
        if np.isfinite(times[t]):
            time = int(times[t]*fr)
            if time+winS[1]<tsep.shape[2]:
                if not downSample:
                    aligned.append(tsep[t,:,time+winS[0]:time+winS[1]])
                else:
                    pre = []
                    for tp in range(time-np.abs(winS[0]), time, downSampleFactor):
                        pre.append(np.nanmean(tsep[t,:,tp:tp+downSampleFactor], axis = 1))
                    post = []
                    for tp in range(time, time+winS[1], downSampleFactor):
                        post.append(np.nanmean(tsep[t,:,tp:tp+downSampleFactor], axis = 1))
                    
                    aligned.append(np.vstack([pre, post]).T)
            else:
                # if time+winS[0]<tsep.shape[2]:
                #     conglom = tsep[t,:,time+winS[0]:]
                #     frMissing = nFr-conglom.shape[1]
                #     if not downSample:
                #         conglom2 = np.hstack([conglom, np.ones([conglom.shape[0], frMissing])*np.nan])
                #     else:
                #         conglom2 = np.hstack([conglom, np.ones([conglom.shape[0], frMissing])*np.nan])
                #     aligned.append(conglom2)
                # else:
                if not downSample:
                    aligned.append(np.ones([tsep.shape[1], nFr])*np.nan)  
                else:
                    aligned.append(np.ones([tsep.shape[1], int(nFr/downSampleFactor)])*np.nan)  
        else:
            if not downSample:
                aligned.append(np.ones([tsep.shape[1], nFr])*np.nan)  
            else:
                aligned.append(np.ones([tsep.shape[1], int(nFr/downSampleFactor)])*np.nan)  
    aligned = np.array(aligned)
    return aligned

def psTsep_func(tsep, times, fr, windowSize=[-1, 1], 
                downSample=False, downSampleFactor=2):

    aligned = []
    if type(windowSize) == list:
        windowSize = np.hstack(windowSize)
    winS = (windowSize * fr).astype(int)
    
    preFrames = abs(winS[0])
    postFrames = winS[1]
    totalFrames = preFrames + postFrames
    
    preBins = len(list(range(-preFrames, 0, downSampleFactor)))
    postBins = len(list(range(0, postFrames, downSampleFactor)))
    expectedBins = preBins + postBins
    
    for t in range(tsep.shape[0]):
        if np.isfinite(times[t]):
            time = int(times[t] * fr)
    
            if time>tsep.shape[2]:
                if not downSample:
                    aligned.append(np.ones([tsep.shape[1], totalFrames])*np.nan)
                else:
                    aligned.append(np.ones([tsep.shape[1], expectedBins])*np.nan)
            else:
                if not downSample:
                    pre = tsep[t,:,time-preFrames:time]
                    post = []
                    for tp in range(time, time+postFrames):
                        if tp < tsep.shape[2]:
                            post.append(tsep[t,:,tp])
                        else:
                            post.append(np.ones(tsep.shape[1]) * np.nan)
                    aligned.append(np.vstack([list(pre.T),post]).T)
                else:
                    with warnings.catch_warnings():
                        # tp:tp+downSampleFactor bins can be entirely NaN near trial edges; expected, not a bug
                        warnings.filterwarnings('ignore', message='Mean of empty slice')
                        pre = []
                        for tp in range(time-preFrames, time, downSampleFactor):
                            pre.append(np.nanmean(tsep[t,:,tp:tp + downSampleFactor], axis=1))
                        post = []
                        for tp in range(time, time+postFrames, downSampleFactor):
                            if tp+downSampleFactor<tsep.shape[2]:
                                post.append(np.nanmean(tsep[t,:,tp:tp + downSampleFactor], axis=1))
                            else:
                                post.append(np.ones(tsep.shape[1])*np.nan)
                    aligned.append(np.vstack([pre, post]).T)
        else:
            if not downSample:
                aligned.append(np.ones((tsep.shape[1], totalFrames))*np.nan)
            else:
                aligned.append(np.ones((tsep.shape[1], expectedBins))*np.nan)
    return np.array(aligned)
 
def get_behaviorShiftMask(tracker, learningPhase, verbose = False):
    relDays = np.unique(tracker[:,0])
    ###
    if np.sum(relDays>=0)>0:
        shiftIdx = np.where(relDays>=0)[0][0]
    else:
        shiftIdx = len(relDays)-1
    ###
    if learningPhase == 'all':
        shiftMask = tracker[:, 0] < np.Inf
    elif learningPhase == 'pre':
        preDayMask = tracker[:, 0] < 0
        shiftDayMask = np.logical_and(tracker[:,0]==0, tracker[:,1]<=99)
        shiftMask = np.logical_or(preDayMask, shiftDayMask)

    ###
    elif learningPhase == 'day-1_beg':
        shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx-1], tracker[:,1]<=99)
    elif learningPhase == 'day-1_mid':
        t1 = tracker[:,1]>99
        t2 = tracker[:,1]<=199
        m1 = np.logical_and(t1,t2)
        shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx-1], m1)
    elif learningPhase == 'day-1_end':
        shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx-1], tracker[:,1]>199)
        
    elif learningPhase == 'preShift_early':
        if len(relDays[relDays<=relDays[shiftIdx]])>2:
            shiftMask = tracker[:, 0] <= relDays[shiftIdx-2]
        else:
            shiftMask = tracker[:,0] == relDays[0]

    elif learningPhase == 'preShift_late':
        m1 = tracker[:,0]<relDays[shiftIdx]
        m2 = tracker[:,0]>=relDays[shiftIdx-2]
        if len(relDays)==2 and relDays[0]<0 and relDays[1]<0:
            m2 = tracker[:,0]>relDays[0]
        shiftMask = np.logical_and(m1, m2)
        shiftDayMask = np.logical_and(tracker[:,0]==0, tracker[:,1]<=99)
        shiftMask = np.logical_or(shiftMask, shiftDayMask)
        #print(np.unique(tracker[shiftMask,0]))
    elif learningPhase == 'shiftDay_pre':
        shiftMask = np.logical_and(tracker[:,0]==0, tracker[:,1]<=99)
    elif learningPhase == 'shiftDay_post':
        shiftMask = np.logical_and(tracker[:,0]==0, tracker[:,1]>=100)
    elif learningPhase == 'shiftDay_post1':
        t1 = tracker[:,1]>99
        t2 = tracker[:,1]<=133
        m1 = np.logical_and(t1,t2)
        shiftMask = np.logical_and(tracker[:,0]==0, m1)
    elif learningPhase == 'shiftDay_post2':
        t1 = tracker[:,1]>133
        t2 = tracker[:,1]<=166
        m1 = np.logical_and(t1,t2)
        shiftMask = np.logical_and(tracker[:,0]==0, m1)
    elif learningPhase == 'shiftDay_post3':
        t1 = tracker[:,1]>166
        t2 = tracker[:,1]<=200
        m1 = np.logical_and(t1,t2)
        shiftMask = np.logical_and(tracker[:,0]==0, m1)
    ###

    elif learningPhase == 'early':
        nDaysPost = len(relDays>=0)
        shiftDayMask = np.logical_and(tracker[:, 0] == 0, tracker[:, 1] > 99) 
        if relDays[shiftIdx] != relDays[-1]:
            print(f'phase: {learningPhase}, days used: {[0, int(relDays[shiftIdx+1])]}')
            shiftMask = np.logical_or(shiftDayMask, tracker[:, 0] == relDays[shiftIdx+1]) 
        else:
            shiftMask = np.logical_or(tracker[:,0]>=0, shiftDayMask)

    elif learningPhase == 'late':
        try:
            shiftMask = tracker[:, 0] >= relDays[shiftIdx+2]
        except:
            shiftMask = np.zeros(tracker.shape[0]).astype(bool)
            print('empty fail no late post')

    elif learningPhase == 'day+1_beg':
        if relDays[shiftIdx] == 0:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx+1], tracker[:,1]<=99)
        else:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx], tracker[:,1]<=99)

    elif learningPhase == 'day+1_mid':
        t1 = tracker[:,1]>99
        t2 = tracker[:,1]<=199
        m1 = np.logical_and(t1,t2)
        if relDays[shiftIdx] == 0:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx+1], m1)
        else:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx], m1)
    elif learningPhase == 'day+1_end':
        if relDays[shiftIdx] == 0:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx+1], tracker[:,1]>199)
        else:
            shiftMask = np.logical_and(tracker[:,0]==relDays[shiftIdx], tracker[:,1]>199)            

    elif learningPhase == 'post':
        try:
            shiftDayMask = np.logical_and(tracker[:,0]==0, tracker[:,1]>99)
            shiftMask = np.logical_or(shiftDayMask, tracker[:,0]>0)
        except:
            shiftMask = np.zeros(tracker.shape[0]).astype(bool)
            print('empty fail no postShift')

    else:
        shiftMask = np.zeros(tracker.shape[0]).astype(bool)
        print(f'<<<Error: {learningPhase} learning phase not found>>>')
     
    if verbose:
        print(np.unique(tracker[shiftMask,0]))
    return shiftMask
       
def get_errors_LDA(azi, ele, tracker, returnLDA = False):
    
    coors = np.vstack([azi,ele]).T
    nanMask = np.any(np.isfinite(coors),axis = 1)
    fCoors = coors[nanMask,:]
    
    preMask = get_behaviorShiftMask(tracker, 'pre')
    preT = tracker[np.logical_and(nanMask, preMask),:]
    preCoors = coors[np.logical_and(nanMask, preMask),:]

    crclMask = np.logical_or(preT[:,2]==0, preT[:,2]==1)

    preppedCoors = preCoors[crclMask,:]
    Y = preT[crclMask,2]

    lda = LDA(solver = 'svd').fit(preppedCoors, Y)
    pred = lda.predict(fCoors)
    
    me = []
    de = []
    predCounter = 0
    for t in range(tracker.shape[0]):
        if np.isnan(azi[t]) or np.isnan(ele[t]):
            me.append(False)
            de.append(False)
        else:
            ttype = tracker[t,2]
            if ttype == 5 or ttype == 15:
                if pred[predCounter] == 0:
                    me.append(True)
                    de.append(False)
                else:
                    me.append(False)
                    de.append(True)
            elif ttype == 6 or ttype == 16:
                if pred[predCounter] == 1:
                    me.append(True)
                    de.append(False)
                else:
                    me.append(False)
                    de.append(True)
            else:
                me.append(False)
                de.append(False)
            predCounter+=1
    me = np.hstack(me)
    de = np.hstack(de)
    
    if returnLDA:
        return me, de, lda
    else:
        return me, de

def get_secondLickSide(beh, returnTime = False):
    bpodLicks = get_lickTimes(beh)
    cTimes = get_cTimes2(beh)
    licks = get_contacts(beh['trajs'])
    azi,ele = get_aziEle(beh, licks)
    tracker = get_tracker(beh)
    goCues = get_goCues(beh)
    sLickTrajs = get_consumption_licks(beh['trajs'],0)
    sLickTimes = sLickTrajs[:,:,3]
    sLick = np.ones(tracker.shape[0])*np.nan
    slTimes = np.ones(tracker.shape[0])*np.nan
    
    for t in range(bpodLicks.shape[2]):
        tLicks = bpodLicks[:,:,t]
        mask = np.isfinite(tLicks).all(axis=1)
        tLicks = tLicks[mask]
    
        slTime = sLickTimes[t,:]
        slTime = slTime[np.isfinite(slTime)]
        if len(slTime)>0:
            lMask = np.logical_and((tLicks[:,0]-goCues[t])>slTime[0], (tLicks[:,0]-goCues[t])<slTime[-1])
            if np.sum(lMask)>0:
                fIdx = np.where(lMask)[0][0]
                sTime, sSide = tLicks[fIdx,:]
                sLick[t] = sSide
                slTimes[t] = sTime-goCues[t]
    if returnTime:
        return sLick, slTimes
    else:
        return sLick

def get_MER_bySecondLick(beh):
    licks = get_contacts(beh['trajs'])
    azi,ele = get_aziEle(beh, licks)
    tracker = get_tracker(beh)
    me, de = get_errors_LDA(azi, ele, tracker)
    mer = np.logical_and(me, tracker[:,2]==5)
    slSides = get_secondLickSide(beh)
    
    afterRight = np.zeros(tracker.shape[0]).astype(bool)
    afterLeft = np.zeros(tracker.shape[0]).astype(bool)
    afterNone = np.zeros(tracker.shape[0]).astype(bool)
    for t in range(tracker.shape[0]):
        if mer[t]:
            sSide = slSides[t]
            if sSide == 0:
                afterRight[t]=True
            elif sSide == 1:
                afterLeft[t]=True
            else:
                afterNone[t]=True
    return afterRight, afterLeft, afterNone

def get_sl_errors(azi, ele, merOrig, lda, slSides):

    coors = np.vstack([azi,ele]).T
    nanMask = np.any(np.isfinite(coors),axis = 1)
    fCoors = coors[nanMask,:]
    
    pred = lda.predict(fCoors)
    
    me = []
    de = []
    predCounter = 0
    for t in range(merOrig.shape[0]):
        if np.isnan(azi[t]) or np.isnan(ele[t]):
            me.append(False)
            de.append(False)
        else:
            if merOrig[t]:
                if slSides[t] == 1:
                    if pred[predCounter] == 0:
                        me.append(True)
                        de.append(False)
                    else:
                        me.append(False)
                        de.append(True)
                else:
                    me.append(False)
                    de.append(False)
            else:
                me.append(False)
                de.append(False)
            predCounter+=1
    me = np.hstack(me)
    de = np.hstack(de)
    
    return me, de

def get_secondLickSide2(beh, returnTime = False, returnN = False, returnNextIdent = False):
    bpodLicks = get_lickTimes(beh)
    cTimes = get_cTimes2(beh)
    licks = get_contacts(beh['trajs'])
    azi,ele = get_aziEle(beh, licks)
    tracker = get_tracker(beh)
    goCues = get_goCues(beh)
    sLickTrajs = get_consumption_licks(beh['trajs'],0)
    sLickTimes = sLickTrajs[:,:,3]
    sLick = np.ones(tracker.shape[0])*np.nan
    slTimes = np.ones(tracker.shape[0])*np.nan
    slNHits = np.ones(tracker.shape[0])*np.nan
    nextIdent = np.ones(tracker.shape[0])*np.nan
    for t in range(bpodLicks.shape[2]):
        tLicks = bpodLicks[:,:,t]
        mask = np.isfinite(tLicks).all(axis=1)
        tLicks = tLicks[mask]
    
        slTime = sLickTimes[t,:]
        slTime = slTime[np.isfinite(slTime)]
        if len(slTime)>0:
            lMask = np.logical_and((tLicks[:,0]-goCues[t])>slTime[0], (tLicks[:,0]-goCues[t])<slTime[-1])
            if np.sum(lMask)>0:
                slNHits[t]=np.sum(lMask)
                fIdx = np.where(lMask)[0][0]
                sTime, sSide = tLicks[fIdx,:]
                sLick[t] = sSide
                slTimes[t] = sTime-goCues[t]

                if np.sum(lMask)>1:
                    nIdx = np.where(lMask)[0][1]
                    sTime, sSide = tLicks[nIdx,:]
                    nextIdent[t] = sSide
            else:
                if len(slTime)>25:
                    slNHits[t] = 0
    if returnTime:
        if returnN:
            if returnNextIdent:
                return sLick, slTimes, slNHits, nextIdent
            else:
                return sLick, slTimes, slNHits
        else:
            return sLick, slTimes
    
    elif returnN:
        if returnNextIdent:
            return sLick, slNHits, nextIdent
        else:
            return sLick, slNHits
    elif returnNextIdent:
        return sLick, nextIdent
    else:
        return sLick

def get_relDistFromLDA(azi, ele, lda):
    
    def proj2d(points, origSlipe, origB):
        xrange = np.linspace(np.nanpercentile(points[:,0], 5), np.nanpercentile(points[:,0], 95),100)
        projs = []
        for p in range(points.shape[0]):
            pt = points[p,:]
            if np.sum(np.isnan(pt))>0:
                projs.append([np.nan, np.nan])
                continue
            recM = -(1/origSlope)
            b2 = pt[1]-(pt[0]*recM)
            y2 = recM*xrange+b2

            xProj = (origB-b2)/(recM-origSlope)
            yProj = xProj*recM+b2
            projs.append(np.hstack([xProj, yProj]))
        return np.vstack(projs)
    
    coors = np.vstack([azi,ele]).T
    # nanMask = np.any(np.isfinite(coors),axis = 1)
    # preppedCorrs = coors[nanMask,:]

    coefs = lda.coef_[0]
    intercept = lda.intercept_
    origSlope = -coefs[0]/coefs[1]
    origB = -intercept/coefs[1]

    projs = proj2d(coors, origSlope, origB)
    #pred = lda.predict(coors)

    dists = []
    for t in range(coors.shape[0]):
        direction = projs[t,:][0]-coors[t,:][0]
        if direction<0:
            dists.append(-np.sqrt(np.sum((projs[t,:]-coors[t,:])**2)))
        else:
            dists.append(np.sqrt(np.sum((projs[t,:]-coors[t,:])**2)))
    dists = np.hstack(dists)

    return dists

def clean_subplots(nrows,ncols,figsize,facecolor='none',constrained_layout = False, sharey = False, sharex = False):
    if figsize:
        fig,ax = plt.subplots(nrows,ncols,figsize=figsize, sharey = sharey, sharex = sharex)
    else:
        fig,ax = plt.subplots(nrows,ncols, sharey = sharey, sharex = sharex)
    if hasattr(ax, 'shape'):
        if len(ax.shape)>1:
            for row in range(ax.shape[0]):
                for col in range(ax.shape[1]):
                    ax[row,col].set_facecolor(facecolor)
                    if facecolor=='none':
                        ax[row,col].patch.set_alpha(0.0) # Axes background
        else:
            for row in range(ax.shape[0]):
                ax[row].set_facecolor(facecolor)
                if facecolor=='none':
                    ax[row].patch.set_alpha(0.0) # Axes background
    else:
        ax.set_facecolor(facecolor)
        if facecolor=='none':
            ax.patch.set_alpha(0.0) # Axes background
    if constrained_layout:
        fig.set_constrained_layout(True)
    fig.set_facecolor(facecolor)
    if facecolor=='none':
        fig.patch.set_alpha(0.0) # Figure background
    return fig,ax

def get_go_contact_irf(data,ctf):
    from scipy.linalg import toeplitz
    from sklearn.linear_model import LinearRegression, RidgeCV, LassoCV, ElasticNetCV
    y=copy.deepcopy(data)
    tlength = y.shape[1]
    y[np.isnan(y)]=0
    iG=np.zeros_like(y)
    iC=np.zeros_like(y)
    iG[:,0]=1
    for i,x in enumerate(iC):
        x[ctf[i]]=1
    iG=iG.flatten()
    iC=iC.flatten()
    y=y.flatten()
    T1 = toeplitz(iG, np.zeros(tlength))
    T2 = toeplitz(iC, np.zeros(tlength))
    A_k = np.hstack([T1, T2])
    # try:
    model = LinearRegression(positive=True,fit_intercept=False)
    # model = ElasticNetCV(fit_intercept=False, positive=True)
    model.fit(A_k, y)
    irfs = model.coef_.reshape(-1,2,order='F')
    # except:
    #     irfs = np.ones([tlength,2])*np.nan
    return irfs

def uniqueSplit(array):
    shuff = np.random.permutation(array)
    split = int(len(array)/ 2)
    return shuff[:split], shuff[split:]

def pullAnmHH(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, alignTime = 'goCues',
                 dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events', factors = [4,2],
              enforceMinTrial = False, returnErrors = False, binary = False, anmsMinTrial = None):
    
    mT_preTType1 = []
    mT_preTType2 = []
    mTs = [mT_preTType1, mT_preTType2]
    fr = frs[ak]
    
    tsepKey = 'dsTsep'
    factor = factors[ak]
    keepRois = keepRoisDS[anat][anm]
        
    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'lastLick':
        alignKey = 'll-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]
        
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays

    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]  
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    FRm = tkr[:,2]==5
    FLm = tkr[:,2]==6

    if enforceMinTrial:
        nTriToDraw = np.floor(anmsMinTrial[anat][anm]/2).astype(int)

    if not returnErrors:
        FRts = uniqueSplit(np.where(FRm)[0])
        CRts = uniqueSplit(np.where(CRm)[0])
        CLts = uniqueSplit(np.where(CLm)[0])
        if np.sum(FLm)>10:
            FLts = uniqueSplit(np.where(FLm)[0])
            hasFL = True
        else:
            hasFL = False
            
        for i in range(2):
            crt = CRts[i]
            clt = CLts[i]
            frt = FRts[i]
            if enforceMinTrial:
                CRt = np.random.choice(crt, nTriToDraw, replace = False)
                FRt = np.random.choice(frt, nTriToDraw)#, replace = False)
                CLt = np.random.choice(clt, nTriToDraw, replace = False)
            else:
                CRt = crt
                FRt = frt
                CLt = clt
                
            if binary:
                CR = np.nanmean(psTsep[CRt,:,:]>0, axis = 0)
                FR = np.nanmean(psTsep[FRt,:,:]>0, axis = 0)
                CL = np.nanmean(psTsep[CLt,:,:]>0, axis = 0)
            else:
                CR = np.nanmean(psTsep[CRt,:,:], axis = 0)
                FR = np.nanmean(psTsep[FRt,:,:], axis = 0)
                CL = np.nanmean(psTsep[CLt,:,:], axis = 0)
    
            if hasFL:
                flt = FLts[i]
                if enforceMinTrial:
                    FLt = np.random.choice(flt, nTriToDraw)#, replace = False)
                else:
                    FLt = flt
                if binary:
                    FL = np.nanmean(psTsep[FLt,:,:]>0, axis = 0)
                else:
                    FL = np.nanmean(psTsep[FLt,:,:], axis = 0)
            else:
                FL = np.ones_like(CL)*np.nan
    
            tsep3 = np.array([CR, CL, FR, FL])

            mTs[i].append(tsep3)
    else:
        FRm = tkr[:,2]==5
        MEm = np.logical_and(FRm,me2)
        DEm = np.logical_and(FRm,de2)

        if np.sum(DEm)>2:
            hasDE = True
        else:
            hasDE = False

        if np.sum(MEm)>2:
            hasME = True
        else:
            hasME = False

        CRts = uniqueSplit(np.where(CRm)[0])
        CLts = uniqueSplit(np.where(CLm)[0])
        
        if hasDE:
            DEts = uniqueSplit(np.where(DEm)[0])
        if hasME:
            MEts = uniqueSplit(np.where(MEm)[0])
        
        for i in range(2):
            if hasDE:
                det = DEts[i]
            if hasME:
                met = MEts[i]
            crt = CRts[i]
            clt = CLts[i]
            
            if enforceMinTrial:
                CRt = np.random.choice(crt, nTriToDraw, replace = False)
                CLt = np.random.choice(clt, nTriToDraw, replace = False)
                if hasDE:
                    DEt = np.random.choice(det, nTriToDraw)#, replace = False)
                if hasME:
                    MEt = np.random.choice(met, nTriToDraw)#, replace = False)
            else:
                CRt = crt
                CLt = clt
                if hasDE:
                    DEt = det
                if hasME:
                    MEt = met

            if binary:
                CR = np.nanmean(psTsep[CRt,:,:]>0, axis = 0)
                CL = np.nanmean(psTsep[CLt,:,:]>0, axis = 0)
            else:
                CR = np.nanmean(psTsep[CRt,:,:], axis = 0)
                CL = np.nanmean(psTsep[CLt,:,:], axis = 0)
            
            if hasDE:
                if binary:
                    DE = np.nanmean(psTsep[DEt,:,:]>0, axis = 0)
                else:
                    DE = np.nanmean(psTsep[DEt,:,:], axis = 0)
            else:
                DE = np.ones_like(CR)*np.nan
            if hasME:
                if binary:
                    ME = np.nanmean(psTsep[MEt,:,:]>0, axis = 0)
                else:
                    ME = np.nanmean(psTsep[MEt,:,:], axis = 0)
            else:
                ME = np.ones_like(CR)*np.nan
            
            tsep3 = np.array([CR, CL, ME, DE])#, FR, FL])
            
            mTs[i].append(tsep3)
                
    mT_preTType1 = np.dstack(mT_preTType1)
    mT_preTType2 = np.dstack(mT_preTType2)
    return mT_preTType1, mT_preTType2

def gather_anmShuff_megaTsepHH(allAnmParams, anat, ak, mod, phases, frs,  keepRoisDS, dffKey, alignTime, anmsToUse, nRois, anmIdxs, enforceMinTrial = False,
                               binary = False, returnErrors = False, anmsMinTrial=None):

    allMega1 = {phase: [] for phase in phases}
    allMega2 = {phase: [] for phase in phases}
    LUT = {phase: [] for phase in phases}
    
    allAnmHH1 = {}
    allAnmHH2 = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmHH1[anm] = {}
        allAnmHH2[anm] = {}
        for phase in phases:
            anm1, anm2 = pullAnmHH(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, alignTime = alignTime,
                     dffKey = dffKey, enforceMinTrial = enforceMinTrial, returnErrors = returnErrors,
                                            binary = binary, anmsMinTrial=anmsMinTrial)
            allAnmHH1[anm][phase] = anm1 # trial, roi, time
            allAnmHH2[anm][phase] = anm2
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmHH1[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega1[phase].append(allAnmHH1[anm][phase][:,roi,:])
                    allMega2[phase].append(allAnmHH2[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega1[phase] = np.dstack(allMega1[phase]) #trial, time, roi
        allMega2[phase] = np.dstack(allMega2[phase])
    return allMega1, allMega2

def pullAnmTsep(allAnmParams, anm, anat, ak, phase, frs, keepRoisDS, alignTime = 'goCues',
                 dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events', factors = [4,2], returnErrors = False, binary = False,
               enforceMinTrial = False, anmsMinTrial = None):
    
    mT_preTType = []
    if anat == 'dendrites':
        ak = 0
    else:
        ak = 1
        
    fr = frs[ak]
    tsepKey = 'dsTsep'
    factor = factors[ak]
    keepRois = keepRoisDS[anat][anm]
        
    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'lastLick':
        alignKey = 'll-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]
        
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays
        
    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]
        
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1

    if enforceMinTrial:
        nTriToDraw = anmsMinTrial[anat][anm]
        
    if not returnErrors:

        FRm = tkr[:,2]==5
        FLm = tkr[:,2]==6

        if enforceMinTrial:
            CRt = np.random.choice(np.where(CRm)[0], nTriToDraw, replace = False)
            CLt = np.random.choice(np.where(CLm)[0], nTriToDraw, replace = False)
            FRt = np.random.choice(np.where(FRm)[0], nTriToDraw)
        else:
            CRt = np.where(CRm)[0]
            CLt = np.where(CLm)[0]
            FRt = np.where(FRm)[0]

        if binary:
            CR = np.nanmean(psTsep[CRt,:,:]>0, axis = 0)
            CL = np.nanmean(psTsep[CLt,:,:]>0, axis = 0)
            FR = np.nanmean(psTsep[FRt,:,:]>0, axis = 0)
        else:
            CR = np.nanmean(psTsep[CRt,:,:], axis = 0)
            CL = np.nanmean(psTsep[CLt,:,:], axis = 0)
            FR = np.nanmean(psTsep[FRt,:,:], axis = 0)

        if np.sum(FLm)>=5:
            if enforceMinTrial:
                FLt = np.random.choice(np.where(FLm)[0], nTriToDraw)#, replace = False)
            else:
                FLt = np.where(FLm)[0]
            
            if binary:
                FL = np.nanmean(psTsep[FLt,:,:]>0, axis = 0)
            else:
                FL = np.nanmean(psTsep[FLt,:,:], axis = 0)
        else:
            FL = np.ones_like(CL)*np.nan

        tsep3 = np.array([CR, CL, FR, FL])

    else:

        rMEm = np.logical_and(tkr[:,2]==5, me2)
        rDEm = np.logical_and(tkr[:,2]==5, de2)

        if binary:
            CR = np.nanmean(psTsep[CRm,:,:]>0, axis = 0)
            CL = np.nanmean(psTsep[CLm,:,:]>0, axis = 0)
        else:
            CR = np.nanmean(psTsep[CRm,:,:], axis = 0)
            CL = np.nanmean(psTsep[CLm,:,:], axis = 0)

        if np.sum(rMEm)>=5:
            if binary:
                ME = np.nanmean(psTsep[rMEm,:,:]>0, axis = 0)
            else:
                ME = np.nanmean(psTsep[rMEm,:,:], axis = 0)
        else:
            ME = np.ones_like(CL)*np.nan

        if np.sum(rDEm)>=5:
            if binary:
                DE = np.nanmean(psTsep[rDEm,:,:]>0, axis = 0)
            else:
                DE = np.nanmean(psTsep[rDEm,:,:], axis = 0)
        else:
            DE = np.ones_like(CL)*np.nan

        tsep3 = np.array([CR, CL, ME, DE])
        
    mT_preTType.append(tsep3)

    mT_preTType = np.dstack(mT_preTType)
    return mT_preTType

def pullAnmHH_ALT(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, keepRoisFull, alignTime = 'goCues',
                 dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events', factors = [4,2],
              useDS = True, enforceMinTrial = False, nTriToDraw = 20, binary = False, giveAll = False):
    
    mT_preTType1 = []
    mT_preTType2 = []
    mTs = [mT_preTType1, mT_preTType2]
    fr = frs[ak]
    if useDS:
        tsepKey = 'dsTsep'
        factor = factors[ak]
        keepRois = keepRoisDS[anat][anm]
    else:
        tsepKey = 'fullTsep'
        factor = 1
        keepRois = keepRoisFull[anat][anm]
        
    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'secondLick':
        alignKey = 'sl-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]
        
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    ar = allAnmParams[anm]['afterRight']
    al_noME = allAnmParams[anm]['MER-L_noME']
    al_ME = allAnmParams[anm]['MER-L_ME']
    an = allAnmParams[anm]['afterNone']
    slN = allAnmParams[anm]['sl_nHits']
    slTracking = allAnmParams[anm]['slSideTracking']
    
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays

    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]  
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    FRm = tkr[:,2]==5
    MEm = np.logical_and(FRm,me2)
    DEm = np.logical_and(FRm,de2)
    ARm = ar[fullMask]
    AL_noMEm = al_noME[fullMask]
    AL_MEm = al_ME[fullMask]
    ANm = an[fullMask]
    noSecondHitMask = np.logical_and(slN[fullMask]==0, MEm)
    hitMask = np.logical_and(slN[fullMask]>0, MEm)
    slAttmptR = np.logical_and(MEm, slTracking[fullMask]==0)
    slAttmptL = np.logical_and(MEm, slTracking[fullMask]==1)

    hasAR = np.sum(ARm)>2
    hasAL_noME = np.sum(AL_noMEm)>2
    hasAL_ME = np.sum(AL_MEm)>2
    hasAN = np.sum(ANm)>2
    hasDE = np.sum(DEm)>2
    hasME = np.sum(MEm)>2
    hasNSH = np.sum(noSecondHitMask)>2
    hasHit = np.sum(hitMask)>2
    
    CRts = uniqueSplit(np.where(CRm)[0])
    CLts = uniqueSplit(np.where(CLm)[0])
    SLRts = uniqueSplit(np.where(slAttmptR)[0])
    SLLts = uniqueSplit(np.where(slAttmptL)[0])

    if hasDE:
        DEts = uniqueSplit(np.where(DEm)[0])
    if hasME:
        MEts = uniqueSplit(np.where(MEm)[0])
    if hasAR:
        ARts = uniqueSplit(np.where(ARm)[0])
    if hasAL_noME:
        AL_noMEts = uniqueSplit(np.where(AL_noMEm)[0])
    if hasAL_ME:
        AL_MEts = uniqueSplit(np.where(AL_MEm)[0])
    if hasAN:
        ANts = uniqueSplit(np.where(ANm)[0])
    if hasNSH:
        NSHts = uniqueSplit(np.where(noSecondHitMask)[0])
    if hasHit:
        Hit_ts = uniqueSplit(np.where(hitMask)[0])

    for i in range(2):
        if hasDE:
            det = DEts[i]
        if hasME:
            met = MEts[i]
        if hasAR:
            art = ARts[i]
        if hasAL_noME:
            al_noMEt = AL_noMEts[i]
        if hasAL_ME:
            al_MEt = AL_MEts[i]
        if hasAN:
            ant = ANts[i]
        if hasNSH:
            nsh_t = NSHts[i]
        if hasHit:
            hit_t = Hit_ts[i]
            
        crt = CRts[i]
        clt = CLts[i]
        slrt = SLRts[i]
        sllt = SLLts[i]
        
        if enforceMinTrial:
            CRt = np.random.choice(crt, nTriToDraw)
            CLt = np.random.choice(clt, nTriToDraw)
            SLRt = np.random.choice(slrt, nTriToDraw)
            SLLt = np.random.choice(sllt, nTriToDraw)
            if hasDE:
                DEt = np.random.choice(det, nTriToDraw)
            if hasME:
                MEt = np.random.choice(met, nTriToDraw)
            if hasAR:
                ARt = np.random.choice(art, nTriToDraw)
            if hasAL_ME:
                AL_MEt = np.random.choice(al_MEt, nTriToDraw)
            if hasAL_noME:
                AL_noMEt = np.random.choice(al_noMEt, nTriToDraw)
            if hasAN:
                ANt = np.random.choice(ant, nTriToDraw)
            if hasNSH:
                NSHt = np.random.choice(nsh_t, nTriToDraw)
            if hasHit:
                hitT = np.random.choice(hit_t, nTriToDraw)
 
        else:
            CRt = crt
            CLt = clt
            SLRt = slrt
            SLLt = sllt
            if hasDE:
                DEt = det
            if hasME:
                MEt = met
            if hasAR:
                ARt = art
            if hasAL_ME:
                AL_MEt = al_MEt
            if hasAL_noME:
                AL_noMEt = al_noMEt
            if hasAN:
                ANt = ant
            if hasNSH:
                NSHt = nsh_t
            if hasHit:
                hitT = hit_t

        if binary:
            CR = np.nanmean(psTsep[CRt,:,:]>0, axis = 0)
        else:
            CR = np.nanmean(psTsep[CRt,:,:], axis = 0)

        if binary:
            SLR = np.nanmean(psTsep[SLRt,:,:]>0, axis = 0)
        else:
            SLR = np.nanmean(psTsep[SLRt,:,:], axis = 0)
            
        if binary:
            SLL = np.nanmean(psTsep[SLLt,:,:]>0, axis = 0)
        else:
            SLL = np.nanmean(psTsep[SLLt,:,:], axis = 0)

        if binary:
            CL = np.nanmean(psTsep[CLt,:,:]>0, axis = 0)
        else:
            CL = np.nanmean(psTsep[CLt,:,:], axis = 0)

        if hasNSH:
            if binary:
                NSH = np.nanmean(psTsep[NSHt,:,:]>0, axis = 0)
            else:
                NSH = np.nanmean(psTsep[NSHt,:,:], axis = 0)
        else:
            NSH = np.ones_like(CR)*np.nan

        if hasHit:
            if binary:
                HIT = np.nanmean(psTsep[hitT,:,:]>0, axis = 0)
            else:
                HIT = np.nanmean(psTsep[hitT,:,:], axis = 0)
        else:
            HIT = np.ones_like(CR)*np.nan
        
        if hasDE:
            if binary:
                DE = np.nanmean(psTsep[DEt,:,:]>0, axis = 0)
            else:
                DE = np.nanmean(psTsep[DEt,:,:], axis = 0)
        else:
            DE = np.ones_like(CR)*np.nan
        
        if hasME:
            if binary:
                ME = np.nanmean(psTsep[MEt,:,:]>0, axis = 0)
            else:
                ME = np.nanmean(psTsep[MEt,:,:], axis = 0)
        else:
            ME = np.ones_like(CR)*np.nan
            
        if hasAR:
            if binary:
                AR = np.nanmean(psTsep[ARt,:,:]>0, axis = 0)
            else:
                AR = np.nanmean(psTsep[ARt,:,:], axis = 0)
        else:
            AR = np.ones_like(CR)*np.nan

        if hasAL_ME:
            if binary:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:]>0, axis = 0)
            else:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:], axis = 0)
        else:
            AL_ME = np.ones_like(CR)*np.nan

        if hasAL_noME:
            if binary:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:]>0, axis = 0)
            else:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:], axis = 0)
        else:
            AL_noME = np.ones_like(CR)*np.nan

        if hasAN:
            if binary:
                AN = np.nanmean(psTsep[ANt,:,:]>0, axis = 0)
            else:
                AN = np.nanmean(psTsep[ANt,:,:], axis = 0)
        else:
            AN = np.ones_like(CR)*np.nan

        if giveAll:
            #tsep3 = np.array([CR,AL_noME,AL_ME,AR,AN,DE,CL,NSH,HIT])
            tsep3 = np.array([CR,AL_noME,AL_ME,AR,ME, AN,DE,CL,SLR,SLL])
        else:
            #tsep3 = np.array([CR, ME, AR, AL])
            tsep3 = np.array([CR, AL_noME, AL_ME, AR])
        
        mTs[i].append(tsep3)
                
    mT_preTType1 = np.dstack(mT_preTType1)
    mT_preTType2 = np.dstack(mT_preTType2)
    return mT_preTType1, mT_preTType2

def gather_anmShuff_megaTsepHH_ALT(allAnmParams, anat, ak, mod, phases, frs, \
    keepRoisDS, keepRoisFull, dffKey, alignTime, anmsToUse, nRois, anmIdxs, \
    enforceMinTrial = False, nTriToDraw = 50,
    binary = False, useDS=True, giveAll = False):

    allMega1 = {phase: [] for phase in phases}
    allMega2 = {phase: [] for phase in phases}
    LUT = {phase: [] for phase in phases}
    
    allAnmHH1 = {}
    allAnmHH2 = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmHH1[anm] = {}
        allAnmHH2[anm] = {}
        for phase in phases:
            anm1, anm2 = pullAnmHH_ALT(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, keepRoisFull, alignTime = alignTime,
                     dffKey = dffKey, useDS = useDS, enforceMinTrial = enforceMinTrial, 
                                               nTriToDraw = nTriToDraw, binary = binary, giveAll = giveAll)
            allAnmHH1[anm][phase] = anm1 # trial, roi, time
            allAnmHH2[anm][phase] = anm2
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmHH1[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega1[phase].append(allAnmHH1[anm][phase][:,roi,:])
                    allMega2[phase].append(allAnmHH2[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega1[phase] = np.dstack(allMega1[phase]) #trial, time, roi
        allMega2[phase] = np.dstack(allMega2[phase])
    return allMega1, allMega2

def pullAnmTsep_ALT(allAnmParams, anm, anat, phase, frs, \
    keepRoisDS, keepRoisFull, alignTime = 'goCues',
    dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events', factors = [4,2], useDS = True, binary = False, giveAll = False):
    
    mT_preTType = []
    if anat == 'dendrites':
        ak = 0
    else:
        ak = 1
        
    fr = frs[ak]
    if useDS:
        tsepKey = 'dsTsep'
        factor = factors[ak]
        keepRois = keepRoisDS[anat][anm]
    else:
        tsepKey = 'fullTsep'
        factor = 1
        keepRois = keepRoisFull[anat][anm]
        
    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'secondLick':
        alignKey = 'sl-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]
        
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    ar = allAnmParams[anm]['afterRight']
    al_noME = allAnmParams[anm]['MER-L_noME']
    al_ME = allAnmParams[anm]['MER-L_ME']
    an = allAnmParams[anm]['afterNone']
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        
        anatMask = ~dendDays
        
    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]
        
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    rMEm = np.logical_and(tkr[:,2]==5, me2)
    rDEm = np.logical_and(tkr[:,2]==5, de2)
    ARm = ar[fullMask]
    ALnoMEm = al_noME[fullMask]
    AL_MEm = al_ME[fullMask]
    ANm = an[fullMask]

    if binary:
        CR = np.nanmean(psTsep[CRm,:,:]>0, axis = 0)
        CL = np.nanmean(psTsep[CLm,:,:]>0, axis = 0)
    else:
        CR = np.nanmean(psTsep[CRm,:,:], axis = 0)
        CL = np.nanmean(psTsep[CLm,:,:], axis = 0)

    if np.sum(rDEm)>=5:
        if binary:
            DE = np.nanmean(psTsep[rDEm,:,:]>0, axis = 0)
        else:
            DE = np.nanmean(psTsep[rDEm,:,:], axis = 0)
    else:
        DE = np.ones_like(CR)*np.nan

    if np.sum(rMEm)>=5:
        if binary:
            ME = np.nanmean(psTsep[rMEm,:,:]>0, axis = 0)
        else:
            ME = np.nanmean(psTsep[rMEm,:,:], axis = 0)
    else:
        ME = np.ones_like(CR)*np.nan

    if np.sum(ARm)>=5:
        if binary:
            AR = np.nanmean(psTsep[ARm,:,:]>0, axis = 0)
        else:
            AR = np.nanmean(psTsep[ARm,:,:], axis = 0)
    else:
        AR = np.ones_like(CR)*np.nan

    if np.sum(ALnoMEm)>=5:
        if binary:
            AL_noME = np.nanmean(psTsep[ALnoMEm,:,:]>0, axis = 0)
        else:
            AL_noME = np.nanmean(psTsep[ALnoMEm,:,:], axis = 0)
    else:
        AL_noME = np.ones_like(CR)*np.nan

    if np.sum(AL_MEm)>=5:
        if binary:
            AL_ME = np.nanmean(psTsep[AL_MEm,:,:]>0, axis = 0)
        else:
            AL_ME = np.nanmean(psTsep[AL_MEm,:,:], axis = 0)
    else:
        AL_ME = np.ones_like(CR)*np.nan

    if np.sum(ANm)>=5:
        if binary:
            AN = np.nanmean(psTsep[ANm,:,:]>0, axis = 0)
        else:
            AN = np.nanmean(psTsep[ANm,:,:], axis = 0)
    else:
        AN = np.ones_like(CR)*np.nan

    if giveAll:
        tsep3 = np.array([CR, AL_noME, AL_ME, AR, ME, AN, DE, CL])
    else:
        #tsep3 = np.array([CR, ME, AR, AL])
        tsep3 = np.array([CR, AL_noME, AL_ME, AR])
        
    mT_preTType.append(tsep3)

    mT_preTType = np.dstack(mT_preTType)
    return mT_preTType

def gather_anmShuff_megaTsepFull_ALT(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, keepRoisFull, dffKey, \
    alignTime, anmsToUse, nRois, anmIdxs, enforceMinTrial = False, nTriToDraw = 50,\
    binary = False, useDS=True, giveAll = False):

    allMega = {phase: [] for phase in phases}
    
    allAnmTsep = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmTsep[anm] = {}
        for phase in phases:
            anmTsep = pullAnmTsep_ALT(allAnmParams, anm, anat, phase, frs, \
                    keepRoisDS, keepRoisFull, alignTime = alignTime,\
                     dffKey = dffKey, useDS = useDS, binary = binary, giveAll = giveAll)
            allAnmTsep[anm][phase] = anmTsep # trial, roi, time
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmTsep[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega[phase].append(allAnmTsep[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega[phase] = np.dstack(allMega[phase]) #trial, time, roi
    return allMega

def gather_anmShuff_megaTsepFull(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, dffKey, alignTime, anmsToUse, nRois, anmIdxs, enforceMinTrial = False,
                               binary = False, returnErrors = False, anmsMinTrial=None):

    allMega = {phase: [] for phase in phases}
    
    allAnmTsep = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmTsep[anm] = {}
        for phase in phases:
            anmTsep = pullAnmTsep(allAnmParams, anm, anat, ak, phase, frs, keepRoisDS, alignTime = alignTime,
                     dffKey = dffKey, returnErrors = returnErrors, binary = binary, enforceMinTrial = enforceMinTrial,
                                 anmsMinTrial=anmsMinTrial)
            allAnmTsep[anm][phase] = anmTsep # trial, roi, time
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmTsep[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega[phase].append(allAnmTsep[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega[phase] = np.dstack(allMega[phase]) #trial, time, roi
    return allMega

def pullAnmTsepComputeSMOXval(allAnmParams, anm, anat, phase, frs, alignTime, factors, minTrial, code, keepRoisDS,  keepRoisFull, \
    ortho=False, orthoOrder=False, fixedCDs=False, dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events',useDS = True, binary = False, returnCDs = False):
    """
    Compute cross-validated (leave-one-out) coding-direction (CD) projections of
    single-trial, time-separated activity for one animal, one anatomy, one task
    phase, and one alignment.

    For each of the four trial types (CR, CL, IR, IL) the CD is the weighted
    contrast of the condition mean traces (weights given by `code`), where the
    held-out trial's own contribution is replaced by the leave-one-out mean of
    its condition. Each single trial is then projected (elementwise multiplied)
    onto that CD and averaged across held-out trials. Optionally the CD is
    Gram-Schmidt orthogonalized (`ortho`/`orthoOrder`) against other CDs, and/or
    collapsed to fixed time windows (`fixedCDs`). If any condition has
    <= `minTrial` trials the output for this animal is all-NaN.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anm = 'JS001'             # animal id key into allAnmParams / keepRois dicts
    anat = 'dendrites'        # 'dendrites' -> frs/factors index 0; anything else -> index 1
    phase = 'pre'             # task phase key into allAnmParams[anm]['shiftMask']
    frs = [frDend, frSoma]    # per-anatomy frame rate (Hz); index 0 dendrites, 1 otherwise
    alignTime = 'goCues'      # 'goCues' | 'cTimes' | 'lastLick'
    factors = [10, 5]         # per-anatomy temporal downsample factor; index matches frs
    minTrial = 5              # min trials per condition; <= this -> all-NaN output
    code = [1, -1, 1, -1]     # 4 CD weights over [CR, CL, IR, IL] (e.g. s/m/o contrasts)
    keepRoisDS = keepRoisDS   # dict[anat][anm] -> ROI indices (downsampled data)
    keepRoisFull = keepRoisFull  # dict[anat][anm] -> ROI indices (full-res data)
    ortho = False             # False -> no orthogonalization; or list of other code vectors to orthogonalize against
    orthoOrder = False        #     > only if ortho set: np.array giving Gram-Schmidt order; the target CD is the entry equal to 0
    fixedCDs = False          # False -> time-resolved CD; or list of [start, end] windows (s, rel. to align) to average the CD over
    dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events'  # key selecting the dff/events trace in allAnmParams[anm]
    useDS = True              # True -> use 'dsTsep' + factors + keepRoisDS; False -> 'fullTsep' + factor 1 + keepRoisFull
    binary = False            # accepted for call-signature compatibility; not used in this function
    returnCDs = False         # False -> return projections only; True -> also return the trial-averaged CD axes
    importlib.reload(js_manuscript)
    mT_preTType = js_manuscript.pullAnmTsepComputeSMOXval(
        allAnmParams, anm=anm, anat=anat, phase=phase, frs=frs, alignTime=alignTime,
        factors=factors, minTrial=minTrial, code=code, keepRoisDS=keepRoisDS,
        keepRoisFull=keepRoisFull, ortho=ortho, orthoOrder=orthoOrder,
        fixedCDs=fixedCDs, dffKey=dffKey, useDS=useDS, binary=binary, returnCDs=returnCDs)
    # when returnCDs=True: mT_preTType, cdAxes = js_manuscript.pullAnmTsepComputeSMOXval(..., returnCDs=returnCDs)

    Parameters
    ----------
    allAnmParams : dict
        Per-animal parameter/data dict. `allAnmParams[anm]` provides the
        time-separated traces (`[dffKey][tsepKey][anat][alignKey]`), plus
        'tracker', 'dendDays', 'cTimes', 'goCues', 'shiftMask', 'de', 'me'.
    anm : str
        Animal id key into `allAnmParams` and the `keepRois` dicts.
    anat : str
        Anatomy key. 'dendrites' selects frame-rate/factor index 0 and the
        dendrite-day mask; any other value selects index 1 and its complement.
    phase : str
        Task phase key into `allAnmParams[anm]['shiftMask']`.
    frs : list of float
        Per-anatomy imaging frame rate (Hz); index 0 for dendrites, 1 otherwise.
    alignTime : str
        Trial-alignment event, one of 'goCues', 'cTimes', 'lastLick', mapped to
        the 'go-psTsep' / 'ct-psTsep' / 'll-psTsep' trace keys.
    factors : list of int
        Per-anatomy temporal downsample factor (used only when `useDS` is True),
        indexed the same way as `frs`.
    minTrial : int
        Minimum trials required in each of the four conditions; if the smallest
        condition count is <= `minTrial` the output is all-NaN.
    code : list of float
        Length-4 CD weight vector over [CR, CL, IR, IL] defining the contrast
        (e.g. sensory [1,-1,1,-1], motor [1,-1,-1,1], outcome [1,1,-1,-1]).
    keepRoisDS : dict
        Nested dict `keepRoisDS[anat][anm]` of ROI indices to keep for the
        downsampled traces (used when `useDS` is True).
    keepRoisFull : dict
        Nested dict `keepRoisFull[anat][anm]` of ROI indices to keep for the
        full-resolution traces (used when `useDS` is False).
    ortho : bool or list
        False for no orthogonalization; otherwise a list of additional `code`
        vectors whose CDs the target CD is Gram-Schmidt orthogonalized against.
    orthoOrder : bool or np.ndarray
        Only used when `ortho` is set. Array giving the orthogonalization order
        of the stacked CDs; the target CD is the position whose value is 0.
    fixedCDs : bool or list
        False for a time-resolved CD; otherwise a list of `[start, end]` windows
        (seconds relative to the alignment event) over which each CD is averaged
        and broadcast across time before projection.
    dffKey : str
        Key selecting the dff/deconvolved-events trace within `allAnmParams[anm]`.
    useDS : bool
        True to use the downsampled traces ('dsTsep', `factors`, `keepRoisDS`);
        False to use full-resolution traces ('fullTsep', factor 1, `keepRoisFull`).
    binary : bool
        Accepted for call-signature compatibility with related functions; not
        referenced in this function.
    returnCDs : bool
        False returns only the projection tensor. True additionally returns
        `cdAxes`, the trial-averaged coding-direction axes (before projection)
        that produced those projections, as a second output.

    Returns
    -------
    mT_preTType : np.ndarray
        Array of shape (4, nRoi, nTime) holding the CD projections for the four
        trial types [CR, CL, IR, IL]; all-NaN if any condition has too few trials.
    cdAxes : np.ndarray
        Only returned when `returnCDs` is True. Array of shape (4, nRoi, nTime)
        holding the trial-averaged coding-direction axis (the `cd` vector, after
        any orthogonalization/fixed-window collapse) for each of the four trial
        types [CR, CL, IR, IL]; all-NaN if any condition has too few trials. The
        stored projections equal these axes multiplied elementwise by each
        single trial and averaged over held-out trials.
    """


    mT_preTType = []
    if anat == 'dendrites':
        ak = 0
    else:
        ak = 1

    fr = frs[ak]
    if useDS:
        tsepKey = 'dsTsep'
        factor = factors[ak]
        keepRois = keepRoisDS[anat][anm]
    else:
        tsepKey = 'fullTsep'
        factor = 1
        keepRois = keepRoisFull[anat][anm]

    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'lastLick':
        alignKey = 'll-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]

    if np.any(fixedCDs):
        nFrames=tsep.shape[2]
        hF = int(nFrames/2)
        fixedCDsF=[]
        for fcd in fixedCDs:
            fCDStart= hF+round(fcd[0]*fr/factor)
            fCDend= hF+round(fcd[1]*fr/factor) 
            fixedCDsF.append([fCDStart,fCDend])
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays
        
    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]
        
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]  
    
    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    IRm = tkr[:,2]==5
    ILm = tkr[:,2]==6

    nTrials = [np.sum(CRm),np.sum(CLm),np.sum(IRm),np.sum(ILm)]

    # rMEm = np.logical_and(tkr[:,2]==5, me2)
    CR = np.nanmean(psTsep[CRm,:,:], axis = 0)

    if np.min(nTrials)<=minTrial:
        CR = np.ones_like(CR)*np.nan
        CL = np.ones_like(CR)*np.nan
        IR = np.ones_like(CR)*np.nan
        IL = np.ones_like(CR)*np.nan
        if returnCDs:
            CRcd = np.ones_like(CR)*np.nan
            CLcd = np.ones_like(CR)*np.nan
            IRcd = np.ones_like(CR)*np.nan
            ILcd = np.ones_like(CR)*np.nan
    else:
        CRm = CRm.nonzero()[0]
        CLm = CLm.nonzero()[0]
        IRm = IRm.nonzero()[0]
        ILm = ILm.nonzero()[0]

        RC=[]
        LC=[]
        RI=[]
        LI=[]
        if np.any(fixedCDs):
            for fcd in fixedCDsF:
                RCi = np.nanmean(psTsep[CRm,:,:], axis = 0)
                RCi = np.nanmean(RCi[:,fcd[0]:fcd[1]],axis=1)
                RC.append(RCi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                LCi = np.nanmean(psTsep[CLm,:,:], axis = 0)
                LCi = np.nanmean(LCi[:,fcd[0]:fcd[1]],axis=1)
                LC.append(LCi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                RIi = np.nanmean(psTsep[IRm,:,:], axis = 0)
                RIi = np.nanmean(RIi[:,fcd[0]:fcd[1]],axis=1)
                RI.append(RIi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                LIi = np.nanmean(psTsep[ILm,:,:], axis = 0)
                LIi = np.nanmean(LIi[:,fcd[0]:fcd[1]],axis=1)
                LI.append(LIi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))          
        else:
            for i in range(3):
                RC.append(np.nanmean(psTsep[CRm,:,:], axis = 0))
                LC.append(np.nanmean(psTsep[CLm,:,:], axis = 0))
                RI.append(np.nanmean(psTsep[IRm,:,:], axis = 0))
                LI.append(np.nanmean(psTsep[ILm,:,:], axis = 0))
        
        CR=[]
        CRcd=[]
        for trial in CRm:
            cdtrials = np.setxor1d(trial,CRm)
            x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
            if np.any(fixedCDs):
                x = np.nanmean(x[:,fixedCDsF[0][0]:fixedCDsF[0][1]],axis=1)
                x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            cd = (x*code[0]+LC[0]*code[1]+RI[0]*code[2]+LI[0]*code[3])/2
            if np.any(ortho):
                ocds = []
                for i,oc in enumerate(ortho):      
                    if np.any(fixedCDs):
                        x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
                        x = np.nanmean(x[:,fixedCDsF[i+1][0]:fixedCDsF[i+1][1]],axis=1)
                        x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    ocds.append((x*oc[0]+LC[i+1]*oc[1]+RI[i+1]*oc[2]+LI[i+1]*oc[3])/2)
                acds = np.dstack([cd]+ocds)
                acds = acds[:,:, orthoOrder]
                cdIdx = (orthoOrder==0).nonzero()[0][0]
                if np.min(acds.shape)>2:
                    cd=[]
                    for tp in range(acds.shape[1]):
                        [ocds,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                        cd.append(ocds[:,cdIdx]*r[cdIdx,cdIdx])
                    cd = np.vstack(cd).T
            if returnCDs:
                CRcd.append(cd)
            CR.append(cd*psTsep[trial,:,:])
        CR=np.nanmean(np.dstack(CR),axis=2)
        if returnCDs:
            CRcd=np.nanmean(np.dstack(CRcd),axis=2)

        CL=[]
        CLcd=[]
        for trial in CLm:
            cdtrials = np.setxor1d(trial,CLm)
            x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
            if np.any(fixedCDs):
                x = np.nanmean(x[:,fixedCDsF[0][0]:fixedCDsF[0][1]],axis=1)
                x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            cd = (RC[0]*code[0]+x*code[1]+RI[0]*code[2]+LI[0]*code[3])/2
            if np.any(ortho):
                ocds = []
                for i,oc in enumerate(ortho):      
                    if np.any(fixedCDs):
                        x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
                        x = np.nanmean(x[:,fixedCDsF[i+1][0]:fixedCDsF[i+1][1]],axis=1)
                        x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    ocds.append((RC[i+1]*oc[0]+x*oc[1]+RI[i+1]*oc[2]+LI[i+1]*oc[3])/2)
                acds = np.dstack([cd]+ocds)
                acds = acds[:,:, orthoOrder]
                cdIdx = (orthoOrder==0).nonzero()[0][0]
                if np.min(acds.shape)>2:
                    cd=[]
                    for tp in range(acds.shape[1]):
                        [ocds,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                        cd.append(ocds[:,cdIdx]*r[cdIdx,cdIdx])
                    cd = np.vstack(cd).T
            if returnCDs:
                CLcd.append(cd)
            CL.append(cd*psTsep[trial,:,:])
        CL=np.nanmean(np.dstack(CL),axis=2)
        if returnCDs:
            CLcd=np.nanmean(np.dstack(CLcd),axis=2)

        IR=[]
        IRcd=[]
        for trial in IRm:
            cdtrials = np.setxor1d(trial,IRm)
            x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
            if np.any(fixedCDs):
                x = np.nanmean(x[:,fixedCDsF[0][0]:fixedCDsF[0][1]],axis=1)
                x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            cd = (RC[0]*code[0]+LC[0]*code[1]+x*code[2]+LI[0]*code[3])/2
            if np.any(ortho):
                ocds = []
                for i,oc in enumerate(ortho):      
                    if np.any(fixedCDs):
                        x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
                        x = np.nanmean(x[:,fixedCDsF[i+1][0]:fixedCDsF[i+1][1]],axis=1)
                        x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    ocds.append((RC[i+1]*oc[0]+LC[i+1]*oc[1]+x*oc[2]+LI[i+1]*oc[3])/2)
                acds = np.dstack([cd]+ocds)
                acds = acds[:,:, orthoOrder]
                cdIdx = (orthoOrder==0).nonzero()[0][0]
                if np.min(acds.shape)>2:
                    cd=[]
                    for tp in range(acds.shape[1]):
                        [ocds,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                        cd.append(ocds[:,cdIdx]*r[cdIdx,cdIdx])
                    cd = np.vstack(cd).T
            if returnCDs:
                IRcd.append(cd)
            IR.append(cd*psTsep[trial,:,:])
        IR=np.nanmean(np.dstack(IR),axis=2)
        if returnCDs:
            IRcd=np.nanmean(np.dstack(IRcd),axis=2)

        IL=[]
        ILcd=[]
        for trial in ILm:
            cdtrials = np.setxor1d(trial,ILm)
            x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
            if np.any(fixedCDs):
                x = np.nanmean(x[:,fixedCDsF[0][0]:fixedCDsF[0][1]],axis=1)
                x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            cd = (RC[0]*code[0]+LC[0]*code[1]+RI[0]*code[2]+x*code[3])/2
            if np.any(ortho):
                ocds = []
                for i,oc in enumerate(ortho):      
                    if np.any(fixedCDs):
                        x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
                        x = np.nanmean(x[:,fixedCDsF[i+1][0]:fixedCDsF[i+1][1]],axis=1)
                        x = (x.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    ocds.append((RC[i+1]*oc[0]+LC[i+1]*oc[1]+RI[i+1]*oc[2]+x*oc[3])/2)
                acds = np.dstack([cd]+ocds)
                acds = acds[:,:, orthoOrder]
                cdIdx = (orthoOrder==0).nonzero()[0][0]
                if np.min(acds.shape)>2:
                    cd=[]
                    for tp in range(acds.shape[1]):
                        [ocds,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                        cd.append(ocds[:,cdIdx]*r[cdIdx,cdIdx])
                    cd = np.vstack(cd).T
            if returnCDs:
                ILcd.append(cd)
            IL.append(cd*psTsep[trial,:,:])
        IL=np.nanmean(np.dstack(IL),axis=2)
        if returnCDs:
            ILcd=np.nanmean(np.dstack(ILcd),axis=2)

    tsep3 = np.array([CR, CL, IR, IL])

    mT_preTType.append(tsep3)

    mT_preTType = np.dstack(mT_preTType)
    if returnCDs:
        cdAxes = np.dstack([np.array([CRcd, CLcd, IRcd, ILcd])])
        return mT_preTType, cdAxes
    return mT_preTType

def pullAnmTsepComputeSMOXvalHH(allAnmParams, anm, anat, phase, frs, alignTime, factors, minTrial, code, keepRoisDS,  keepRoisFull, \
    ortho=False, orthoOrder=False, fixedCDs=False, dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events',useDS = True, binary = False, returnCDs = False):
    """
    Split-half ("half-half") version of pullAnmTsepComputeSMOXval. Each condition's trials are split
    50:50 with uniqueSplit, and the full leave-one-out cross-validated SMO coding-direction (CD)
    projection is computed independently within each half. Returns two tensors (one per half), each
    shaped like pullAnmTsepComputeSMOXval's single return (4-conditions=[CR,CL,IR,IL], nRoi, nTime).
    The 50:50 split is drawn once per call, so downstream bootstrap variance comes only from animal
    resampling (matching how the non-split SMO is precomputed). See also [[gather_anmShuff_megaTsepHHSMO]].

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anm = 'JS041'                                       # animal id key into allAnmParams
    anat = 'somas'                                      # 'somas' | 'dendrites'
    phase = 'pre'                                        # behavioral phase key in shiftMask (e.g. 'pre' | 'post')
    frs = [30, 15]                                       # [dendrite_fr, soma_fr] frame rates (Hz)
    alignTime = 'goCues'                                 # 'goCues' | 'cTimes' | 'lastLick'
    factors = [4, 2]                                     # [dend, soma] temporal downsample factors for dsTsep
    minTrial = 5                                         # NaN both halves if min FULL condition count <= this (matches non-split SMO)
    code = [1, -1, 1, -1]                                # length-4 CD weights over [CR, CL, IR, IL]
    ortho = False                                        # False -> no orthogonalization; or list of length-4 weight vectors
    orthoOrder = False                                   #     > only if ortho set: np.array reorder for QR (index of primary CD == 0)
    fixedCDs = False                                     # False -> time-resolved CD; or list of [start_s, end_s] windows (rel. to align)
    dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events'    # key selecting the dff/events trace in allAnmParams
    useDS = True                                         # True -> dsTsep (downsampled); False -> fullTsep
    binary = False                                       # True | False  (binary currently unused in SMO math, kept for signature parity)
    returnCDs = False                                    # False -> return projections only; True -> also return per-half CD axes
    importlib.reload(js_manuscript)
    half1, half2 = js_manuscript.pullAnmTsepComputeSMOXvalHH(
        allAnmParams, anm=anm, anat=anat, phase=phase, frs=frs, alignTime=alignTime, factors=factors,
        minTrial=minTrial, code=code, keepRoisDS=keepRoisDS, keepRoisFull=keepRoisFull, ortho=ortho,
        orthoOrder=orthoOrder, fixedCDs=fixedCDs, dffKey=dffKey, useDS=useDS, binary=binary, returnCDs=returnCDs)
    # when returnCDs=True: half1, half2, cdAxes1, cdAxes2 = js_manuscript.pullAnmTsepComputeSMOXvalHH(..., returnCDs=returnCDs)

    Parameters
    ----------
    allAnmParams : dict
        Master per-animal parameter/data dict (tsep tensors, tracker, masks, etc.).
    anm : str
        Animal id key into allAnmParams.
    anat : str
        Anatomy key, 'somas' or 'dendrites'; selects the day mask and roi set.
    phase : str
        Behavioral phase key indexing allAnmParams[anm]['shiftMask'].
    frs : list of int
        [dendrite_fr, soma_fr] frame rates in Hz.
    alignTime : str
        Alignment event: 'goCues', 'cTimes', or 'lastLick'.
    factors : list of int
        [dend, soma] temporal downsample factors used to convert seconds to frames for fixedCDs.
    minTrial : int
        Minimum trials per condition, evaluated on the FULL (pre-split) counts to match
        pullAnmTsepComputeSMOXval; both halves are NaN-filled if the animal's smallest full
        condition count is <= this. This keeps the included-animal set identical to the non-split
        SMO (each surviving condition has >=6 trials, so every half has >=3).
    code : list or np.ndarray
        Length-4 CD weight vector over conditions [CR, CL, IR, IL].
    keepRoisDS : dict
        keepRoisDS[anat][anm] boolean roi mask for the downsampled tsep.
    keepRoisFull : dict
        keepRoisFull[anat][anm] boolean roi mask for the full-resolution tsep (used when useDS=False).
    ortho : bool or list
        False for no orthogonalization, or a list of length-4 weight vectors defining the CDs to
        orthogonalize against the primary CD (via per-timepoint QR).
    orthoOrder : bool or np.ndarray
        Only used when ortho is set: reorder array for the QR stack; the primary CD's slot must be 0.
    fixedCDs : bool or list
        False for a time-resolved CD, or a list of [start_s, end_s] windows (relative to alignTime)
        over which the CD template is averaged into a fixed vector.
    dffKey : str
        Key selecting the trace type inside allAnmParams[anm].
    useDS : bool
        True to use 'dsTsep' (downsampled), False to use 'fullTsep'.
    binary : bool
        Kept for signature parity with pullAnmTsepComputeSMOXval; not used by the SMO math.
    returnCDs : bool
        False returns only the two per-half projection tensors. True additionally
        returns `cdAxes1, cdAxes2`, the trial-averaged coding-direction axes
        (before projection) for each half, as third and fourth outputs.

    Returns
    -------
    half1, half2 : np.ndarray
        Per-half SMO CD projection tensors, each shaped (4, nRoi, nTime) over the
        four trial types [CR, CL, IR, IL]; all-NaN when the animal fails the
        minimum-trial guard.
    cdAxes1, cdAxes2 : np.ndarray
        Only returned when `returnCDs` is True. Per-half trial-averaged
        coding-direction axes (the `cd` vector, after any orthogonalization/
        fixed-window collapse), each shaped (4, nRoi, nTime) and keyed the same
        as half1/half2.
    """

    if anat == 'dendrites':
        ak = 0
    else:
        ak = 1

    fr = frs[ak]
    if useDS:
        tsepKey = 'dsTsep'
        factor = factors[ak]
        keepRois = keepRoisDS[anat][anm]
    else:
        tsepKey = 'fullTsep'
        factor = 1
        keepRois = keepRoisFull[anat][anm]

    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'lastLick':
        alignKey = 'll-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]

    if np.any(fixedCDs):
        nFrames=tsep.shape[2]
        hF = int(nFrames/2)
        fixedCDsF=[]
        for fcd in fixedCDs:
            fCDStart= hF+round(fcd[0]*fr/factor)
            fCDend= hF+round(fcd[1]*fr/factor)
            fixedCDsF.append([fCDStart,fCDend])
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays

    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]

    tkr = tracker[fullMask,:]

    psTsep = tsep2[:,keepRois,:]

    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    IRm = tkr[:,2]==5
    ILm = tkr[:,2]==6

    # 50:50 split of each condition's trial indices, drawn once (see pullAnmHH / uniqueSplit).
    CRts = uniqueSplit(CRm.nonzero()[0])
    CLts = uniqueSplit(CLm.nonzero()[0])
    IRts = uniqueSplit(IRm.nonzero()[0])
    ILts = uniqueSplit(ILm.nonzero()[0])

    nanShape = (psTsep.shape[1], psTsep.shape[2])  # (nRoi, nTime)

    def _cond_smo(slot, condIdx, RC, LC, RI, LI):
        # Leave-one-out cross-validated CD projection for one condition within one half.
        # `slot` selects which of [RC,LC,RI,LI] is replaced by the held-out-mean `x` (0=CR,1=CL,2=IR,3=IL).
        # When returnCDs is set, also returns the trial-averaged CD axis (the `cd` vector before projection).
        out = []
        cds = []
        for trial in condIdx:
            cdtrials = np.setxor1d(trial, condIdx)
            x = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
            if np.any(fixedCDs):
                xw = np.nanmean(x[:,fixedCDsF[0][0]:fixedCDsF[0][1]],axis=1)
                xw = (xw.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            else:
                xw = x
            comps = [RC[0], LC[0], RI[0], LI[0]]
            comps[slot] = xw
            cd = (comps[0]*code[0]+comps[1]*code[1]+comps[2]*code[2]+comps[3]*code[3])/2
            if np.any(ortho):
                ocds = []
                for i,oc in enumerate(ortho):
                    if np.any(fixedCDs):
                        xo = np.nanmean(psTsep[cdtrials,:,:], axis = 0)
                        xo = np.nanmean(xo[:,fixedCDsF[i+1][0]:fixedCDsF[i+1][1]],axis=1)
                        xo = (xo.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    else:
                        xo = x
                    ocomps = [RC[i+1], LC[i+1], RI[i+1], LI[i+1]]
                    ocomps[slot] = xo
                    ocds.append((ocomps[0]*oc[0]+ocomps[1]*oc[1]+ocomps[2]*oc[2]+ocomps[3]*oc[3])/2)
                acds = np.dstack([cd]+ocds)
                acds = acds[:,:, orthoOrder]
                cdIdx = (orthoOrder==0).nonzero()[0][0]
                if np.min(acds.shape)>2:
                    cd=[]
                    for tp in range(acds.shape[1]):
                        [ocds2,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                        cd.append(ocds2[:,cdIdx]*r[cdIdx,cdIdx])
                    cd = np.vstack(cd).T
            if returnCDs:
                cds.append(cd)
            out.append(cd*psTsep[trial,:,:])
        proj = np.nanmean(np.dstack(out), axis=2)
        if returnCDs:
            return proj, np.nanmean(np.dstack(cds), axis=2)
        return proj

    # Inclusion guard on FULL condition counts (matching pullAnmTsepComputeSMOXval), so the same
    # animals survive here as in the non-split SMO; each half then uses whatever trials its 50:50
    # split received. Guard requires min full count > minTrial (>=6), so every half has >=3 trials.
    nTrialsFull = [int(np.sum(CRm)), int(np.sum(CLm)), int(np.sum(IRm)), int(np.sum(ILm))]

    halves = []
    halvesCD = []
    for h in range(2):
        crH, clH, irH, ilH = CRts[h], CLts[h], IRts[h], ILts[h]

        if np.min(nTrialsFull) <= minTrial:
            CR = np.full(nanShape, np.nan)
            CL = np.full(nanShape, np.nan)
            IR = np.full(nanShape, np.nan)
            IL = np.full(nanShape, np.nan)
            if returnCDs:
                CRcd = np.full(nanShape, np.nan)
                CLcd = np.full(nanShape, np.nan)
                IRcd = np.full(nanShape, np.nan)
                ILcd = np.full(nanShape, np.nan)
        else:
            # Per-half CD templates (RC/LC/RI/LI), built only from this half's trials.
            RC=[]; LC=[]; RI=[]; LI=[]
            if np.any(fixedCDs):
                for fcd in fixedCDsF:
                    RCi = np.nanmean(psTsep[crH,:,:], axis = 0)
                    RCi = np.nanmean(RCi[:,fcd[0]:fcd[1]],axis=1)
                    RC.append(RCi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    LCi = np.nanmean(psTsep[clH,:,:], axis = 0)
                    LCi = np.nanmean(LCi[:,fcd[0]:fcd[1]],axis=1)
                    LC.append(LCi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    RIi = np.nanmean(psTsep[irH,:,:], axis = 0)
                    RIi = np.nanmean(RIi[:,fcd[0]:fcd[1]],axis=1)
                    RI.append(RIi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
                    LIi = np.nanmean(psTsep[ilH,:,:], axis = 0)
                    LIi = np.nanmean(LIi[:,fcd[0]:fcd[1]],axis=1)
                    LI.append(LIi.reshape(-1,1)@np.ones([1,psTsep.shape[2]]))
            else:
                for i in range(3):
                    RC.append(np.nanmean(psTsep[crH,:,:], axis = 0))
                    LC.append(np.nanmean(psTsep[clH,:,:], axis = 0))
                    RI.append(np.nanmean(psTsep[irH,:,:], axis = 0))
                    LI.append(np.nanmean(psTsep[ilH,:,:], axis = 0))

            if returnCDs:
                CR, CRcd = _cond_smo(0, crH, RC, LC, RI, LI)
                CL, CLcd = _cond_smo(1, clH, RC, LC, RI, LI)
                IR, IRcd = _cond_smo(2, irH, RC, LC, RI, LI)
                IL, ILcd = _cond_smo(3, ilH, RC, LC, RI, LI)
            else:
                CR = _cond_smo(0, crH, RC, LC, RI, LI)
                CL = _cond_smo(1, clH, RC, LC, RI, LI)
                IR = _cond_smo(2, irH, RC, LC, RI, LI)
                IL = _cond_smo(3, ilH, RC, LC, RI, LI)

        tsep3 = np.array([CR, CL, IR, IL])
        halves.append(np.dstack([tsep3]))
        if returnCDs:
            halvesCD.append(np.dstack([np.array([CRcd, CLcd, IRcd, ILcd])]))

    if returnCDs:
        return halves[0], halves[1], halvesCD[0], halvesCD[1]
    return halves[0], halves[1]

def gather_anmShuff_megaTsepFullSMO(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, cdKey, dffKey, alignTime, anmsToUse, nRois, anmIdxs, preComp, enforceMinTrial = False,
                               binary = False, returnErrors = False, anmsMinTrial=None):

    allMega = {phase: [] for phase in phases}
    allAnmTsep = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmTsep[anm] = {}
        for phase in phases:
            # anmTsep = pullAnmTsep(allAnmParams, anm, anat, ak, phase, frs, keepRoisDS, alignTime = alignTime,
            #          dffKey = dffKey, returnErrors = returnErrors, binary = binary, enforceMinTrial = enforceMinTrial,
            #                      anmsMinTrial=anmsMinTrial)
            # allAnmTsep[anm][phase] = anmTsep # trial, roi, time
            anmTsep = preComp[anat][anm][phase][alignTime][cdKey]
            allAnmTsep[anm][phase] = anmTsep # trial, roi, time
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmTsep[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega[phase].append(allAnmTsep[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega[phase] = np.dstack(allMega[phase]) #trial, time, roi
    return allMega

def gather_anmShuff_megaTsepHHSMO(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, cdKey, dffKey, alignTime, anmsToUse, nRois, anmIdxs, preCompHH, enforceMinTrial = False,
                               binary = False, returnErrors = False, anmsMinTrial=None):
    """
    Split-half ("half-half") counterpart of gather_anmShuff_megaTsepFullSMO. Reads the two
    precomputed 50:50 SMO tensors stored per animal in preCompHH[anat][anm][phase][alignTime][cdKey]
    (each entry a (half1, half2) tuple produced by [[pullAnmTsepComputeSMOXvalHH]]) and gathers them,
    roi-by-roi, into two parallel mega tensors — mirroring the two-half return of
    gather_anmShuff_megaTsepHH while keeping the SMO CD projection of FullSMO. Returns
    (allMega1, allMega2), each a {phase: (4-conditions, time, roi)} dict.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anat = 'somas'                                       # 'somas' | 'dendrites'
    ak = 1                                               # anat index: 0 dendrites, 1 somas
    mod = 'GO'                                           # alignment/error shorthand tag ('GO' | 'CT' | 'errGO' ...)
    phases = ['pre']                                     # list of phase keys to gather
    frs = [30, 15]                                       # [dendrite_fr, soma_fr] frame rates (Hz)
    cdKey = 's'                                           # SMO axis key: 's' | 'm' | 'o'
    dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events'    # trace key (unused here; parity with FullSMO)
    alignTime = 'goCues'                                 # 'goCues' | 'cTimes' | 'lastLick'
    nRois = 250                                          # total rois to gather (pooled across animals)
    anmIdxs = np.arange(len(anmsToUse['somas']))         # animal indices (bootstrap-resampled upstream)
    enforceMinTrial = False                              # True | False (unused here; parity with FullSMO)
    binary = False                                       # True | False (unused here; parity with FullSMO)
    returnErrors = False                                 # True | False (unused here; parity with FullSMO)
    anmsMinTrial = None                                  #     > only if enforceMinTrial: dict of per-anm min trial counts
    importlib.reload(js_manuscript)
    allMega1, allMega2 = js_manuscript.gather_anmShuff_megaTsepHHSMO(
        allAnmParams, anat=anat, ak=ak, mod=mod, phases=phases, frs=frs, keepRoisDS=keepRoisDS,
        cdKey=cdKey, dffKey=dffKey, alignTime=alignTime, anmsToUse=anmsToUse, nRois=nRois,
        anmIdxs=anmIdxs, preCompHH=preCompHH, enforceMinTrial=enforceMinTrial, binary=binary,
        returnErrors=returnErrors, anmsMinTrial=anmsMinTrial)

    Parameters
    ----------
    allAnmParams : dict
        Master per-animal parameter/data dict (accepted for signature parity; unused here).
    anat : str
        Anatomy key, 'somas' or 'dendrites'.
    ak : int
        Anatomy index (0 dendrites, 1 somas); accepted for parity, unused here.
    mod : str
        Alignment/error shorthand tag; accepted for parity, unused here.
    phases : list of str
        Phase keys to gather (e.g. ['pre']).
    frs : list of int
        [dendrite_fr, soma_fr] frame rates; accepted for parity, unused here.
    keepRoisDS : dict
        keepRoisDS[anat][anm] roi masks; accepted for parity, unused here.
    cdKey : str
        SMO axis key ('s', 'm', or 'o') selecting which precomputed CD projection to gather.
    dffKey : str
        Trace key; accepted for parity, unused here.
    alignTime : str
        Alignment event key indexing preCompHH ('goCues' | 'cTimes' | 'lastLick').
    anmsToUse : dict
        anmsToUse[anat] -> list of animal ids; indexed by anmIdxs.
    nRois : int
        Total number of rois to pool across animals into the mega tensor.
    anmIdxs : np.ndarray
        Animal indices into anmsToUse[anat] (bootstrap-resampled by the caller).
    preCompHH : dict
        preCompHH[anat][anm][phase][alignTime][cdKey] -> (half1, half2) SMO tensors, each
        (4-conditions, nRoi, nTime), from pullAnmTsepComputeSMOXvalHH.
    enforceMinTrial : bool
        Accepted for parity with FullSMO; unused here (split/minTrial handled in preCompHH build).
    binary : bool
        Accepted for parity; unused here.
    returnErrors : bool
        Accepted for parity; unused here.
    anmsMinTrial : dict or None
        Accepted for parity; unused here.
    """

    allMega1 = {phase: [] for phase in phases}
    allMega2 = {phase: [] for phase in phases}
    allAnmTsep1 = {}
    allAnmTsep2 = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmTsep1[anm] = {}
        allAnmTsep2[anm] = {}
        for phase in phases:
            half1, half2 = preCompHH[anat][anm][phase][alignTime][cdKey]
            allAnmTsep1[anm][phase] = half1 # trial, roi, time
            allAnmTsep2[anm][phase] = half2

    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]

            nRoi = allAnmTsep1[anm][phases[0]].shape[1]

            for roi in range(nRoi):
                for phase in phases:
                    allMega1[phase].append(allAnmTsep1[anm][phase][:,roi,:])
                    allMega2[phase].append(allAnmTsep2[anm][phase][:,roi,:])

                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break

    for phase in phases:
        allMega1[phase] = np.dstack(allMega1[phase]) #trial, time, roi
        allMega2[phase] = np.dstack(allMega2[phase])
    return allMega1, allMega2

def get_SMO_fast(megaTsepsS, megaTsepsM, megaTsepsO, nBoot, anat, mod, calcPhase='pre1', center = False):
    """
    Collapse the bootstrapped SMO coding-direction (CD) tensors to population-mean
    projection traces per trial type, stacked over the three CDs.

    For each bootstrap iteration the `calcPhase` slice of the sensory, motor, and
    outcome mega-tensors (each shape (4, nTime, nRoi)) is taken for the given
    `anat` and `mod`, NaNs are zeroed, and the ROI axis is averaged out to give a
    (4, nTime) population mean per trial type [CR, CL, IR, IL]. The three CD means
    are stacked along a new last axis and accumulated across bootstraps.

    Example
    -------
    megaTsepsS = megaTsepsS     # sensory CD mega-tensor from full_get_SMO, keyed [phase][mod][anat] -> (4, nTime, nRoi, nBoot)
    megaTsepsM = megaTsepsM     # motor CD mega-tensor, keyed the same
    megaTsepsO = megaTsepsO     # outcome CD mega-tensor, keyed the same
    nBoot = 1000                # number of bootstrap iterations to read (final axis of each tensor)
    anat = 'dendrites'          # anatomy key into the mega-tensors, e.g. 'dendrites' | 'somas'
    mod = 'GO'                  # alignment/mod key into the mega-tensors, e.g. 'GO' | 'CT'
    calcPhase = 'pre'           # phase key into the mega-tensors and the output 'projs' dict, e.g. 'pre' | 'post'
    center = False              # if True, intended to mean-center each CD over time first (see Notes)
    importlib.reload(js_manuscript)
    allData = js_manuscript.get_SMO_fast(
        megaTsepsS, megaTsepsM=megaTsepsM, megaTsepsO=megaTsepsO, nBoot=nBoot,
        anat=anat, mod=mod, calcPhase=calcPhase, center=center)

    Parameters
    ----------
    megaTsepsS : dict
        Sensory-CD bootstrapped tensor keyed `[phase][mod][anat]` with array shape
        (4, nTime, nRoi, nBoot); the `calcPhase` phase is read.
    megaTsepsM : dict
        Motor-CD bootstrapped tensor, keyed and shaped like `megaTsepsS`.
    megaTsepsO : dict
        Outcome-CD bootstrapped tensor, keyed and shaped like `megaTsepsS`.
    nBoot : int
        Number of bootstrap iterations to iterate over (indexes the final axis).
    anat : str
        Anatomy key into the mega-tensors (e.g. 'dendrites', 'somas').
    mod : str
        Alignment/mod key into the mega-tensors (e.g. 'GO', 'CT').
    calcPhase : str
        Phase key selecting which slice of the mega-tensors to read (e.g. 'pre',
        'post'); also used as the key under `allData['projs']` in the output.
    center : bool
        If True, intended to mean-center each CD tensor over time before averaging.
        See Notes: as written this has no effect.

    Returns
    -------
    allData : dict
        `{'projs': {calcPhase: ndarray}}` where the array has shape (nBoot, 4, nTime, 3):
        bootstrap x trial type [CR, CL, IR, IL] x time x CD [sensory, motor, outcome].

    Notes
    -----
    The `center` branch rebinds a local loop variable rather than the array, so it
    currently does not modify the data; the NaN-zeroing in the following loop does
    mutate the tensors in place.
    """
    allData = {'projs': {calcPhase: []}}

    for i in range(nBoot):
        s  = megaTsepsS[calcPhase][mod][anat][:,:,:,i]
        m  = megaTsepsM[calcPhase][mod][anat][:,:,:,i]
        o  = megaTsepsO[calcPhase][mod][anat][:,:,:,i]

        nTime = s.shape[1]

        allD = [s, m, o]
        if center:
            # BUG (not fixed): `d = d - ...` rebinds the local loop variable instead of
            # mutating the array, so centering has no effect on s/m/o. Left as-is for now.
            for d in allD:
                d = d-np.nanmean(d, axis = 1, keepdims=True)

        for d in allD:
            d[np.isnan(d)] = 0

        data = [s,m,o]
        projs = []
        for a in range(3):
            axProjs = {
                calcPhase: np.nanmean(data[a],axis=2)
            }
            projs.append(axProjs)
        
        for key in list(axProjs.keys()):
            allData['projs'][key].append(np.dstack([projs[a][key] for a in range(3)]))

    for key in list(allData['projs'].keys()):
        allData['projs'][key] = np.stack(allData['projs'][key])
    
    return allData

def get_SMO_CD2_fast(megaTsepsHH, anat, mod, cdType = 's', orth=False,calcPhase='post1', square = False, cdTrials = [0,3], nBoot = 1000, normalize = True):
    """
    Build bootstrapped SMO coding-direction (CD) axis vectors from the half-split
    (HH) mega-tensor for a chosen split phase.

    For each bootstrap iteration the four half-split slices (pre1, pre2, post1,
    post2), each shape (4, nTime, nRoi), are read for the given `anat` and `mod`
    and NaN-zeroed. The slice named by `calcPhase` is selected, its four trial-type
    groups (indexed by `cdTrials`, order [CR, CL, IR, IL]) are combined into the
    sensory/motor/outcome contrast chosen by `cdType`. When `orth` is set, the
    chosen contrast is orthogonalized against the other two (QR, rescaled by the
    diagonal R term); otherwise the raw contrast is used. The resulting per-time
    ROI vector is stored as `axesOrig`, and — when `normalize` is True — an
    L2-normalized-across-ROIs (unit) copy is stored as `axes`. Both are accumulated
    across bootstraps.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    megaTsepsHH = megaTsepsHH   # combined half-split tensor from full_get_SMO, keyed [pre1|pre2|post1|post2][f'err{mod}'][anat] -> (4, nTime, nRoi, nBoot)
    anat = 'dendrites'          # anatomy key into the tensor, e.g. 'dendrites' | 'somas'
    mod = 'GO'                  # alignment/mod key into the tensor, e.g. 'GO' | 'CT' (read as f'err{mod}')
    cdType = 's'                # 's' sensory | 'm' motor | anything else -> outcome contrast
    orth = False                # True -> QR-orthogonalize chosen contrast vs the other two; False -> raw contrast
    calcPhase = 'post1'         # which half-split slice to build the axis from: 'pre1' | 'pre2' | 'post1' | 'post2'
    square = False              # accepted for signature compatibility; not used in this function
    cdTrials = [0, 1, 2, 3]     # length-4 trial-type indices [CR, CL, IR, IL] selecting the contrast groups
    nBoot = 1000                # number of bootstrap iterations to read (final axis of the tensor)
    normalize = True            # True -> 'axes' is per-timepoint L2-normalized; False -> 'axes' is a raw copy of 'axesOrig'
    importlib.reload(js_manuscript)
    allData = js_manuscript.get_SMO_CD2_fast(
        megaTsepsHH, anat=anat, mod=mod, cdType=cdType, orth=orth, calcPhase=calcPhase,
        square=square, cdTrials=cdTrials, nBoot=nBoot, normalize=normalize)

    Parameters
    ----------
    megaTsepsHH : dict
        Combined half-split bootstrapped tensor keyed
        `['pre1' | 'pre2' | 'post1' | 'post2'][f'err{mod}'][anat]` with array shape
        (4, nTime, nRoi, nBoot); all four phase slices are read.
    anat : str
        Anatomy key into the tensor (e.g. 'dendrites', 'somas').
    mod : str
        Alignment/mod key into the tensor (e.g. 'GO', 'CT'); indexed as `f'err{mod}'`.
    cdType : str
        CD contrast to build: 's' sensory ((CR+IR)-(CL+IL))/2, 'm' motor
        ((CR+IL)-(CL+IR))/2, otherwise outcome ((CR+CL)-(IR+IL))/2.
    orth : bool
        If True, orthogonalize the chosen contrast against the other two contrasts
        per timepoint via QR (rescaled by the diagonal R element). If False, use the
        raw contrast vector directly.
    calcPhase : str
        Name of the half-split slice to build the axis from, one of 'pre1',
        'pre2', 'post1', 'post2' (resolved by `eval`).
    square : bool
        Accepted for call-signature compatibility; not referenced in this function.
    cdTrials : list of int
        Length-4 trial-type indices in [CR, CL, IR, IL] order used to pick the four
        contrast groups from the selected slice. (The default `[0, 3]` is a stale
        2-element value; callers pass a full 4-element list.)
    nBoot : int
        Number of bootstrap iterations to iterate over (indexes the final axis).
    normalize : bool
        If True (default), `axes` is `axesOrig` L2-normalized across ROIs at each
        timepoint (unit vectors). If False, `axes` is an unnormalized copy of
        `axesOrig`.

    Returns
    -------
    allData : dict
        `{'axesOrig': ndarray, 'axes': ndarray}`, each of shape (nBoot, nTime, nRoi):
        the raw contrast axis and, in `axes`, its per-timepoint L2-normalized (unit)
        version when `normalize` is True, or an unnormalized copy when False.
    """
    allData = {'axes': [], 'axesOrig': []}
    
    for i in tqdm(range(nBoot)):
        post1 = megaTsepsHH["post1"][f'{mod}'][anat][:,:,:,i].copy()
        post2 = megaTsepsHH["post2"][f'{mod}'][anat][:,:,:,i].copy()
        pre1  = megaTsepsHH["pre1"][f'{mod}'][anat][:,:,:,i].copy()
        pre2  = megaTsepsHH["pre2"][f'{mod}'][anat][:,:,:,i].copy()

        allD = [pre1,pre2,post1,post2]
        for d in allD:
            d[np.isnan(d)] = 0
            
        res = eval(calcPhase).copy()
        group1 = res[cdTrials[0],:,:] # CR
        group2 = res[cdTrials[1],:,:] # CL
        group3 = res[cdTrials[2],:,:] # IR
        group4 = res[cdTrials[3],:,:] # IL

        s = ((group1 + group3) - (group2 + group4))/2
        m = ((group1 + group4) - (group2 + group3))/2
        o = ((group1 + group2) - (group3 + group4))/2
        orthOrder=np.array([2,1,0])
        if cdType == 's':
            cdIdx = 0
        elif cdType == 'm':
            cdIdx = 1
        else:
            cdIdx = 2
        cds = [s,m,o]
        acds = copy.deepcopy(cds)
        acds.pop(cdIdx)
        acds = [cds[cdIdx]]+acds
        acds = np.dstack(acds)
        acds = acds[:,:,orthOrder]
        cdOut = (orthOrder==0).nonzero()[0][0]
        if orth:
            cd=[]
            for tp in range(acds.shape[1]):
                [ocds,r] = np.linalg.qr(acds[:,tp,:],mode='complete')
                cd.append(ocds[:,cdOut]*r[cdOut,cdOut]) 
            axesOrig = np.vstack(cd).T
        else:
            axesOrig = acds[:,:,cdOut]
        if normalize:
            norms = np.linalg.norm(axesOrig, axis=1, keepdims=True)
            axes = axesOrig / norms         # normalized per timepoint
        else:
            axes = axesOrig.copy()          # normalization disabled: pass through raw axis

        allData['axesOrig'].append(axesOrig)
        allData['axes'].append(axes)       
    allData['axesOrig'] = np.stack(allData['axesOrig'])
    allData['axes'] = np.stack(allData['axes'])
    
    return allData

def pullAnmHH_ALT_SMO(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, keepRoisFull, alignTime = 'goCues',
                 dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events', factors = [4,2],
              useDS = True, enforceMinTrial = False, nTriToDraw = 20, binary = False, giveAll = False):
    
    mT_preTType1 = []
    mT_preTType2 = []
    mTs = [mT_preTType1, mT_preTType2]
    fr = frs[ak]
    if useDS:
        tsepKey = 'dsTsep'
        factor = factors[ak]
        keepRois = keepRoisDS[anat][anm]
    else:
        tsepKey = 'fullTsep'
        factor = 1
        keepRois = keepRoisFull[anat][anm]
        
    if alignTime == 'goCues':
        alignKey = 'go-psTsep'
    elif alignTime == 'cTimes':
        alignKey = 'ct-psTsep'
    elif alignTime == 'secondLick':
        alignKey = 'sl-psTsep'
    else:
        print('incorrect alignKey passed')

    tsep = allAnmParams[anm][dffKey][tsepKey][anat][alignKey]
        
    tracker = allAnmParams[anm]['tracker']
    dendDays = allAnmParams[anm]['dendDays']
    cTimes = allAnmParams[anm]['cTimes']
    goCues = allAnmParams[anm]['goCues']
    phaseMask = allAnmParams[anm]['shiftMask'][phase]
    de = allAnmParams[anm]['de']
    me = allAnmParams[anm]['me']
    ar = allAnmParams[anm]['afterRight']
    al_noME = allAnmParams[anm]['MER-L_noME']
    al_ME = allAnmParams[anm]['MER-L_ME']
    an = allAnmParams[anm]['afterNone']
    slN = allAnmParams[anm]['sl_nHits']
    slTracking = allAnmParams[anm]['slSideTracking']
    
    if anat == 'dendrites':
        anatMask = dendDays
    else:
        anatMask = ~dendDays

    fullMask = np.logical_and(phaseMask, anatMask)
    tsep2 = tsep[phaseMask[anatMask],:,:]  
    tkr = tracker[fullMask,:]
    
    psTsep = tsep2[:,keepRois,:]
    
    de2 = de[fullMask]
    me2 = me[fullMask]

    CRm = tkr[:,2]==0
    CLm = tkr[:,2]==1
    IRm = tkr[:,2]==5
    ILm = tkr[:,2]==6
    
    FRm = tkr[:,2]==5
    MEm = np.logical_and(FRm,me2)
    DEm = np.logical_and(FRm,de2)
    ARm = ar[fullMask]
    AL_noMEm = al_noME[fullMask]
    AL_MEm = al_ME[fullMask]
    ANm = an[fullMask]
    noSecondHitMask = np.logical_and(slN[fullMask]==0, MEm)
    hitMask = np.logical_and(slN[fullMask]>0, MEm)
    slAttmptR = np.logical_and(MEm, slTracking[fullMask]==0)
    slAttmptL = np.logical_and(MEm, slTracking[fullMask]==1)

    hasIR = np.sum(IRm)>2
    hasIL = np.sum(ILm)>2
    hasAR = np.sum(ARm)>2
    hasAL_noME = np.sum(AL_noMEm)>2
    hasAL_ME = np.sum(AL_MEm)>2
    hasAN = np.sum(ANm)>2
    hasDE = np.sum(DEm)>2
    hasME = np.sum(MEm)>2
    hasNSH = np.sum(noSecondHitMask)>2
    hasHit = np.sum(hitMask)>2
    
    CRts = uniqueSplit(np.where(CRm)[0])
    CLts = uniqueSplit(np.where(CLm)[0])
    SLRts = uniqueSplit(np.where(slAttmptR)[0])
    SLLts = uniqueSplit(np.where(slAttmptL)[0])

    if hasIR:
        IRts = uniqueSplit(np.where(IRm)[0])
    if hasIL:
        ILts = uniqueSplit(np.where(ILm)[0])
        
    if hasDE:
        DEts = uniqueSplit(np.where(DEm)[0])
    if hasME:
        MEts = uniqueSplit(np.where(MEm)[0])
    if hasAR:
        ARts = uniqueSplit(np.where(ARm)[0])
    if hasAL_noME:
        AL_noMEts = uniqueSplit(np.where(AL_noMEm)[0])
    if hasAL_ME:
        AL_MEts = uniqueSplit(np.where(AL_MEm)[0])
    if hasAN:
        ANts = uniqueSplit(np.where(ANm)[0])
    if hasNSH:
        NSHts = uniqueSplit(np.where(noSecondHitMask)[0])
    if hasHit:
        Hit_ts = uniqueSplit(np.where(hitMask)[0])

    for i in range(2):
        if hasIR:
            irt = IRts[i]
        if hasIL:
            ilt = ILts[i]
        if hasDE:
            det = DEts[i]
        if hasME:
            met = MEts[i]
        if hasAR:
            art = ARts[i]
        if hasAL_noME:
            al_noMEt = AL_noMEts[i]
        if hasAL_ME:
            al_MEt = AL_MEts[i]
        if hasAN:
            ant = ANts[i]
        if hasNSH:
            nsh_t = NSHts[i]
        if hasHit:
            hit_t = Hit_ts[i]
            
        crt = CRts[i]
        clt = CLts[i]
        slrt = SLRts[i]
        sllt = SLLts[i]
        
        if enforceMinTrial:
            CRt = np.random.choice(crt, nTriToDraw)
            CLt = np.random.choice(clt, nTriToDraw)
            SLRt = np.random.choice(slrt, nTriToDraw)
            SLLt = np.random.choice(sllt, nTriToDraw)
            if hasIR:
                IRt = np.random.choice(irt, nTriToDraw)
            if hasIL:
                ILt = np.random.choice(ilt, nTriToDraw)
            if hasDE:
                DEt = np.random.choice(det, nTriToDraw)
            if hasME:
                MEt = np.random.choice(met, nTriToDraw)
            if hasAR:
                ARt = np.random.choice(art, nTriToDraw)
            if hasAL_ME:
                AL_MEt = np.random.choice(al_MEt, nTriToDraw)
            if hasAL_noME:
                AL_noMEt = np.random.choice(al_noMEt, nTriToDraw)
            if hasAN:
                ANt = np.random.choice(ant, nTriToDraw)
            if hasNSH:
                NSHt = np.random.choice(nsh_t, nTriToDraw)
            if hasHit:
                hitT = np.random.choice(hit_t, nTriToDraw)
 
        else:
            CRt = crt
            CLt = clt
            SLRt = slrt
            SLLt = sllt
            if hasIR:
                IRt = irt
            if hasIL:
                ILt = ilt
            if hasDE:
                DEt = det
            if hasME:
                MEt = met
            if hasAR:
                ARt = art
            if hasAL_ME:
                AL_MEt = al_MEt
            if hasAL_noME:
                AL_noMEt = al_noMEt
            if hasAN:
                ANt = ant
            if hasNSH:
                NSHt = nsh_t
            if hasHit:
                hitT = hit_t

        if binary:
            CR = np.nanmean(psTsep[CRt,:,:]>0, axis = 0)
        else:
            CR = np.nanmean(psTsep[CRt,:,:], axis = 0)

        if binary:
            SLR = np.nanmean(psTsep[SLRt,:,:]>0, axis = 0)
        else:
            SLR = np.nanmean(psTsep[SLRt,:,:], axis = 0)
            
        if binary:
            SLL = np.nanmean(psTsep[SLLt,:,:]>0, axis = 0)
        else:
            SLL = np.nanmean(psTsep[SLLt,:,:], axis = 0)

        if binary:
            CL = np.nanmean(psTsep[CLt,:,:]>0, axis = 0)
        else:
            CL = np.nanmean(psTsep[CLt,:,:], axis = 0)

        if hasIR:
            if binary:
                IR = np.nanmean(psTsep[IRt,:,:]>0, axis = 0)
            else:
                IR = np.nanmean(psTsep[IRt,:,:], axis = 0)
        else:
            IR = np.ones_like(CR)*np.nan

        if hasIL:
            if binary:
                IL = np.nanmean(psTsep[ILt,:,:]>0, axis = 0)
            else:
                IL = np.nanmean(psTsep[ILt,:,:], axis = 0)
        else:
            IL = np.ones_like(CR)*np.nan
        
        if hasNSH:
            if binary:
                NSH = np.nanmean(psTsep[NSHt,:,:]>0, axis = 0)
            else:
                NSH = np.nanmean(psTsep[NSHt,:,:], axis = 0)
        else:
            NSH = np.ones_like(CR)*np.nan

        if hasHit:
            if binary:
                HIT = np.nanmean(psTsep[hitT,:,:]>0, axis = 0)
            else:
                HIT = np.nanmean(psTsep[hitT,:,:], axis = 0)
        else:
            HIT = np.ones_like(CR)*np.nan
        
        if hasDE:
            if binary:
                DE = np.nanmean(psTsep[DEt,:,:]>0, axis = 0)
            else:
                DE = np.nanmean(psTsep[DEt,:,:], axis = 0)
        else:
            DE = np.ones_like(CR)*np.nan
        
        if hasME:
            if binary:
                ME = np.nanmean(psTsep[MEt,:,:]>0, axis = 0)
            else:
                ME = np.nanmean(psTsep[MEt,:,:], axis = 0)
        else:
            ME = np.ones_like(CR)*np.nan
            
        if hasAR:
            if binary:
                AR = np.nanmean(psTsep[ARt,:,:]>0, axis = 0)
            else:
                AR = np.nanmean(psTsep[ARt,:,:], axis = 0)
        else:
            AR = np.ones_like(CR)*np.nan

        if hasAL_ME:
            if binary:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:]>0, axis = 0)
            else:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:], axis = 0)
        else:
            AL_ME = np.ones_like(CR)*np.nan

        if hasAL_noME:
            if binary:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:]>0, axis = 0)
            else:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:], axis = 0)
        else:
            AL_noME = np.ones_like(CR)*np.nan

        if hasAN:
            if binary:
                AN = np.nanmean(psTsep[ANt,:,:]>0, axis = 0)
            else:
                AN = np.nanmean(psTsep[ANt,:,:], axis = 0)
        else:
            AN = np.ones_like(CR)*np.nan

        if giveAll:
            #tsep3 = np.array([CR,AL_noME,AL_ME,AR,AN,DE,CL,NSH,HIT])
            tsep3 = np.array([CR,AL_noME,AL_ME,AR,ME, AN,DE,CL,IR,IL])
        else:
            #tsep3 = np.array([CR, ME, AR, AL])
            tsep3 = np.array([CR, AL_noME, AL_ME, AR])
        
        mTs[i].append(tsep3)
                
    mT_preTType1 = np.dstack(mT_preTType1)
    mT_preTType2 = np.dstack(mT_preTType2)
    return mT_preTType1, mT_preTType2

def gather_anmShuff_megaTsepHH_ALT_SMO(allAnmParams, anat, ak, mod, phases, frs, \
    keepRoisDS, keepRoisFull, dffKey, alignTime, anmsToUse, nRois, anmIdxs, \
    enforceMinTrial = False, nTriToDraw = 50,
    binary = False, useDS=True, giveAll = False):

    allMega1 = {phase: [] for phase in phases}
    allMega2 = {phase: [] for phase in phases}
    LUT = {phase: [] for phase in phases}
    
    allAnmHH1 = {}
    allAnmHH2 = {}
    for anmIdx in anmIdxs:
        anm = anmsToUse[anat][anmIdx]
        allAnmHH1[anm] = {}
        allAnmHH2[anm] = {}
        for phase in phases:
            anm1, anm2 = pullAnmHH_ALT_SMO(allAnmParams, anm, anat, ak, phase, frs,  keepRoisDS, keepRoisFull, alignTime = alignTime,
                     dffKey = dffKey, useDS = useDS, enforceMinTrial = enforceMinTrial, 
                                               nTriToDraw = nTriToDraw, binary = binary, giveAll = giveAll)
            allAnmHH1[anm][phase] = anm1 # trial, roi, time
            allAnmHH2[anm][phase] = anm2
    
    roiCounter = 0
    while roiCounter < nRois:
        for anmIdx in anmIdxs:
            anm = anmsToUse[anat][anmIdx]
        
            nRoi = allAnmHH1[anm][phases[0]].shape[1]
    
            for roi in range(nRoi):
                for phase in phases:
                    allMega1[phase].append(allAnmHH1[anm][phase][:,roi,:])
                    allMega2[phase].append(allAnmHH2[anm][phase][:,roi,:])
                
                roiCounter += 1
                if roiCounter >= nRois:
                    break
            if roiCounter >= nRois:
                break
                
    for phase in phases:
        allMega1[phase] = np.dstack(allMega1[phase]) #trial, time, roi
        allMega2[phase] = np.dstack(allMega2[phase])
    return allMega1, allMega2

def get_CRMER_CD2_fast(megaTsepsHH, anat, mod, calcPhase='post1', square = False, cdTrials = [0,3], nBoot = 1000):
    def proj(mat, axes, sqr = True):
        prjn = np.einsum('gtr,tr->gt', mat, axes)
        if sqr:
            prjn = np.sqrt(abs(prjn))*np.sign(prjn)
        return prjn
    allData = {'projs': {'pre1': [], 'pre2': [], 'post1': [], 'post2': []}, 
               'axes': [], 'axesOrig': []}
    
    for i in range(nBoot):
        post1 = megaTsepsHH["post1"][f'err{mod}'][anat][:,:,:,i].copy()
        post2 = megaTsepsHH["post2"][f'err{mod}'][anat][:,:,:,i].copy()
        pre1  = megaTsepsHH["pre1"][f'err{mod}'][anat][:,:,:,i].copy()
        pre2  = megaTsepsHH["pre2"][f'err{mod}'][anat][:,:,:,i].copy()

        allD = [pre1,pre2,post1,post2]#, late1BS, late2BS]
        for d in allD:
            d[np.isnan(d)] = 0
            
        res = eval(calcPhase).copy()
        group1 = res[cdTrials[0],:,:]  # time, roi
        # group2 = res[1,:,:] # ME-L_noME
        #group2 = res[2,:,:]  # ME-L_ME
        group2 = res[cdTrials[1],:,:] # ME-R
        
        axesOrig = (group1 - group2)  # time, roi
        norms = np.linalg.norm(axesOrig, axis=1, keepdims=True)
        axes = axesOrig / norms         # normalized per timepoint
        
        allData['axesOrig'].append(axesOrig)
        allData['axes'].append(axes)    
        
        # projsOrig = {
        #     'pre1': proj(pre1, axesOrig, sqr = square),
        #     'pre2': proj(pre2, axesOrig, sqr = square),
        #     'post1': proj(post1, axesOrig, sqr = square),
        #     'post2': proj(post2, axesOrig, sqr = square)
        # }
        
        projs = {
            'pre1': proj(pre1, axes, sqr = square),
            'pre2': proj(pre2, axes, sqr = square),
            'post1': proj(post1, axes, sqr = square),
            'post2': proj(post2, axes, sqr = square)
        }
        
        for key, prj in projs.items():
            allData['projs'][key].append(prj)
        # for key, prj in projsOrig.items():
        #     allData['projsOrig'][key].append(prj)

    for key in ['pre1', 'pre2', 'post1', 'post2']:
        allData['projs'][key] = np.stack(allData['projs'][key])
    # for key in ['pre1', 'pre2', 'post1', 'post2']:
    #     allData['projsOrig'][key] = np.stack(allData['projsOrig'][key])
    
    allData['axesOrig'] = np.stack(allData['axesOrig'])
    allData['axes'] = np.stack(allData['axes'])
    
    return allData

def get_CRCL_CD_fast(megaTsepsHH, anat, mod, calcPhase='post1', nBoot = 1000):
    allData = {'projs': {'pre1': [], 'pre2': [], 'post1': [], 'post2': []}, 
               'axes': [], 'axesOrig': []}
    
    for i in range(nBoot):
        post1 = megaTsepsHH["post1"][mod][anat][:,:,:,i].copy()
        post2 = megaTsepsHH["post2"][mod][anat][:,:,:,i].copy()
        pre1  = megaTsepsHH["pre1"][mod][anat][:,:,:,i].copy()
        pre2  = megaTsepsHH["pre2"][mod][anat][:,:,:,i].copy()
        # post1 = megaTsepsHH["early1"][mod][anat][:,:,:,i].copy()
        # post2 = megaTsepsHH["early2"][mod][anat][:,:,:,i].copy()

        allD = [pre1, pre2, post1, post2]
        for d in allD:
            d[np.isnan(d)] = 0
        
        res = eval(calcPhase).copy()
        group1 = res[0,:,:]  # time, roi
        group2 = res[1,:,:]
        
        axesOrig = (group1 - group2)  # time, roi
        norms = np.linalg.norm(axesOrig, axis=1, keepdims=True)
        axes = axesOrig / norms         # normalized per timepoint
        
        # Store these
        allData['axesOrig'].append(axesOrig)
        allData['axes'].append(axes)    
        
        projs = {
            'pre1': np.einsum('gtr,tr->gt', pre1, axes),
            'pre2': np.einsum('gtr,tr->gt', pre2, axes),
            'post1': np.einsum('gtr,tr->gt', post1, axes),
            'post2': np.einsum('gtr,tr->gt', post2, axes)
        }
        
        # store all projections
        for key, prj in projs.items():
            allData['projs'][key].append(prj)
        # allData['projs']['pre2'].append(pre2P)
        # allData['projs']['post1'].append(post1P)
        # allData['projs']['post2'].append(post2P)

    # ---- stack bootstraps ----
    for key in ['pre1', 'pre2', 'post1', 'post2']:
        allData['projs'][key] = np.stack(allData['projs'][key])
    
    allData['axesOrig'] = np.stack(allData['axesOrig'])
    allData['axes'] = np.stack(allData['axes'])
    
    return allData

def get_CRCL_CD_fast_earlyLate(megaTsepsHH, anat, mod, calcPhase='pre1', sqrt = False, nBoot = 1000):
    allData = {'projs': {'pre1': [], 'pre2': [], 'early1': [], 'early2': [], 'late1': [], 'late2': []}, 
               'axes': [], 'axesOrig': []}
    
    for i in range(nBoot):
        early1 = megaTsepsHH["early1"][mod][anat][:,:,:,i].copy()
        early2 = megaTsepsHH["early2"][mod][anat][:,:,:,i].copy()
        late1 = megaTsepsHH["late1"][mod][anat][:,:,:,i].copy()
        late2 = megaTsepsHH["late2"][mod][anat][:,:,:,i].copy()
        # pre1 = megaTsepsHH["pre1"][mod][anat][:,:,:,i].copy()
        # pre2 = megaTsepsHH["pre2"][mod][anat][:,:,:,i].copy()
        pre1 = megaTsepsHH["preShift_late1"][mod][anat][:,:,:,i].copy()
        pre2 = megaTsepsHH["preShift_late2"][mod][anat][:,:,:,i].copy()

        allD = [pre1, pre2, early1, early2, late1, late2]
        for d in allD:
            d[np.isnan(d)] = 0

        res = eval(calcPhase).copy()
        group1 = res[0,:,:]  # time, roi
        group2 = res[1,:,:]
        
        axesOrig = (group1 - group2)  # time, roi
        norms = np.linalg.norm(axesOrig, axis=1, keepdims=True)
        axes = axesOrig / norms         # normalized per timepoint
        
        allData['axesOrig'].append(axesOrig)
        allData['axes'].append(axes)    
        
        projs = {
            'pre1': np.einsum('gtr,tr->gt', pre1, axes),
            'pre2': np.einsum('gtr,tr->gt', pre2, axes),
            'early1': np.einsum('gtr,tr->gt', early1, axes),
            'early2': np.einsum('gtr,tr->gt', early2, axes),
            'late1': np.einsum('gtr,tr->gt', late1, axes),
            'late2': np.einsum('gtr,tr->gt', late2, axes)
        }
        
        for key, prj in projs.items():
            if sqrt:
                allData['projs'][key].append(np.sqrt(abs(prj))*np.sign(prj))
            else:
                allData['projs'][key].append(prj)
        # allData['projs']['pre2'].append(pre2P)
        # allData['projs']['post1'].append(post1P)
        # allData['projs']['post2'].append(post2P)

    for key in list(projs.keys()):
        allData['projs'][key] = np.stack(allData['projs'][key])
    
    allData['axesOrig'] = np.stack(allData['axesOrig'])
    allData['axes'] = np.stack(allData['axes'])
    
    return allData

def get_mag_corr(CD1, CD2):
    norm1 = np.linalg.norm(CD1, axis=2)
    norm2 = np.linalg.norm(CD2, axis=2)
    dot = np.sum(CD1 * CD2, axis=2)
    cos_sim = dot / (norm1 * norm2)
    return norm1, cos_sim
    
def get_CRCL_CD_full(megaTseps, anat, mod, calcPhase='post', nBoot = 1000):
    allData = {'projs': {'pre': [], 'early': [], 'late': []}, 
               'axes': [], 'axesOrig': []}
    
    for i in range(nBoot):
        pre = megaTseps["preShift_late"][mod][anat][:,:,:,i]
        early  = megaTseps["early"][mod][anat][:,:,:,i]
        late  = megaTseps["late"][mod][anat][:,:,:,i]

        allD = [pre, early, late]#, post1BS, post2BS]
        for d in allD:
            d[np.isnan(d)] = 0
        
        res = eval(calcPhase).copy()
        group1 = res[0,:,:] 
        group2 = res[1,:,:]
        
        axesOrig = (group1 - group2)  
        norms = np.linalg.norm(axesOrig, axis=1, keepdims=True)
        axes = axesOrig / norms   
        
        allData['axesOrig'].append(axesOrig)
        allData['axes'].append(axes)    
        
        projs = {
            'pre': np.einsum('gtr,tr->gt', pre, axes),
            'early': np.einsum('gtr,tr->gt', early, axes),
            'late': np.einsum('gtr,tr->gt', late, axes)
        }
        
        for key, prj in projs.items():
            allData['projs'][key].append(prj)

    for key in ['pre', 'early','late']:
        allData['projs'][key] = np.stack(allData['projs'][key])
    
    allData['axesOrig'] = np.stack(allData['axesOrig'])
    allData['axes'] = np.stack(allData['axes'])
    
    return allData

def get_cTimes(mat, sub_goCue = False):
    cTimes = []
    nTrials = mat['nTrials']
    for t in range(nTrials):
        goQue = mat['RawEvents']['Trial'][t]['States']['AnswerPeriod'][0]
        oStateData = mat['RawData']['OriginalStateData'][t]
        if 12 in oStateData: # no respose
            respTime = np.nan
        else:
            if 10 in oStateData: #correct
                respIdx = (mat['RawData']['OriginalStateData'][t]==10).nonzero()[0][0]
                respTime = mat['RawData']['OriginalStateTimestamps'][t][respIdx]
            elif 13 in oStateData: #incorrect
                respIdx = (mat['RawData']['OriginalStateData'][t]==13).nonzero()[0][0]
                respTime = mat['RawData']['OriginalStateTimestamps'][t][respIdx]
            else:
                print('<<< problem in get_cTimes >>>')
        if sub_goCue:
            respTime = respTime-goQue
        cTimes.append(respTime)
    cTimes = np.array(cTimes)
    return cTimes
    
def get_cTimes2(behavior, sub_goCue = False):
    cTimes = np.hstack([get_cTimes(mat, sub_goCue = sub_goCue) for mat in behavior['mats']])
    return cTimes

def get_lickTimes(behavior):
    ports = ['Port2In', 'Port1In'] # IDX 0 = right, 1 = left. reorder these if they change
    bpodLicks = []

    sub_goCue = False
    for mat in behavior['mats']:
        nTrials = mat['nTrials']
        for t in range(nTrials):
            tLicks = np.ones([2,500])*np.nan #array to fill with lick timestamps and lick side
            lenS0 = 0
            trialStart_noMansLand = np.nanmax(mat['RawEvents']['Trial'][t]['States']['TrigTrialStart'][1].flatten())
            for s in range(2):  # side loop
                keys = mat['RawEvents']['Trial'][t]['Events'].keys()
                if ports[s] in keys:
                    sideLicks = mat['RawEvents']['Trial'][t]['Events'][ports[s]]
                    sideLicks = sideLicks
                    slType = type(sideLicks)
                    if s == 0:
                        if slType == list or slType == np.ndarray:
                            sideLicks = sideLicks[sideLicks>trialStart_noMansLand]
                            tLicks[0, :len(sideLicks)] = sideLicks  # timestamps
                            tLicks[1, :len(sideLicks)] = s  # side
                            lenS0 = len(sideLicks)
                        else:
                            if sideLicks<trialStart_noMansLand:
                                continue
                            else:
                                tLicks[0, 0] = sideLicks
                                tLicks[1, 0] = s
                                lenS0 = 1    
                    else:
                        if slType == list or slType == np.ndarray:
                            sideLicks = sideLicks[sideLicks>trialStart_noMansLand]
                            tLicks[0, lenS0:lenS0+len(sideLicks)] = sideLicks  # timestamps
                            tLicks[1, lenS0:lenS0+len(sideLicks)] = s  # side
                        else:
                            if sideLicks<trialStart_noMansLand:
                                continue
                            else:
                                tLicks[0, lenS0] = sideLicks
                                tLicks[1, lenS0] = s

            tLicks = tLicks[np.isfinite(tLicks)]
            tLicks = tLicks.reshape(int(tLicks.shape[0]/2),2, order = 'F')
            tLicks2 = np.ones([500,2])*np.nan
            tLicks2[:tLicks.shape[0],0] = tLicks[:,0]
            tLicks2[:tLicks.shape[0],1] = tLicks[:,1]
            tLicks2 = tLicks2[tLicks2[:,0].argsort()]
            bpodLicks.append(tLicks2)   

            goCue = mat['RawEvents']['Trial'][t]['States']['AnswerPeriod'][0]

    bpodLicks = np.dstack(bpodLicks) #lick (in order) x lick or lick code code (0,1) x trial
    return bpodLicks

def get_goCues(behavior):
    goCues = []
    for mat in behavior['mats']:
        nTrials = mat['nTrials']
        for t in range(nTrials):
            if 6 in mat['RawData']['OriginalStateData'][t]:
                tsIdx = np.where(6 == mat['RawData']['OriginalStateData'][t])[0][0]
                ts = mat['RawData']['OriginalStateTimestamps'][t][tsIdx]
                goCues.append(ts)
            else:
                goCues.append(np.nan)
    goCues = np.array(goCues)
    return goCues

def get_lastLick(behavior, firstBout = True, freqThresh = 0.25):
    goCues = get_goCues(behavior)
    lickTimes = get_lickTimes(behavior)
    lastLick = []
    for t in range(lickTimes.shape[2]):
        tlt = lickTimes[np.isfinite(lickTimes[:,0,t]),0,t]
        if np.sum(tlt>goCues[t])==0:
            lastLick.append(np.nan)
        else:
            if not firstBout:
                lastLick.append(tlt[-1])
            else:
                tlt = tlt[tlt>goCues[t]]
                boutBreaks = tlt[np.where(np.diff(tlt)>freqThresh)[0]]
                if len(boutBreaks)==0:
                    lastLick.append(tlt[-1])
                else:
                    firstBoutEnd = boutBreaks[0]
                    lastLick.append(firstBoutEnd)

    lastLick = np.hstack(lastLick)
    return lastLick
    
def get_autoWaterTrials(behavior):
    awKeyCode = [True, False]
    aw = []
    mats = behavior['mats']
    for mat in mats:
        aw_ = [awKeyCode[mat['TrialSettings'][i]['GUI']['Autowater']-1] for i in range(mat['nTrials'])]
        aw.append(aw_)
    aw = np.hstack(aw)
    return aw

def extract_camIDs(files):
    dirStructure = files[0]
    if '/' in dirStructure:
        camIDs = np.array([files[i].split('/')[-1].split('_')[0] for i in range(len(files))])
    else:
        camIDs = np.array([files[i].split('_')[0] for i in range(len(files))])
    camIDs = list(np.unique(camIDs))
    for i, camID in enumerate(camIDs):
        if 'B' in camID:
            camIDs.pop(i)
    return np.array(camIDs)
   
# Wrapper functions

def get_allAnmParams(masterData, anms, frs, factors, Z, phases2, dffKeys, maskKey = 'consensus_NMFtc_mask', windowSize = [-3.6,3.6]):
    """
    Build the per-animal parameter dictionary (behavior + event-aligned dF/F) used across the
    Scheib 2026 figures.

    For each animal this gathers trial timing (go cues, contact times, second licks), the tracker,
    lick trajectories / azimuth-elevation / LDA distances, motor-error masks, and — for every
    dffKey and anatomy (dendrites/somas) — the full and downsampled trial-separated traces plus
    contact-, go-cue-, and second-lick-aligned peri-stimulus windows. It also builds a boolean
    trial mask per learning phase in `phases2` under params['shiftMask'].

    VERBOSE / EXTENDED-PHASE VARIANT. Differs from get_allAnmParams2() in exactly two ways:
      1) It prints progress (per-anat/dffKey line, the relDays list, and the days chosen for the
         'late' / 'preShift_late' / 'preShift_early' phases). get_allAnmParams2() is quiet.
      2) It resolves eight extra phase names INLINE instead of deferring to get_behaviorShiftMask():
         'd1', 'd2', 'd3'                     -> tracker relDay == 1 / 2 / 3
         'shiftDay_pre-1'                     -> relDay 0, trials 51-99
         'shiftDay_post1' ... 'shiftDay_post5'-> trials 101-200 on relDay 0 / 1 / 2 / 3 / 4
         Note that 'shiftDay_post1/2/3' ALSO exist in get_behaviorShiftMask() but mean something
         different there (thirds of the trials within relDay 0 only), so those three names give
         DIFFERENT masks in the two functions. The other five names are unknown to
         get_behaviorShiftMask() and come back as all-False (with an error print) in
         get_allAnmParams2().
    Everything else — all returned keys and their values — is identical between the two. For the
    standard manuscript `phases2` list (which contains none of those eight names) the two
    functions return the same result; use this one if you need the extra phases or the printouts,
    and get_allAnmParams2() otherwise. See also update_allAnmParams() to refresh only the traces.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anms = list(masterData.keys())          # animal IDs; keys into masterData
    frs = [9.35, 4.68]                      # frame rate in Hz, ordered [dendrites, somas]
    factors = [2, 1]                        # downsample factor (int >= 1), ordered [dendrites, somas]
    Z = False                               # True | False  (z-score traces inside get_tsep)
    phases2 = ['pre', 'preShift_early', 'preShift_late', 'post', 'late', 'day-1_beg', 'day-1_mid', 'day-1_end', 'shiftDay_post', 'shiftDay_pre', 'early']  # any get_behaviorShiftMask() phase, plus 'd1'|'d2'|'d3'|'shiftDay_pre-1'|'shiftDay_post1'-'shiftDay_post5' handled inline here
    dffKeys = ['dff_NMF']                   # dF/F variants to pull from masterData[anm][anat]
    maskKey = 'consensus_NMFtc_mask'        # frame-mask key passed to get_tsep(nanMaskFrs=True)
    windowSize = [-3.6, 3.6]                # peri-event window in seconds, [pre, post]
    importlib.reload(js_manuscript_final)
    allAnmParams = js_manuscript_final.get_allAnmParams(
        masterData, anms=anms, frs=frs, factors=factors, Z=Z, phases2=phases2, dffKeys=dffKeys, maskKey=maskKey, windowSize=windowSize)

    Parameters
    ----------
    masterData : dict
        Master data dictionary keyed by animal ID; each entry holds 'behavior', 'dendrites', 'somas'.
    anms : list of str
        Animal IDs to process (keys of masterData).
    frs : list of float
        Frame rates in Hz, ordered [dendrites, somas].
    factors : list of int
        Temporal downsample factors, ordered [dendrites, somas].
    Z : bool
        Whether get_tsep() z-scores the traces.
    phases2 : list of str
        Learning-phase names to build boolean trial masks for; stored in params['shiftMask'].
    dffKeys : list of str
        dF/F keys to extract trial-separated traces for.
    maskKey : str
        Frame-mask key used by get_tsep() to NaN bad frames.
    windowSize : list of float
        Peri-event window [pre, post] in seconds for the event-aligned traces.

    Returns
    -------
    allAnmParams : dict
        {animal ID: params dict} where params holds behavior arrays, 'shiftMask', and one entry per
        dffKey containing 'fullTsep' and 'dsTsep' sub-dicts per anatomy.
    """
    allAnmParams = {}
    for anmIdx,anm in enumerate(anms):
        print("=========================================================================================================")
        print(f'{anmIdx} {anm}')

        params = {}
        anmData = masterData[anm]
        beh = anmData['behavior']
        
        params['zavail'] = np.array(anmData['behavior']['twoCams'])
        dendDays = get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = get_goCues(anmData['behavior'])
        cTimes = get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = get_lastLick(anmData['behavior'])
        params['lastLicks'] = lastLicks
        slSides, slTimes = get_secondLickSide(beh, returnTime = True)
        params['sTimes'] = slTimes
        slTimes2 = slTimes+goCues

        tracker = get_tracker(anmData['behavior'])
        params['tracker'] = tracker
        
        try:
            params['aw'] = get_autoWaterTrials(anmData['behavior'])
        except:
            params['aw'] = np.zeros(tracker.shape[0]).astype(bool)

        for dffKey in dffKeys:
            params2 = {}
            params2['dsTsep'] = {}
            params2['fullTsep'] = {}
            for anat in anmA:
                print(f'Getting params for {anat} {dffKey}')
                if anat == 'dendrites':
                    anatMask = dendDays
                    ak = 0
                else:
                    anatMask = ~dendDays
                    ak = 1
                fr = frs[ak]
                factor = factors[ak]
                
                preTsep = get_tsep(anmData[anat], dffKey, maskKey = maskKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
                preTsep = preTsep[:,:,int(fr*3):]
                params2['fullTsep'][anat] = {}
                params2['fullTsep'][anat]['full'] = preTsep
                
                ctps = psTsep_func(preTsep, cTimes[anatMask]+1, fr, windowSize = windowSize)#, downSample = True, downSampleFactor = factors[ak])
                gops = psTsep_func(preTsep, goCues[anatMask]+1, fr, windowSize = windowSize)#, downSample = True, downSampleFactor = factors[ak])
                sps = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, windowSize = windowSize)
                params2['fullTsep'][anat]['ct-psTsep'] = ctps
                params2['fullTsep'][anat]['go-psTsep'] = gops
                params2['fullTsep'][anat]['sl-psTsep'] = sps

                dsPre = downSampleTsep(preTsep, factor = factor)
                ctps2 = psTsep_func(preTsep, cTimes[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                gops2 = psTsep_func(preTsep, goCues[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                stps2 = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                
                params2['dsTsep'][anat] = {}
                params2['dsTsep'][anat]['full'] = dsPre
                params2['dsTsep'][anat]['ct-psTsep'] = ctps2
                params2['dsTsep'][anat]['go-psTsep'] = gops2
                params2['dsTsep'][anat]['sl-psTsep'] = stps2
                # params2['dsTsep'][anat]['preTsep'] = dsPre
                
            params[dffKey] = params2

        cLicks = get_contacts(anmData['behavior']['trajs'])
        cazi,cele = get_aziEle(anmData['behavior'], cLicks)
        cdists, lda = get_relDistances(cazi, cele, tracker, returnLDA = True)
        params['cLicks'] = cLicks
        params['cDists'] = cdists
        params['cAzi'] = cazi
        params['cEle'] = cele
        params['lda'] = lda
        

        fLicks = get_firsts(anmData['behavior']['trajs'])
        fazi,fele = get_aziEle(anmData['behavior'], fLicks)
        fdists2 = get_relDistances(fazi, fele, tracker) ########################## need to check if lda is redefined here
        fdists = get_relDistFromLDA(fazi, fele, lda)
        params['fLicks'] = fLicks
        params['fDists_orig'] = fdists2
        params['fDists'] = fdists
        params['fAzi'] = fazi
        params['fEle'] = fele



        params['shiftMask'] = {}
        # for phase in phases2:
        #     params['shiftMask'][phase] = get_behaviorShiftMask(tracker, phase)
        relDays = np.unique(tracker[:,0])
        print(f'relDays: {relDays}')
        for phase in phases2:
            if phase == 'd1':
                shiftMask = tracker[:, 0] == 1
                params['shiftMask'][phase] = shiftMask
            elif phase == 'd2':
                shiftMask = tracker[:, 0] == 2
                params['shiftMask'][phase] = shiftMask
            elif phase == 'd3':
                shiftMask = tracker[:, 0] == 3
                params['shiftMask'][phase] = shiftMask

            elif phase == 'shiftDay_pre-1':
                t1 = tracker[:,1]>50
                t2 = tracker[:,1]<=99
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==0, m1)
                params['shiftMask'][phase] = shiftMask
            elif phase == 'shiftDay_post1':
                t1 = tracker[:,1]>100
                t2 = tracker[:,1]<=200
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==0, m1)
                params['shiftMask'][phase] = shiftMask
            elif phase == 'shiftDay_post2':
                t1 = tracker[:,1]>100
                t2 = tracker[:,1]<=200
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==1, m1)
                params['shiftMask'][phase] = shiftMask
            elif phase == 'shiftDay_post3':
                t1 = tracker[:,1]>100
                t2 = tracker[:,1]<=200
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==2, m1)
                params['shiftMask'][phase] = shiftMask
            elif phase == 'shiftDay_post4':
                t1 = tracker[:,1]>100
                t2 = tracker[:,1]<=200
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==3, m1)
                params['shiftMask'][phase] = shiftMask
            elif phase == 'shiftDay_post5':
                t1 = tracker[:,1]>100
                t2 = tracker[:,1]<=200
                m1 = np.logical_and(t1,t2)
                shiftMask = np.logical_and(tracker[:,0]==4, m1)
                params['shiftMask'][phase] = shiftMask

            elif phase == 'late':
                try:
                    relDays = np.unique(tracker[:,0])
                    shiftIdx = np.where(relDays>=0)[0][0]
                    shiftMask = tracker[:, 0] >= relDays[shiftIdx+2]
                    if len(np.unique(tracker[shiftMask,0]))>2:
                        keepDays = np.unique(tracker[shiftMask,0])[:2]
                        print(f'phase: {phase}, days used: {keepDays}')

                        shiftMask = np.isin(tracker[:,0], keepDays)
                except:
                    shiftMask = np.zeros(tracker.shape[0]).astype(bool)
                    print('empty fail no late post')
                params['shiftMask'][phase] = shiftMask
                
            elif phase == 'preShift_late':
                availDays = relDays[relDays<=0]
                if len(availDays)>2:
                    days = availDays[-2:]
                    print(f'phase: {phase}, days used: {days}')

                    if 0 in days:
                        dayBeforeMask = tracker[:,0]==availDays[-2]
                        shiftDayMask = np.logical_and(tracker[:,0]==0, tracker[:,1]<100)
                        shiftMask = np.logical_or(dayBeforeMask, shiftDayMask)
                    else:
                        shiftMask = np.isin(tracker[:,0], days)
                else:
                    shiftMask = np.isin(tracker[:,0], availDays)
                params['shiftMask'][phase] = shiftMask
                
            elif phase == 'preShift_early':
                availDays = relDays[relDays<=0]
                if len(availDays)>2:
                    availDays2 = availDays[:-2]
                    if len(availDays2)>1:
                        days = availDays2[-2:]
                        print(f'phase: {phase}, days used: {days}')

                        shiftMask = np.isin(tracker[:,0], days)
                    else:
                        shiftMask = np.isin(tracker[:,0], availDays2)
                else:
                    print('not enough days for preShift_early')
                    shiftMask = np.zeros(tracker.shape[0]).astype(bool)
                params['shiftMask'][phase] = shiftMask
                    
            else:
                params['shiftMask'][phase] = get_behaviorShiftMask(tracker, phase)

        
        me, de = get_errors_LDA(cazi, cele, tracker, returnLDA=False)
        minDists = get_minDistToPort(beh, 'contacts', 'post_shift', 'left')
        params['me'] = me
        params['de'] = de
        params['minDist'] = minDists

        ar, al, an = get_MER_bySecondLick(beh)
        params['afterRight'] = ar
        params['afterLeft'] = al
        params['afterNone'] = an
        
        sLicks = get_consumption_licks(beh['trajs'],0)
        aziS,eleS = get_aziEle(beh, sLicks)
        params['sLicks'] = sLicks
        params['sAzi'] = aziS
        params['sEle'] = eleS
        meS, deS = get_sl_errors(aziS, eleS, np.logical_and(me, tracker[:,2]==5), lda, slSides)
        params['MER-L_ME'] = np.logical_and(np.logical_and(me, tracker[:,2]==5), meS)
        params['MER-L_noME'] = np.logical_and(np.logical_and(me, tracker[:,2]==5), deS)

        slSide, slN = get_secondLickSide2(beh, returnTime = False, returnN = True, returnNextIdent = False)
        params['sl_nHits'] = slN
        params['slSideBpod'] = slSide

        coorsS = np.vstack([aziS,eleS]).T
        nanMask = np.any(np.isfinite(coorsS),axis = 1)
        sCoors = coorsS[nanMask,:]
        
        pred = lda.predict(sCoors)

        actual = np.ones(nanMask.shape[0])*np.nan
        c = 0
        for t in range(len(actual)):
            if nanMask[t]:
                actual[t] = pred[c]
                c+=1

        params['slSideTracking'] = actual
        sdists = get_relDistFromLDA(aziS, eleS, lda)
        params['sDists'] = sdists
        
        allAnmParams[anm] = params

    return allAnmParams

def get_allAnmParams2(masterData, anms, frs, factors, Z, phases2, dffKeys, maskKey = 'consensus_NMFtc_mask', windowSize = [-3.6,3.6]):
    """
    Build the per-animal parameter dictionary (behavior + event-aligned dF/F) used across the
    Scheib 2026 figures.

    For each animal this gathers trial timing (go cues, contact times, second licks), the tracker,
    lick trajectories / azimuth-elevation / LDA distances, motor-error masks, and — for every
    dffKey and anatomy (dendrites/somas) — the full and downsampled trial-separated traces plus
    contact-, go-cue-, and second-lick-aligned peri-stimulus windows. It also builds a boolean
    trial mask per learning phase in `phases2` under params['shiftMask'].

    QUIET / STANDARD-PHASE VARIANT. Differs from get_allAnmParams() in exactly two ways:
      1) It does not print per-anat/dffKey progress, the relDays list, or the days chosen for the
         'late' / 'preShift_late' / 'preShift_early' phases. get_allAnmParams() prints all of those.
      2) It handles only 'late', 'preShift_late' and 'preShift_early' inline and sends every other
         phase name straight to get_behaviorShiftMask(). get_allAnmParams() additionally resolves
         'd1', 'd2', 'd3', 'shiftDay_pre-1' and 'shiftDay_post1'-'shiftDay_post5' inline. Passing
         any of those five first names here yields an all-False mask plus a
         '<<<Error: ... learning phase not found>>>' print, and passing 'shiftDay_post1/2/3' here
         yields get_behaviorShiftMask()'s definition (thirds of the trials within relDay 0), which
         is NOT the same mask get_allAnmParams() would build for those names.
    Everything else — all returned keys and their values — is identical between the two. For the
    standard manuscript `phases2` list (which contains none of those eight names) the two
    functions return the same result; prefer this one unless you need the extra phases or the
    printouts. See also update_allAnmParams() to refresh only the traces.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anms = list(masterData.keys())          # animal IDs; keys into masterData
    frs = [9.35, 4.68]                      # frame rate in Hz, ordered [dendrites, somas]
    factors = [2, 1]                        # downsample factor (int >= 1), ordered [dendrites, somas]
    Z = False                               # True | False  (z-score traces inside get_tsep)
    phases2 = ['pre', 'preShift_early', 'preShift_late', 'post', 'late', 'day-1_beg', 'day-1_mid', 'day-1_end', 'shiftDay_post', 'shiftDay_pre', 'early']  # any get_behaviorShiftMask() phase; the extra 'd1'|'d2'|'d3'|'shiftDay_pre-1'|'shiftDay_post4'|'shiftDay_post5' names are NOT supported here
    dffKeys = ['dff_NMF']                   # dF/F variants to pull from masterData[anm][anat]
    maskKey = 'consensus_NMFtc_mask'        # frame-mask key passed to get_tsep(nanMaskFrs=True)
    windowSize = [-3.6, 3.6]                # peri-event window in seconds, [pre, post]
    importlib.reload(js_manuscript_final)
    allAnmParams = js_manuscript_final.get_allAnmParams2(
        masterData, anms=anms, frs=frs, factors=factors, Z=Z, phases2=phases2, dffKeys=dffKeys, maskKey=maskKey, windowSize=windowSize)

    Parameters
    ----------
    masterData : dict
        Master data dictionary keyed by animal ID; each entry holds 'behavior', 'dendrites', 'somas'.
    anms : list of str
        Animal IDs to process (keys of masterData).
    frs : list of float
        Frame rates in Hz, ordered [dendrites, somas].
    factors : list of int
        Temporal downsample factors, ordered [dendrites, somas].
    Z : bool
        Whether get_tsep() z-scores the traces.
    phases2 : list of str
        Learning-phase names to build boolean trial masks for; stored in params['shiftMask'].
    dffKeys : list of str
        dF/F keys to extract trial-separated traces for.
    maskKey : str
        Frame-mask key used by get_tsep() to NaN bad frames.
    windowSize : list of float
        Peri-event window [pre, post] in seconds for the event-aligned traces.

    Returns
    -------
    allAnmParams : dict
        {animal ID: params dict} where params holds behavior arrays, 'shiftMask', and one entry per
        dffKey containing 'fullTsep' and 'dsTsep' sub-dicts per anatomy.
    """
    allAnmParams = {}
    for anmIdx,anm in enumerate(anms):
        print("=========================================================================================================")
        print(f'{anmIdx} {anm}')

        params = {}
        anmData = masterData[anm]
        beh = anmData['behavior']
        
        params['zavail'] = np.array(anmData['behavior']['twoCams'])
        dendDays = get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = get_goCues(anmData['behavior'])
        cTimes = get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = get_lastLick(anmData['behavior'])
        params['lastLicks'] = lastLicks
        slSides, slTimes = get_secondLickSide(beh, returnTime = True)
        params['sTimes'] = slTimes
        slTimes2 = slTimes+goCues

        tracker = get_tracker(anmData['behavior'])
        params['tracker'] = tracker

        try:
            params['aw'] = get_autoWaterTrials(anmData['behavior'])
        except:
            params['aw'] = np.zeros(tracker.shape[0]).astype(bool)

        for dffKey in dffKeys:
            params2 = {}
            params2['dsTsep'] = {}
            params2['fullTsep'] = {}
            for anat in anmA:
                if anat == 'dendrites':
                    anatMask = dendDays
                    ak = 0
                else:
                    anatMask = ~dendDays
                    ak = 1
                fr = frs[ak]
                factor = factors[ak]
                
                preTsep = get_tsep(anmData[anat], dffKey, maskKey = maskKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
                preTsep = preTsep[:,:,int(fr*3):]
                params2['fullTsep'][anat] = {}
                params2['fullTsep'][anat]['full'] = preTsep
                
                ctps = psTsep_func(preTsep, cTimes[anatMask]+1, fr, windowSize = windowSize)#, downSample = True, downSampleFactor = factors[ak])
                gops = psTsep_func(preTsep, goCues[anatMask]+1, fr, windowSize = windowSize)#, downSample = True, downSampleFactor = factors[ak])
                sps = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, windowSize = windowSize)
                params2['fullTsep'][anat]['ct-psTsep'] = ctps
                params2['fullTsep'][anat]['go-psTsep'] = gops
                params2['fullTsep'][anat]['sl-psTsep'] = sps

                dsPre = downSampleTsep(preTsep, factor = factor)
                ctps2 = psTsep_func(preTsep, cTimes[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                gops2 = psTsep_func(preTsep, goCues[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                stps2 = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                
                params2['dsTsep'][anat] = {}
                params2['dsTsep'][anat]['full'] = dsPre
                params2['dsTsep'][anat]['ct-psTsep'] = ctps2
                params2['dsTsep'][anat]['go-psTsep'] = gops2
                params2['dsTsep'][anat]['sl-psTsep'] = stps2
                # params2['dsTsep'][anat]['preTsep'] = dsPre
                
            params[dffKey] = params2

        cLicks = get_contacts(anmData['behavior']['trajs'])
        cazi,cele = get_aziEle(anmData['behavior'], cLicks)
        cdists, lda = get_relDistances(cazi, cele, tracker, returnLDA = True)
        params['cLicks'] = cLicks
        params['cDists'] = cdists
        params['cAzi'] = cazi
        params['cEle'] = cele
        params['lda'] = lda
        

        fLicks = get_firsts(anmData['behavior']['trajs'])
        fazi,fele = get_aziEle(anmData['behavior'], fLicks)
        fdists2 = get_relDistances(fazi, fele, tracker) ########################## need to check if lda is redefined here
        fdists = get_relDistFromLDA(fazi, fele, lda)
        params['fLicks'] = fLicks
        params['fDists_orig'] = fdists2
        params['fDists'] = fdists
        params['fAzi'] = fazi
        params['fEle'] = fele



        params['shiftMask'] = {}
        # for phase in phases2:
        #     params['shiftMask'][phase] = get_behaviorShiftMask(tracker, phase)
        relDays = np.unique(tracker[:,0])
        for phase in phases2:
            if phase == 'late':
                try:
                    relDays = np.unique(tracker[:,0])
                    shiftIdx = np.where(relDays>=0)[0][0]
                    shiftMask = tracker[:, 0] >= relDays[shiftIdx+2]
                    if len(np.unique(tracker[shiftMask,0]))>2:
                        keepDays = np.unique(tracker[shiftMask,0])[:2]
                        shiftMask = np.isin(tracker[:,0], keepDays)
                except:
                    shiftMask = np.zeros(tracker.shape[0]).astype(bool)
                    print('empty fail no late post')
                params['shiftMask'][phase] = shiftMask
                
            elif phase == 'preShift_late':
                availDays = relDays[relDays<=0]
                if len(availDays)>2:
                    days = availDays[-2:]
                    if 0 in days:
                        dayBeforeMask = tracker[:,0]==availDays[-2]
                        shiftDayMask = np.logical_and(tracker[:,0]==0, tracker[:,1]<100)
                        shiftMask = np.logical_or(dayBeforeMask, shiftDayMask)
                    else:
                        shiftMask = np.isin(tracker[:,0], days)
                else:
                    shiftMask = np.isin(tracker[:,0], availDays)
                params['shiftMask'][phase] = shiftMask
                
            elif phase == 'preShift_early':
                availDays = relDays[relDays<=0]
                if len(availDays)>2:
                    availDays2 = availDays[:-2]
                    if len(availDays2)>1:
                        days = availDays2[-2:]
                        shiftMask = np.isin(tracker[:,0], days)
                    else:
                        shiftMask = np.isin(tracker[:,0], availDays2)
                else:
                    print('not enough days for preShift_early')
                    shiftMask = np.zeros(tracker.shape[0]).astype(bool)
                params['shiftMask'][phase] = shiftMask
                    
            else:
                params['shiftMask'][phase] = get_behaviorShiftMask(tracker, phase)

        
        me, de = get_errors_LDA(cazi, cele, tracker, returnLDA=False)
        minDists = get_minDistToPort(beh, 'contacts', 'post_shift', 'left')
        params['me'] = me
        params['de'] = de
        params['minDist'] = minDists

        ar, al, an = get_MER_bySecondLick(beh)
        params['afterRight'] = ar
        params['afterLeft'] = al
        params['afterNone'] = an
        
        sLicks = get_consumption_licks(beh['trajs'],0)
        aziS,eleS = get_aziEle(beh, sLicks)
        params['sLicks'] = sLicks
        params['sAzi'] = aziS
        params['sEle'] = eleS
        meS, deS = get_sl_errors(aziS, eleS, np.logical_and(me, tracker[:,2]==5), lda, slSides)
        params['MER-L_ME'] = np.logical_and(np.logical_and(me, tracker[:,2]==5), meS)
        params['MER-L_noME'] = np.logical_and(np.logical_and(me, tracker[:,2]==5), deS)

        slSide, slN = get_secondLickSide2(beh, returnTime = False, returnN = True, returnNextIdent = False)
        params['sl_nHits'] = slN
        params['slSideBpod'] = slSide

        coorsS = np.vstack([aziS,eleS]).T
        nanMask = np.any(np.isfinite(coorsS),axis = 1)
        sCoors = coorsS[nanMask,:]
        
        pred = lda.predict(sCoors)

        actual = np.ones(nanMask.shape[0])*np.nan
        c = 0
        for t in range(len(actual)):
            if nanMask[t]:
                actual[t] = pred[c]
                c+=1

        params['slSideTracking'] = actual
        sdists = get_relDistFromLDA(aziS, eleS, lda)
        params['sDists'] = sdists
        
        allAnmParams[anm] = params

    return allAnmParams

def update_allAnmParams(allAnmParams, masterData, anms, frs, factors, Z, dffKeys, maskKey = 'consensus_NMFtc_mask',  windowSize = [-3.6,3.6]):

    for anm in anms:
        
        params = {}
        anmData = masterData[anm]

        dendDays = get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = get_goCues(anmData['behavior'])
        cTimes = get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = get_lastLick(anmData['behavior'])
        params['lastLicks'] = lastLicks
        slSides, slTimes = get_secondLickSide(anmData['behavior'], returnTime = True)
        params['sTimes'] = slTimes
        slTimes2 = slTimes+goCues

        for dffKey in dffKeys:
            params2 = {}
            params2['dsTsep'] = {}
            #params2['fullTsep'] = {}
            for anat in anmA:
                if anat == 'dendrites':
                    anatMask = dendDays
                    ak = 0
                else:
                    anatMask = ~dendDays
                    ak = 1
                fr = frs[ak]
                factor = factors[ak]
                
                preTsep = get_tsep(anmData[anat], dffKey, maskKey = maskKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
                preTsep = preTsep[:,:,int(fr*3):]
                dsPre = downSampleTsep(preTsep, factor = factor)
                ctps2 = psTsep_func(preTsep, cTimes[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                gops2 = psTsep_func(preTsep, goCues[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                stps2 = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, downSample = True, downSampleFactor = factor, windowSize = windowSize)
                
                params2['dsTsep'][anat] = {}
                params2['dsTsep'][anat]['full'] = dsPre
                params2['dsTsep'][anat]['ct-psTsep'] = ctps2
                params2['dsTsep'][anat]['go-psTsep'] = gops2
                params2['dsTsep'][anat]['sl-psTsep'] = stps2
                # params2['dsTsep'][anat]['preTsep'] = dsPre
                
            params[dffKey] = params2
        allAnmParams[anm][dffKey]['dsTsep'] = params[dffKey]['dsTsep']
    return allAnmParams

def check_anms_ROIs(allAnmParams, anms, anatKeys, dffKey, thrsh = 0.5):

    ##################################################################################################################################

    aAnms = {}
    for anat in anatKeys:
        aAnms[anat] = []
        for anm in anms:
            if anat in list(allAnmParams[anm][dffKey]['dsTsep'].keys()):
                aAnms[anat].append(anm)
    ##################################################################################################################################

    keepRoisFull = {}
    keepRoisDS = {}
    allNaNs = {}
    #ax = plt.figure().add_subplot()
    for anat in anatKeys:
        keepRoisFull[anat] = {}
        keepRoisDS[anat] = {}
        allNaNs[anat] = []
        for anm in aAnms[anat]:
            # anmData = masterData[anm][anat]
            # goodRois = getGoodRois(anmData, dffKey = 'NMF')
            # masks = anmData['sessions'][0]['consensus_NMF_masks'][:,:,goodRois]
            #ftsep = allAnmParams[anm][dffKey]['fullTsep'][anat]['full']
            dsTsep = allAnmParams[anm][dffKey]['dsTsep'][anat]['full']

            for ts, tsep in enumerate([dsTsep]):
                kr = []
                for roi in range(tsep.shape[1]):
                    nNaN = np.sum(np.isnan(tsep[:,roi,:]))
                    nTime = tsep.shape[0]*tsep.shape[2]
                    if nNaN/nTime<thrsh:
                        #kr.append(True)
                        ###################################################################################
                        # if anat == 'dendrites' and anm == 'B00002213997' and roi == 1:
                        #     kr.append(False)
                        # elif anat == 'dendrites' and anm == 'B00002213997' and roi == 0:
                        #     kr.append(False)
                        # else:
                            kr.append(True)
                        ###################################################################################
                    else:
                        kr.append(False)
                # if ts == 0:
                #     keepRoisFull[anat][anm] = np.hstack(kr)
                #else:
                keepRoisDS[anat][anm] = np.hstack(kr)

    for i, j in enumerate([keepRoisDS]):
        print(f'{["full", "downSampled"][i]}')
        for ak, anat in enumerate(anatKeys):
            sm = np.sum(np.hstack([j[anat][anm] for anm in list(j[anat].keys())]))
            total = np.sum(np.hstack([len(j[anat][anm]) for anm in list(j[anat].keys())]))
            print(anat, sm,'/',total, np.round(sm/total, 3))
        print('')

    ##################################################################################################################################
    return aAnms, keepRoisFull, keepRoisDS, allNaNs

def check_trials(allAnmParams, anms, anatKeys, dffKey, minTrials = 5):


    iAnms = {}
    aAnms = {}
    
    allMinTrials = []
    anmsMinTrial = {}
    for anat in anatKeys:
        iAnms[anat] = []
        aAnms[anat] = []
        anmsMinTrial[anat] = {}
        for anm in anms:
            if anat in list(allAnmParams[anm][dffKey]['dsTsep'].keys()):
                aAnms[anat].append(anm)
                tracker = allAnmParams[anm]['tracker']
                dendDays = allAnmParams[anm]['dendDays']
                preMask = allAnmParams[anm]['shiftMask']['pre']
                
                if anat == 'dendrites':
                    anatMask = dendDays
                else:
                    anatMask = ~dendDays

                preAnatMask = np.logical_and(anatMask, preMask)
                preCR = np.sum(tracker[preAnatMask,2]==0)
                preCL = np.sum(tracker[preAnatMask,2]==1)
                preFR = np.sum(tracker[preAnatMask,2]==5)
                preFL = np.sum(tracker[preAnatMask,2]==6)
                dsTsep = allAnmParams[anm][dffKey]['dsTsep'][anat]['full']

                tsep = dsTsep[0]

                if np.all(np.hstack([item>minTrials for item in [preCR, preCL, preFR, preFL]])):
                    if anm != 'B00002121777':
                        if tsep.shape[1]>=10:
                            iAnms[anat].append(anm)
                            print(anat, anm, preCR, preCL, preFR, preFL)
                            allMinTrials.append(np.hstack([preCR, preCL, preFR, preFL]))
                            anmsMinTrial[anat][anm] = np.min(np.hstack([preCR, preCL, preFR, preFL]))
    print(np.min(np.hstack(allMinTrials)))
    print("Finished processing trials.")
    ##################################################################################################################################

    return aAnms, iAnms, anmsMinTrial, allMinTrials

def full_get_SMO(anatKeys, mods, phases, allAnmParams, anmsToUse, dffKey, alignTimes, frs, keepRoisDS, label, preComp, nBoot = 1000, shuffAnms = True,
             binary = False, enforceMinTrial = True, anmsMinTrial = None, preCompHH = None):
    """
    Bootstrap-assemble the population sensory/motor/outcome (SMO) coding-direction
    (CD) projection tensors across animals, anatomies, task phases, and alignments.

    For each bootstrap iteration animals are optionally resampled (with
    replacement) per anatomy, their per-animal precomputed CD projections
    (`preComp`, produced by `pullAnmTsepComputeSMOXval`) are concatenated across
    ROIs via `gather_anmShuff_megaTsepFullSMO`, and the result is written into
    per-CD "mega" tensors. If `preCompHH` is supplied, the matched half-split
    (50:50) projections are also gathered with the same animal resampling.

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    anatKeys = ['dendrites', 'somas']   # anatomy keys; index 0 dendrites, 1 somas
    mods = ['GO', 'CT']                 # per-alignment (and error) shorthand; 'err' in a mod -> returnErrors
    phases = ['pre']                    # task phase keys to gather
    anmsToUse = anmsToUse               # dict[anat] -> list of animal ids to include
    dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events'  # key selecting the trace in allAnmParams[anm]
    alignTimes = ['goCues', 'cTimes']   # alignment event per mod (same length/order as mods)
    frs = [frDend, frSoma]              # per-anatomy frame rate (Hz), indexed like anatKeys
    keepRoisDS = keepRoisDS             # dict[anat][anm] -> ROI indices (downsampled data)
    label = ['s', 'm', 'o']            # CD keys: sensory / motor / outcome
    preComp = preComp                   # nested precomputed CD projections keyed [anat][anm][phase][alignT][cdKey]
    nBoot = 1000                        # number of bootstrap iterations
    shuffAnms = True                    # True -> resample animals with replacement each boot; False -> fixed order
    binary = False                      # forwarded to the gather funcs (binarize events)
    enforceMinTrial = True              # True -> drop animals below the min-trial threshold
    anmsMinTrial = None                 #     > only if enforceMinTrial=True: dict[anat] -> animals passing the min-trial cutoff
    preCompHH = None                    # None -> skip half-split; or HH precomputed projections keyed like preComp with (half1, half2)
    importlib.reload(js_manuscript)
    megaTsepsS, megaTsepsM, megaTsepsO, megaTsepsHH, megaTsepsHHS, megaTsepsHHM, megaTsepsHHO = js_manuscript.full_get_SMO(
        anatKeys, mods=mods, phases=phases, allAnmParams=allAnmParams, anmsToUse=anmsToUse,
        dffKey=dffKey, alignTimes=alignTimes, frs=frs, keepRoisDS=keepRoisDS, label=label,
        preComp=preComp, nBoot=nBoot, shuffAnms=shuffAnms, binary=binary,
        enforceMinTrial=enforceMinTrial, anmsMinTrial=anmsMinTrial, preCompHH=preCompHH)

    Parameters
    ----------
    anatKeys : list of str
        Anatomy keys to loop over (e.g. ['dendrites', 'somas']); the index into
        this list selects the matching `frs` entry and holder tensor.
    mods : list of str
        Per-alignment shorthand labels (e.g. ['GO', 'CT']); a mod containing
        'err' turns on error-trial output (`returnErrors=True`) in the gather.
    phases : list of str
        Task phase keys to gather; also used to key the returned dicts (and the
        `{phase}1`/`{phase}2` half-split keys).
    allAnmParams : dict
        Per-animal parameter/data dict providing the time-separated traces used
        to size the holder tensors.
    anmsToUse : dict
        `anmsToUse[anat]` -> list of animal ids contributing to that anatomy.
    dffKey : str
        Key selecting the dff/deconvolved-events trace within `allAnmParams[anm]`.
    alignTimes : list of str
        Alignment event per mod (e.g. 'goCues', 'cTimes'); same length/order as
        `mods`.
    frs : list of float
        Per-anatomy imaging frame rate (Hz), indexed like `anatKeys`.
    keepRoisDS : dict
        `keepRoisDS[anat][anm]` -> ROI indices (boolean or integer) to keep for
        the downsampled traces; used to size and concatenate the ROI axis.
    label : list of str
        CD keys to gather, e.g. ['s', 'm', 'o'] for sensory/motor/outcome.
    preComp : dict
        Precomputed per-animal CD projections keyed `[anat][anm][phase][alignT][cdKey]`
        (from `pullAnmTsepComputeSMOXval`).
    nBoot : int
        Number of bootstrap iterations (final axis length of every holder tensor).
    shuffAnms : bool
        True to resample animals with replacement each iteration; False to use a
        fixed animal order (no resampling).
    binary : bool
        Forwarded to the gather functions (binarize events before projecting).
    enforceMinTrial : bool
        True to drop animals that fail the per-animal minimum-trial threshold.
    anmsMinTrial : dict or None
        Only used when `enforceMinTrial` is True. `anmsMinTrial[anat]` -> the
        animals passing the min-trial cutoff.
    preCompHH : dict or None
        None to skip the half-split output; otherwise HH precomputed projections
        keyed like `preComp` with each leaf a (half1, half2) pair, gathered with
        the same animal resampling as the full-SMO path.

    Returns
    -------
    megaTsepsS, megaTsepsM, megaTsepsO : dict
        Bootstrapped SMO tensors for the sensory / motor / outcome CDs, each keyed
        `[phase][mod][anat]` with array shape (4, nTime, nRoi, nBoot) over the four
        trial types [CR, CL, IR, IL].
    megaTsepsHH : dict
        Combined half-split tensors keyed `[f'{phase}1' | f'{phase}2'][mod][anat]`
        (empty/unfilled when `preCompHH` is None).
    megaTsepsHHS, megaTsepsHHM, megaTsepsHHO : dict
        Per-CD (s/m/o) half-split tensors, keyed the same as `megaTsepsHH`.
    """
    print("Preparing SMO holders...")
    megaTsepsS = {}
    megaTsepsM = {}
    megaTsepsO = {}

    megaTsepsHH = {}
    # half-split SMO holders, one per cdKey (s/m/o), each keyed {phase}1/{phase}2
    megaTsepsHHS = {}
    megaTsepsHHM = {}
    megaTsepsHHO = {}
    megaTsepsHHByKey = {'s': megaTsepsHHS, 'm': megaTsepsHHM, 'o': megaTsepsHHO}
    emptyMegaTsep = {'dendrites': [], 'somas': []}

    for p, phase in enumerate(phases):
        megaTsepsS[phase] = {}
        megaTsepsM[phase] = {}
        megaTsepsO[phase] = {}
        megaTsepsHH[f'{phase}1'] = {}
        megaTsepsHH[f'{phase}2'] = {}
        for hh in megaTsepsHHByKey.values():
            hh[f'{phase}1'] = {}
            hh[f'{phase}2'] = {}
        for mod in mods:
            megaTsepsS[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsM[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsO[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsHH[f'{phase}1'][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsHH[f'{phase}2'][mod] = copy.deepcopy(emptyMegaTsep)
            for hh in megaTsepsHHByKey.values():
                hh[f'{phase}1'][mod] = copy.deepcopy(emptyMegaTsep)
                hh[f'{phase}2'][mod] = copy.deepcopy(emptyMegaTsep)

    print("Making holder matrices...")
    somaTime = allAnmParams[anmsToUse['somas'][0]][dffKey]['dsTsep']['somas']['go-psTsep'].shape[2]
    nSoma = 0
    for anm in anmsToUse['somas']:
        nSoma+=allAnmParams[anm][dffKey]['dsTsep']['somas']['go-psTsep'][:,keepRoisDS['somas'][anm],:].shape[1]
        
    dendTime = allAnmParams[anmsToUse['dendrites'][0]][dffKey]['dsTsep']['dendrites']['go-psTsep'].shape[2]
    nDend = 0
    for anm in anmsToUse['dendrites']:
        nDend+=allAnmParams[anm][dffKey]['dsTsep']['dendrites']['go-psTsep'][:,keepRoisDS['dendrites'][anm],:].shape[1]

    somaEmpty = np.empty([4, somaTime, nSoma, nBoot])
    dendEmpty = np.empty([4, dendTime, nDend, nBoot])
    empties = [dendEmpty, somaEmpty]
    for phase in phases:
        for mod in mods:
            for ak, anat in enumerate(anatKeys):
                megaTsepsS[phase][mod][anat] = copy.deepcopy(empties[ak])
                megaTsepsM[phase][mod][anat] = copy.deepcopy(empties[ak])
                megaTsepsO[phase][mod][anat] = copy.deepcopy(empties[ak])
                megaTsepsHH[f'{phase}1'][mod][anat] = copy.deepcopy(empties[ak])
                megaTsepsHH[f'{phase}2'][mod][anat] = copy.deepcopy(empties[ak])
                for hh in megaTsepsHHByKey.values():
                    hh[f'{phase}1'][mod][anat] = copy.deepcopy(empties[ak])
                    hh[f'{phase}2'][mod][anat] = copy.deepcopy(empties[ak])
        
    print(f"Starting bootstrapping {nBoot} iterations...")
    for i in tqdm(range(nBoot)):
        for ak, anat in enumerate(anatKeys):
            if shuffAnms:
                anmIdxs = np.random.choice(np.arange(len(anmsToUse[anat])), len(anmsToUse[anat]))
            else:
                anmIdxs = np.arange(len(anmsToUse[anat]))
            nRois = np.sum(np.hstack([np.sum(keepRoisDS[anat][anm]) for anm in anmsToUse[anat]]))
            for m, mod in enumerate(mods):
                if 'err' in mod:
                    returnErrors = True
                else:
                    returnErrors = False
                mT={}
                mTHH1={}
                mTHH2={}
                for cdKey in label:
                    mT[cdKey] = gather_anmShuff_megaTsepFullSMO(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, cdKey, dffKey, alignTimes[m], anmsToUse, nRois,
                                                    anmIdxs, preComp, enforceMinTrial = enforceMinTrial,
                                                    binary = binary, returnErrors = returnErrors, anmsMinTrial = anmsMinTrial)
                    if preCompHH is not None:
                        # same anmIdxs shuffle so animal resampling matches the FullSMO path
                        mTHH1[cdKey], mTHH2[cdKey] = gather_anmShuff_megaTsepHHSMO(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, cdKey, dffKey, alignTimes[m], anmsToUse, nRois,
                                                    anmIdxs, preCompHH, enforceMinTrial = enforceMinTrial,
                                                    binary = binary, returnErrors = returnErrors, anmsMinTrial = anmsMinTrial)
                for phase in phases:
                    megaTsepsS[phase][mod][anat][:,:,:,i] = mT['s'][phase]
                    megaTsepsM[phase][mod][anat][:,:,:,i] = mT['m'][phase]
                    megaTsepsO[phase][mod][anat][:,:,:,i] = mT['o'][phase]
                    if preCompHH is not None:
                        for cdKey in label:
                            hh = megaTsepsHHByKey[cdKey]
                            hh[f'{phase}1'][mod][anat][:,:,:,i] = mTHH1[cdKey][phase]
                            hh[f'{phase}2'][mod][anat][:,:,:,i] = mTHH2[cdKey][phase]

    print("Returning SMO results...")
    return megaTsepsS, megaTsepsM, megaTsepsO, megaTsepsHH, megaTsepsHHS, megaTsepsHHM, megaTsepsHHO

def plot_CD_projections(nRows,nCols,mod,colors,labels,niceAnatLabels,plot_scaleBar, plotTtypes, plotTtypeLabels, dataSource, normMethod, colScaling, normFactors, ratio, CR_idx, CA_idx, AP_idx, CL_idx, CRCA_HHscroll, APCA_HHscroll, CRCL_HHscroll, anatKeys, frs, factors, window, xspace, xlim, error, CI, nBoots, sharex, sharey, gridspec_kw, vscalar, hscalar):
    nGroups = len(plotTtypes)
    delCol = 0
    if not nCols == nGroups:
        delCol=nCols-nGroups
    figSize = [nCols*hscalar,nRows*vscalar]
    fig,ax=clean_subplots(nRows,nCols,figsize=figSize,sharex=sharex,sharey=sharey,gridspec_kw=gridspec_kw,verbose=True)
    plotData = {}
    ylim = [0,0]
    for col in range(nGroups):
        plotData[col] = {}
        for ak, anat in enumerate(anatKeys):
            row=ak
            F = frs[ak]/factors[ak]
            plotData[col][anat] = {}
            winI = [int(F*window[0]),int(F*window[1])]
            if dataSource[col] == 'CRCA_HHscroll':
                data = copy.deepcopy(CRCA_HHscroll[anat]['post1'][mod]['projs']['post2']) #boot, trial, time
                if normMethod[col] == 'CR+CA':
                    denom = (np.nanmean(abs(data[:,CR_idx,winI[0]:winI[1]]), axis = 1)+np.nanmean(abs(data[:,CA_idx,winI[0]:winI[1]]), axis = 1))
                elif normMethod[col] == 'none':
                    denom = np.ones((nBoots))
                else:
                    raise Exception("missing normMethod")
            elif dataSource[col] == 'APCA_HHscroll':
                data = copy.deepcopy(APCA_HHscroll[anat]['post1'][mod]['projs']['post2']) #boot, trial, time
                if normMethod[col] == 'AP+CA':
                    denom = (np.nanmean(abs(data[:,AP_idx,winI[0]:winI[1]]), axis = 1)+np.nanmean(abs(data[:,CA_idx,winI[0]:winI[1]]), axis = 1))
                elif normMethod[col] == 'none':
                    denom = np.ones((nBoots))
                else:
                    raise Exception("missing normMethod")
            elif dataSource[col] == 'CRCL_HHscroll':
                data = copy.deepcopy(CRCL_HHscroll[anat]['post1'][mod]['projs']['post2']) #boot, trial, time
                if normMethod[col] == 'CR+CL':
                    denom = (np.nanmean(abs(data[:,CR_idx,winI[0]:winI[1]]), axis = 1)+np.nanmean(abs(data[:,CL_idx,winI[0]:winI[1]]), axis = 1))
                elif normMethod[col] == 'none':
                    denom = np.ones((nBoots))
                else:
                    raise Exception("missing normMethod")
            else:
                raise Exception("missing dataSource")
            for t, ttype in enumerate(plotTtypes[col]):
                trialTraces = data[:,ttype,:]/np.expand_dims(denom, axis = 1)
                # mean = np.nanmean(data[:,ttype,:], axis = 0)
                if colScaling[col]:
                    print("colScaling "+anat +" "+str(plotTtypeLabels[col])+ " x"+str(normFactors[anat][ratio]))
                    trialTraces = trialTraces * normFactors[anat][ratio]
                mean = np.nanmean(trialTraces, axis = 0)
                med = np.nanmedian(trialTraces, axis = 0)
                std = np.nanstd(trialTraces, axis = 0, ddof = 1)
                sem = std / np.sqrt(trialTraces.shape[0])
                low,high = np.nanpercentile(trialTraces, [100-(CI+((100-CI)/2)),(CI+((100-CI)/2))], axis = 0)
                x = np.linspace(xspace[0],xspace[1],len(mean))
                plotIdx = np.arange(np.argmin(np.absolute(x-xlim[0])),np.argmin(np.absolute(x-xlim[1]))+1)
                fr = x[1]-x[0]
                if error == 'CI':
                    mainPlot = med
                    errorLow = low
                    errorHigh = high
                elif error == 'sem':
                    mainPlot = mean
                    errorLow = mean-sem
                    errorHigh = mean+sem
                elif error == 'std':
                    mainPlot = mean
                    errorHigh = mean+std
                    errorLow = mean-std
                ylim[0] = np.nanmin([ylim[0],np.nanmin(errorLow)])
                ylim[1] = np.nanmax([ylim[1],np.nanmax(errorHigh)])
                winData = trialTraces[:,int(np.ceil(F*window[0])):int(np.ceil(F*window[1]))] #np.ceil
                if t == 0:
                    pval = np.sum(winData.flatten()<0)/len(winData.flatten())
                else:
                    pval = np.sum(winData.flatten()>0)/len(winData.flatten())
                pText = pstar(pVal,inclP = False)
                print(anat, plotTtypeLabels[col][t]+"("+labels[ttype]+")", np.round(pval,3))
                plotData[col][anat][t] = {'fr':fr,'x':x,'plotIdx':plotIdx, 'mainPlot':mainPlot,'errorHigh':errorHigh,'errorLow':errorLow,'pval':pval,'color':colors[ttype],'label':plotTtypeLabels[col][t]}
            
    ylim[0] = ylim[0] + np.absolute(ylim[1]-ylim[0])*-0.02
    ylim[1] = ylim[1] + np.absolute(ylim[1]-ylim[0])*0.02
    if manualYlim:
        ylim = copy.deepcopy(ylims)
    for col in range(nGroups):
        for ak, anat in enumerate(anatKeys):
            row = ak
            ax[row,col].axhline(0,lw = markerLW, color = markerColor, linestyle = markerLS, alpha = markerAlpha)
            F = frs[ak]/factors[ak]
            winF = x[int(F*window[0]):int(np.ceil(F*window[1]))]
            if winF[0]<0:
                winF = x[int(F*window[0])+1:int(np.ceil(F*window[1]))+1]
            ax[row,col].fill_between([winF[0],winF[-1]], [ylim[0]]*2, [ylim[1]]*2, color = (0,0,0), alpha = 0.1,linewidth=0,edgecolor='none')
            ax[row,col].axvline(0,lw = markerLW, color = markerColor, linestyle = markerLS, alpha = markerAlpha)
            if col == 0:
                ax[row,col].set_ylabel(f'{niceAnatLabels[ak]}\n{yLabel}',fontsize=fontsize,labelpad=0)
            if mod == 'GO':
                alignTime = 'Go Cue'
                markers = [0,-1.25,-2.45]
            elif mod == 'LL':
                alignTime = 'Last Lick'
                markers = [0]
            else:
                alignTime = 'Contact'
                markers = [0]
            for mark in markers:
                ax[row,col].plot([mark]*2,ylim,ls='-',alpha=0.5,color=(0,0,0),linewidth=0.5)
            if row == nRows-1:
                ax[row,col].set_xlabel(f'Time (s; {alignTime})',fontsize=fontsize,labelpad=0)
            ax[row,col].set_xlim(xlim)
            ax[row,col].set_xticks(xticks)
            ax[row,col].set_yticks(yticks)
            ax[row,col].spines[['top','right']].set_visible(False)
            ax[row,col],plot_scaleBar = add_plot_scaleBar(ax[row,col],plot_scaleBar,(0,0,0),True,vertLabel)
        for ak, anat in enumerate(anatKeys):
            row = ak
            for t, ttype in enumerate(plotTtypes[col]):
                ax[row,col].fill_between(plotData[col][anat][t]['x'][plotData[col][anat][t]['plotIdx']], \
                                     plotData[col][anat][t]['errorLow'][plotData[col][anat][t]['plotIdx']], \
                                     plotData[col][anat][t]['errorHigh'][plotData[col][anat][t]['plotIdx']], \
                                     color = plotData[col][anat][t]['color'], alpha = 0.5,linewidth=0,edgecolor='none')
                if lineStyles[col][t] == '-':
                    ax[row,col].plot(plotData[col][anat][t]['x'][plotData[col][anat][t]['plotIdx']], plotData[col][anat][t]['mainPlot'][plotData[col][anat][t]['plotIdx']],\
                                 color = plotData[col][anat][t]['color'], label = plotData[col][anat][t]['label'],lw = linewidth, ls = lineStyles[col][t])
                else:
                    ax[row,col].plot(plotData[col][anat][t]['x'][plotData[col][anat][t]['plotIdx']], plotData[col][anat][t]['mainPlot'][plotData[col][anat][t]['plotIdx']],\
                                 color = plotData[col][anat][t]['color'], label = plotData[col][anat][t]['label'],lw = linewidth, ls = lineStyles[col][t], dashes = lineDashes[col][t])
            ax[row,col].set_ylim(ylim)
            ax[row,col].tick_params(axis='both', which='major', labelsize=fontsize) #
            ax[row,col].tick_params(axis='x', which='major', labelsize=fontsize,pad= 2) #
            ax[row,col].tick_params(axis='y', which='major', labelsize=fontsize,pad= 0) #
            ax[row,col].legend(frameon=False,reverse=reverseLeg,fontsize=fontsize,loc=loc,handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 
            if clean:
                ax[row,col].spines[['bottom','left']].set_visible(False)
                ax[row,col].set_xticks([])
                ax[row,col].set_yticks([])
                ax[row,col].set_xlabel('')
                ax[row,col].set_ylabel('')
    if delCol>0:
        col = nGroups-1
        for c in range(delCol):
            col+=1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
    if PDF_export_active:
        with PdfPages(os.path.join(figSaveDir,figName)) as pdf:
            pdf.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600) 
    display_clean_subplots(fig,ax)

    ##################################################################################################

##################################################################################################
#Clean Matplotlib Figures
def enforce_equal_axes(fig, axes):
    """
    Force all axes to have identical positions after layout/rendering.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing the axes.
    axes : np.ndarray of matplotlib.axes.Axes
        Array of axes (e.g., returned by clean_subplots/plt.subplots) to
        force into identical positions (intersection of all axes boxes).

    Example
    -------
    importlib.reload(jsm)
    fig, axes = enforce_equal_axes(fig, axes)
    """
    # Ensure layout is fully computed
    fig.canvas.draw()

    # Get positions of all axes
    positions = [ax.get_position() for ax in axes.flat]

    # Use the smallest common box (intersection-like behavior)
    left = max(p.x0 for p in positions)
    bottom = max(p.y0 for p in positions)
    right = min(p.x1 for p in positions)
    top = min(p.y1 for p in positions)

    unified_pos = [left, bottom, right - left, top - bottom]

    # Apply to all axes
    for ax in axes.flat:
        ax.set_position(unified_pos)
    return fig,axes
def enforce_axes(fig, axes, mode="all"):
    """
    Enforce equal axes positions.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axes : array-like of Axes
    mode : str
        "all"  -> all axes identical
        "row"  -> equalize within each row
        "col"  -> equalize within each column

    Example
    -------
    mode = 'all'    # 'all' | 'row' | 'col'
    importlib.reload(jsm)
    fig, axes = enforce_axes(fig, axes, mode=mode)
    """

    import numpy as np

    # Ensure layout is finalized
    fig.canvas.draw()

    axes = np.asarray(axes)

    # Helper: filter valid axes (skip deleted/invisible)
    def valid_axes(ax_list):
        return [ax for ax in ax_list if ax is not None and ax.get_visible()]

    # Helper: compute intersection box
    def intersect_positions(ax_list):
        positions = [ax.get_position() for ax in ax_list]

        left = max(p.x0 for p in positions)
        bottom = max(p.y0 for p in positions)
        right = min(p.x1 for p in positions)
        top = min(p.y1 for p in positions)

        return [left, bottom, right - left, top - bottom]

    # --- Mode: all ---
    if mode == "all":
        ax_list = valid_axes(axes.flat)
        if not ax_list:
            return

        pos = intersect_positions(ax_list)
        for ax in ax_list:
            ax.set_position(pos)

    # --- Mode: row ---
    elif mode == "row":
        for row in axes:
            ax_list = valid_axes(row)
            if not ax_list:
                continue

            pos = intersect_positions(ax_list)
            for ax in ax_list:
                ax.set_position(pos)

    # --- Mode: col ---
    elif mode == "col":
        for col in axes.T:
            ax_list = valid_axes(col)
            if not ax_list:
                continue

            pos = intersect_positions(ax_list)
            for ax in ax_list:
                ax.set_position(pos)

    else:
        raise ValueError("mode must be one of: 'all', 'row', 'col'")

    return fig,axes
def figure_grid(
    nrows,
    ncols,
    figsize,
    gridspec_kw=None,
    sharex=False,
    sharey=False):
    """
    Create a subplot grid with fixed overall figure size and GridSpec margins.

    Parameters
    ----------
    nrows, ncols : int
        Grid shape
    figsize : tuple
        (width, height) in inches
    gridspec_kw : dict
        GridSpec kwargs including margins and spacing:
        {
            "left": ..., "right": ..., "top": ..., "bottom": ...,
            "wspace": ..., "hspace": ...
        }
        Margins are normalized (0–1), same as matplotlib defaults.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    axes : np.ndarray of matplotlib.axes.Axes
        Axes array reshaped to (nrows, ncols).
    info : dict
        Computed figure/axes sizes, margins, and spacing in inches (see
        keys 'figure_size_in', 'axes_size_in', 'margins_in', 'spacing_in').

    Example
    -------
    nrows = 2                             # number of rows of axes
    ncols = 3                             # number of columns of axes
    figsize = (8, 5)                      # (width, height) in inches
    gridspec_kw = {'left': 0.1, 'right': 0.95, 'top': 0.9, 'bottom': 0.1,
                   'wspace': 0.3, 'hspace': 0.3}  # GridSpec margins/spacing (normalized 0-1); None -> matplotlib defaults
    sharex = False                        # True | False
    sharey = False                        # True | False
    importlib.reload(jsm)
    fig, axes, info = figure_grid(
        nrows, ncols, figsize, gridspec_kw=gridspec_kw, sharex=sharex, sharey=sharey)
    """
    fig_width, fig_height = figsize

    if gridspec_kw is None:
        gridspec_kw = {}

    fig = plt.figure(figsize=figsize)

    gs = fig.add_gridspec(nrows, ncols, **gridspec_kw)

    axes = gs.subplots(sharex=sharex, sharey=sharey)
    # axes = np.atleast_2d(axes)
    axes = np.array(axes, dtype=object).reshape(nrows, ncols)


    # --- Extract normalized margins (with defaults) ---
    left = gridspec_kw.get("left", 0.125)
    right = gridspec_kw.get("right", 0.9)
    bottom = gridspec_kw.get("bottom", 0.11)
    top = gridspec_kw.get("top", 0.88)
    wspace = gridspec_kw.get("wspace", 0.2)
    hspace = gridspec_kw.get("hspace", 0.2)

    precision = 3

    # --- Convert margins to inches ---
    margins_in = {
        "left": left * fig_width,
        "right": (1 - right) * fig_width,
        "bottom": bottom * fig_height,
        "top": (1 - top) * fig_height,
    }

    # --- Available plotting area ---
    usable_width = (right - left) * fig_width
    usable_height = (top - bottom) * fig_height

    # --- Axes size calculation ---
    # Matplotlib spacing is fraction of average axis size
    ax_width = usable_width / (ncols + (ncols - 1) * wspace)
    ax_height = usable_height / (nrows + (nrows - 1) * hspace)

    # --- Spacing in inches ---
    spacing_in = {
        "wspace": wspace * ax_width,
        "hspace": hspace * ax_height
    }

    info = {
        "figure_size_in": {
            "width": fig_width,
            "height": fig_height
        },
        "axes_size_in": {
            "width": ax_width,
            "height": ax_height
        },
        "margins_in": margins_in,
        "spacing_in": spacing_in
    }


    return fig, axes, info
def clean_subplots(nrows,ncols,figsize = False, facecolor='none', constrained_layout = True, tight_layout = False, \
    sharey = False, sharex = False, gridspec_kw = None, verbose = False, units = 'in'):
    """
    Create a subplot grid with a transparent/clean background, to be used in
    place of plt.subplots() throughout the project for consistent figure styling.

    Parameters
    ----------
    nrows : int
        Number of rows of axes.
    ncols : int
        Number of columns of axes.
    figsize : tuple or False
        (width, height) in `units`. False -> let matplotlib choose the
        default figure size.
    facecolor : str or tuple
        Face color for the figure and each axes, e.g. 'none' for transparent
        or an RGB tuple such as (0.2, 0.4, 0.8).
    constrained_layout : bool
        If True (and gridspec_kw is None), use matplotlib's constrained
        layout engine.
    tight_layout : bool
        Passed through to plt.subplots() as tight_layout.
    sharey : bool
        Share the y-axis across subplots.
    sharex : bool
        Share the x-axis across subplots.
    gridspec_kw : dict or None
        GridSpec kwargs (margins/spacing) forwarded to figure_grid(); when
        provided, figsize must also be given. None -> use plt.subplots()
        directly instead of figure_grid().
    verbose : bool
        If True (and gridspec_kw is given), print the computed figure/axes
        sizing info from figure_grid().
    units : str
        Units for `figsize`: 'in', 'cm', or 'mm'.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    ax : matplotlib.axes.Axes or np.ndarray of Axes
        The created axes.

    Example
    -------
    nrows = 2                  # number of rows of axes
    ncols = 2                  # number of columns of axes
    figsize = (8, 6)           # (width, height) in `units`; False -> matplotlib default size
    facecolor = 'none'         # 'none' | any color string | RGB tuple e.g. (0.2, 0.4, 0.8)
    constrained_layout = True  # True | False   (only applied when gridspec_kw is None)
    tight_layout = False       # True | False
    sharey = False             # True | False
    sharex = False             # True | False
    gridspec_kw = None         # None -> use plt.subplots(); or dict of GridSpec margins/spacing
    verbose = False            # True | False   (prints figure_grid sizing info)
    units = 'in'               # 'in' | 'cm' | 'mm'
    importlib.reload(jsm)
    fig, ax = clean_subplots(
        nrows, ncols, figsize=figsize, facecolor=facecolor, constrained_layout=constrained_layout,
        tight_layout=tight_layout, sharey=sharey, sharex=sharex, gridspec_kw=gridspec_kw,
        verbose=verbose, units=units)
    """
    _to_in = {'in': 1.0, 'cm': 1/2.54, 'mm': 1/25.4}
    if units not in _to_in:
        raise ValueError(f"units must be 'in', 'cm', or 'mm', got '{units}'")
    if figsize:
        scale = _to_in[units]
        figsize = (figsize[0] * scale, figsize[1] * scale)

    if gridspec_kw is None:
        if figsize:
            fig,ax = plt.subplots(nrows,ncols,figsize=figsize,tight_layout=tight_layout,sharey=sharey,sharex=sharex)
        else:
            fig,ax = plt.subplots(nrows,ncols,tight_layout=tight_layout,sharey=sharey,sharex=sharex)
    else:
        if figsize:
            fig,ax,info = figure_grid(nrows,ncols,figsize,gridspec_kw,sharex=sharex,sharey=sharey)
        else:
            raise Exception("gridspec_kw is not None, but figsize is False, this is not currently supported")


    if hasattr(ax, 'shape'):
        if len(ax.shape)>1:
            for row in range(ax.shape[0]):
                for col in range(ax.shape[1]):
                    ax[row,col].set_facecolor(facecolor)
                    if facecolor=='none':
                        ax[row,col].patch.set_alpha(0.0) # Axes background
        else:
            for row in range(ax.shape[0]):
                ax[row].set_facecolor(facecolor)
                if facecolor=='none':
                    ax[row].patch.set_alpha(0.0) # Axes background
    else:
        ax.set_facecolor(facecolor)
        if facecolor=='none':
            ax.patch.set_alpha(0.0) # Axes background
    if constrained_layout and gridspec_kw is None:
        # fig.set_constrained_layout(True)
        fig.set_layout_engine('constrained')
    fig.set_facecolor(facecolor)
    if facecolor=='none':
        fig.patch.set_alpha(0.0) # Figure background
    return fig,ax
def noclip(ax):
    """
    Turn off all clipping in axes ax; call immediately before drawing/showing.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes whose artists (collections, patches, lines, texts) should
        have clipping disabled.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The same axes, with clipping disabled.

    Example
    -------
    importlib.reload(jsm)
    ax = noclip(ax)
    """
    ax.set_clip_on(False)
    artists = []
    artists.extend(ax.collections)
    artists.extend(ax.patches)
    artists.extend(ax.lines)
    artists.extend(ax.texts)
    artists.extend(ax.artists)
    for a in artists:
        a.set_clip_on(False)
    return ax
def remove_all_clipping(fig,ax):
    """
    Disable clipping on every axes in a figure (single axes or an array of
    axes from clean_subplots/plt.subplots).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing the axes (unused for the disabling logic
        itself, but returned alongside ax for chaining).
    ax : matplotlib.axes.Axes or np.ndarray of Axes
        Single axes, or 1D/2D array of axes, to disable clipping on.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes or np.ndarray of Axes
        The same ax(es), with clipping disabled.

    Example
    -------
    importlib.reload(jsm)
    fig, ax = remove_all_clipping(fig, ax)
    """
    if hasattr(ax, 'shape'):
        if len(ax.shape)>1:
            for row in range(ax.shape[0]):
                for col in range(ax.shape[1]):
                    ax[row,col] = noclip(ax[row,col])
        else:
            for a in ax:
                a = noclip(a)
    else:
        ax.set_facecolor((1,1,1))
        ax.patch.set_alpha(1) # Axes background
        ax = noclip(ax)
    return fig,ax
def display_clean_subplots(fig,ax):
    """
    Temporarily restore an opaque white background on a clean_subplots()
    figure/axes and show() it, for previewing before final (transparent)
    export.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to make opaque and display.
    ax : matplotlib.axes.Axes or np.ndarray of Axes
        Single axes, or 1D/2D array of axes, to make opaque.

    Example
    -------
    importlib.reload(jsm)
    display_clean_subplots(fig, ax)
    """
    if hasattr(ax, 'shape'):
        if len(ax.shape)>1:
            for row in range(ax.shape[0]):
                for col in range(ax.shape[1]):
                    ax[row,col].set_facecolor((1,1,1))
                    ax[row,col].patch.set_alpha(1) # Axes background
        else:
            for a in ax:
                a.set_facecolor((1,1,1))
                a.patch.set_alpha(1) # Axes background
    else:
        ax.set_facecolor((1,1,1))
        ax.patch.set_alpha(1) # Axes background
    fig.set_facecolor((1,1,1))
    fig.patch.set_alpha(1) # Figure background
    plt.show()
def set_dynamic_suptitle(fig, text, margin = 0.02, **kwargs):
    """
    Add a suptitle above all existing graphics in the figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure object.
    text : str
        The suptitle text.
    **kwargs :
        Additional keyword arguments passed to fig.suptitle().
    """
    # Draw once to get accurate positions
    fig.canvas.draw()

    # Find the topmost extent of all axes elements
    max_y = 0
    for ax in fig.get_axes():
        # Get the bounding box of the axis including labels/titles
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        max_y = max(max_y, bbox.y1 / fig.bbox.y1)

    # Add a small margin above the highest element
    
    y_position = min(1.0, max_y + margin)

    # Place the suptitle
    return fig.suptitle(text, y=y_position, **kwargs)
def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
##################################################################################################
##################################################################################################    
def summary_stats(tempData,nBins=100,binWidth=0,temp_minVal=np.nan,temp_maxVal=np.nan,\
    calcHist=True,verbose=False, allow_overages = False, splitZeroBin = False, \
    calc_hist_log10 = True, calc_cumulative_hist = True, calc_norm_hist = True, calc_cumulative_freq = True, calc_pdf = True, 
    save_input = False, force_dType=False, dType='float32', exclude_NaNs = True, ddof = 0, ci_level = 95, warningsOn = True,simple_verbose=False):
    """
    Calculate summary statistics for the provided data.

    Parameters:
    tempData: input data to analyze
    nBins: number of bins for histogram
    binWidth: width of each bin for histogram
    temp_minVal: minimum value for histogram
    temp_maxVal: maximum value for histogram
    calcHist: if True, calculate histogram and related statistics
    verbose: if True, print additional information
    allow_overages: if True, allow values greater than the histogram limits
    splitZeroBin: if True, split the zero bin in the histogram
    calc_hist_log10: if True, calculate logarithmic histogram
    calc_cumulative_hist: if True, calculate cumulative histogram
    calc_norm_hist: if True, calculate normalized histogram
    calc_cumulative_freq: if True, calculate cumulative frequency
    calc_pdf: if True, calculate probability density function
    save_input: if True, save the input data in the stats dictionary
    force_dType: if True, force the data to a specific data type
    dType: data type to force the data to
    exclude_NaNs: if True, exclude NaN values from the calculations
    ci_level: confidence level (percent) for the percentile-based confidence
        interval reported as 'ci_lower' and 'ci_upper'. Defaults to 95.
    warningsOn: if True, print warnings

    Returns:
    stats: a dictionary containing the calculated statistics
    """

    stats={}
    tempData = np.array(tempData)
    tempData=tempData.flatten()
    if len(tempData) == 0:
        if warningsOn:
            print("<<WARNING>> No data to analyze")
        stats['success'] = False
        stats['n'] = np.nan
        stats['mean'] = np.nan
        stats['gmean'] = np.nan
        stats['max'] = np.nan
        stats['min'] = np.nan
        stats['med'] = np.nan
        stats['std'] = np.nan
        stats['sem'] = np.nan
        stats['cv'] = np.nan
        stats['sum'] = np.nan
        stats['mode']=np.nan
        stats['q1']=np.nan
        stats['q3']=np.nan
        stats['ci_lower']=np.nan
        stats['ci_upper']=np.nan
        stats['iqr']=np.nan
        stats['np_skew']=np.nan
        stats['fmc_skew']=np.nan
    else:
        if force_dType:
            try:
                tempData = np.array(tempData).astype(dType)
            except:
                pass
        ##################
        #Basic Stats
        if save_input:
            stats['data'] = copy.deepcopy(tempData)
        stats['success'] = False
        if exclude_NaNs:
            stats['n']=sum(~np.isnan(tempData))
            stats['mean']=np.nanmean(tempData)
            stats['max']=np.nanmax(tempData)
            stats['min']=np.nanmin(tempData)
            try:
                stats['minPos']=np.nanmin(tempData[tempData>0])
            except:
                stats['minPos']=np.nan
            stats['med']=np.nanmedian(tempData,axis=0)
            stats['std']=np.nanstd(tempData,ddof=ddof)
            stats['sum']=np.nansum(tempData)
            try:
                stats['q1']=np.nanpercentile(tempData,25)
                stats['q3']=np.nanpercentile(tempData,75)
            except:
                stats['q1']=np.nan
                stats['q3']=np.nan
            try:
                stats['ci_lower']=np.nanpercentile(tempData,(100-ci_level)/2)
                stats['ci_upper']=np.nanpercentile(tempData,100-(100-ci_level)/2)
            except:
                stats['ci_lower']=np.nan
                stats['ci_upper']=np.nan
            # stats['mode']=mode(tempData, nan_policy='omit')
            stats['gmean']=gmean(tempData[np.isfinite(tempData)])

        else:
            stats['n']=sum(tempData)
            stats['mean']=np.mean(tempData)
            stats['max']=np.max(tempData)
            stats['min']=np.min(tempData)
            try:
                stats['minPos']=np.min(tempData[tempData>0])
            except:
                stats['minPos']=np.nan
            stats['med']=np.median(tempData,axis=0)
            stats['std']=np.std(tempData,ddof=ddof)
            stats['sum']=np.sum(tempData)
            stats['mode']=mode(tempData)
            try:
                stats['q1']=np.percentile(tempData,25)
                stats['q3']=np.percentile(tempData,75)
            except:
                stats['q1']=np.nan
                stats['q3']=np.nan
            try:
                stats['ci_lower']=np.percentile(tempData,(100-ci_level)/2)
                stats['ci_upper']=np.percentile(tempData,100-(100-ci_level)/2)
            except:
                stats['ci_lower']=np.nan
                stats['ci_upper']=np.nan
            # stats['mode']=mode(tempData, nan_policy='propogate')
            stats['gmean']=gmean(tempData)
        try:
            stats['iqr']=stats['q3']-stats['q1']
        except:
            stats['iqr']=np.nan
        try:
            stats['sem']=stats['std']/np.sqrt(stats['n'])
        except:
            stats['sem']=np.nan
        try:
            stats['cv']=stats['std']/stats['mean']
        except:
            stats['cv']=np.nan
        try:
            stats['np_skew']=(stats['mean']-stats['med'])/stats['std'] # nonparametric skew
        except:
            stats['np_skew']=np.nan
        try:
            if exclude_NaNs:
                stats['fmc_skew']=float(skew(tempData, axis=0, bias=True, nan_policy='omit')) #Fisher's moment coefficient of skewness
            else:
                stats['fmc_skew']=float(skew(tempData, axis=0, bias=True, nan_policy='propogate')) #Fisher's moment coefficient of skewness

        except:
            stats['fmc_skew'] = np.nan
        ##################
        if calcHist:
            stats['binWidth']=binWidth
            stats['minVal']=temp_minVal
            stats['maxVal']=temp_maxVal
            stats['nBins']=nBins
            try:
                # try:
                if nBins and ('float' in str(type(stats['minVal'])) or 'int' in str(type(stats['minVal']))) and \
                    ('float' in str(type(stats['maxVal'])) or 'int' in str(type(stats['maxVal']))):
                    binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                elif binWidth>0 and ('float' in str(type(stats['minVal'])) or 'int' in str(type(stats['minVal']))) and \
                    ('float' in str(type(stats['maxVal'])) or 'int' in str(type(stats['maxVal']))):
                    binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                elif nBins:
                    if exclude_NaNs:
                        stats['minVal'] = np.floor(np.nanmin(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.nanmax(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                    else:
                        stats['minVal'] = np.floor(np.min(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.max(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.linspace(stats['minVal'], stats['maxVal'], nBins)
                elif binWidth>0:
                    if exclude_NaNs:
                        stats['minVal'] = np.floor(np.nanmin(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.nanmax(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                    else:
                        stats['minVal'] = np.floor(np.min(tempData)*(1/binWidth))/(1/binWidth)
                        stats['maxVal'] = np.ceil(np.max(tempData)*(1/binWidth))/(1/binWidth)
                        binEdges=np.arange(stats['minVal'], stats['maxVal'], binWidth)
                binCenters=binEdges[0:len(binEdges)-1]+(binEdges[1]-binEdges[0])/2
                zeroMin = False
                if exclude_NaNs:
                    if np.nanmin(tempData) == 0 and splitZeroBin:
                        zeroMin = True
                else:
                    if np.min(tempData) == 0 and splitZeroBin:
                        zeroMin = True
                if zeroMin:
                    if verbose:
                        print("NOTE: adding true zero bin")
                    binEdges[0] = 0.0000001
                    binEdges = np.insert(binEdges,0,0)
                    binCenters = np.insert(binCenters,0,0)

                hist,BinCenters=np.histogram(tempData,bins=binEdges)
                # except:
                #     hist,BinCenters=np.histogram(tempData)
                if allow_overages and np.any(tempData>binEdges[-1]):
                    stats['overage'] = True
                    if verbose:
                        print("NOTE: adding an extra bin for values greater than the limits")
                    overage = np.sum(tempData>binEdges[-1])
                    hist = np.append(hist,overage)
                    binCenters = np.append(binCenters,binCenters[-1]+(binCenters[-1]-binCenters[-3]))
                    if exclude_NaNs:
                        binEdges = np.append(binEdges,np.nanmax(tempData))
                    else:
                        binEdges = np.append(binEdges,np.max(tempData))

                hist_log10 = copy.deepcopy(hist)
                # hist_log10[hist_log10<1] = 0.1
                hist_log10 = np.log10(hist_log10)
                hist_log10[np.isinf(hist_log10)] = 0
                norm_hist=hist/np.nanmax(hist)
                cumulative_hist=np.cumsum(hist)
                pdf = hist/stats['n']
                if np.nanmax(cumulative_hist) == np.nan:
                    cumulative_freq=cumulative_hist
                else:
                    cumulative_freq=cumulative_hist/np.nanmax(cumulative_hist)
                stats['binEdges']=binEdges.astype('float32')
                stats['binCenters']=binCenters.astype('float32')
                stats['hist']=hist.astype('int32')
                if calc_hist_log10:
                    stats['hist_log10']=hist_log10.astype('float32')
                if calc_cumulative_hist:
                    stats['cumulative_hist']=cumulative_hist.astype('float32')
                if calc_norm_hist:
                    stats['norm_hist']=norm_hist.astype('float32')
                if calc_cumulative_freq:
                    stats['cumulative_freq']=cumulative_freq.astype('float32')
                if calc_pdf:
                    stats['pdf']=pdf.astype('float32')
                stats['success'] = True
                if not stats['nBins']:
                    stats['nBins'] = len(hist)
            except:
                if warningsOn:
                    print("<<WARNING>> Unable to calculate histograms...")
                # stats['hist'] = [np.nan]
                # stats['cumulative_hist'] = [np.nan]
                # stats['norm_hist'] = [np.nan]
                # stats['cumulative_freq'] = [np.nan]
        else:
            if verbose:
                print("Skipping Histograms")
    ##################
    if verbose:
        if not simple_verbose:
            print("n    = "+str(stats['n']))
            print("mean = "+str(stats['mean']))
            print("max  = "+str(stats['max']))
            print("min  = "+str(stats['min']))
            print("med  = "+str(stats['med']))
            print("std  = "+str(stats['std']))
            print("sem  = "+str(stats['sem']))
            print("cv   = "+str(stats['cv']))
            print("sum  = "+str(stats['sum']))
            print("ci_lower ("+str(ci_level)+"%) = "+str(stats.get('ci_lower', np.nan)))
            print("ci_upper ("+str(ci_level)+"%) = "+str(stats.get('ci_upper', np.nan)))
        # Formatted summary lines (significant-digit rounded)
        def _sig(x, sig=3):
            if x is None or not np.isfinite(x):
                return str(x)
            return str(float(f"{x:.{sig}g}"))
        print("mean ± SEM:   "+_sig(stats['mean'])+" ± "+_sig(stats['sem'])+
              " (n = "+str(stats['n'])+")")
        print("median ± CI:   "+_sig(stats['med'])+
              " + "+_sig(stats.get('ci_upper', np.nan)-stats['med'])+
              " - "+_sig(stats['med']-stats.get('ci_lower', np.nan))+
              " (n = "+str(stats['n'])+")")
    return stats
def pstar(pval,inclP=True,simpleStar=True):
    pstr = str(np.round(pval,4))
    if simpleStar:
        if inclP:
            if pval < 0.05:
                return '* (p = '+pstr+')'
            else:
                return 'n.s. (p = '+pstr+')'
        else:
            if pval < 0.05:
                return '*'
            else:
                return 'n.s.'
    else:
        if inclP:
            if pstr == '0.0' and pval < 0.0001:
                return '**** (p < 0.0001)'
            elif pval < 0.0001:
                return '**** (p = '+pstr+')'
            elif pval < 0.001:
                return '*** (p = '+pstr+')'
            elif pval < 0.01:
                return '** (p = '+pstr+')'
            elif pval < 0.05:
                return '* (p = '+pstr+')'
            else:
                return 'n.s. (p = '+pstr+')'
        else:
            if pstr == '0.0' and pval < 0.0001:
                return '****'
            elif pval < 0.0001:
                return '****'
            elif pval < 0.001:
                return '***'
            elif pval < 0.01:
                return '**'
            elif pval < 0.05:
                return '*'
            else:
                return 'n.s.'
def multi_test_correction(pVals,method='holm',alpha = 0.05,simpleStar = True, verbose=True):
    pVals = np.hstack(pVals)
    multCorr = {}
    multCorr['fail'],multCorr['pvals_corrected'],multCorr['alphacSidak'],multCorr['alphacBonf'] = \
        multipletests(pVals, method=method, alpha = alpha)
    if verbose:
        print(f'MultCompar Correction method = {method} alpha = {alpha}')
        print(f'pVals      = {pVals}')
        print(f'pVals_Corr = {multCorr["pvals_corrected"]}')
    pText = []
    for p in range(len(pVals)):
        if multCorr['fail'][p]:
            pText.append(pstar(multCorr['pvals_corrected'][p],inclP=False, simpleStar = simpleStar))
        else:
            pText.append('n.s.')

    return multCorr, pText    
def smart_floor(val):
    precision = precision_lookup(val)
    val = np.floor(val * 10**precision) / 10**precision
    if np.absolute(val) == 0:
        val = 0
    return val
def smart_ceil(val):
    precision = precision_lookup(val)
    val = np.ceil(val * 10**precision) / 10**precision
    if np.absolute(val) == 0:
        val = 0
    return val
def smart_round(val):
    precision = precision_lookup(val) + 1
    val = np.round(val * 10**precision) / 10**precision
    if np.absolute(val) == 0:
        val = 0
    return val
def smart_float2str(val,decimals=1):
    val = str(np.round(val,decimals = decimals))
    if decimals > 0:
        if not '.' in val:
            val = val+'.'
        while len(val.split('.')[-1]) < decimals:
            val = val+'0'
    return val
def precision_lookup(val):
    if np.absolute(val) < 0.00001:
        precision = 6
    elif np.absolute(val) < 0.0001:
        precision = 5
    elif np.absolute(val) < 0.001:
        precision = 4
    elif np.absolute(val) < 0.01:
        precision = 3
    elif np.absolute(val) < 0.1:
        precision = 2
    elif np.absolute(val) < 1:
        precision = 1
    elif np.absolute(val) < 10:
        precision = 0
    elif np.absolute(val) < 100:
        precision = -1
    elif np.absolute(val) < 1000:
        precision = -2
    elif np.absolute(val) < 10000:
        precision = -3
    elif np.absolute(val) < 100000:
        precision = -4
    elif np.absolute(val) < 1000000:
        precision = -5
    else:
        precision = -6
    return precision
##################################################################################################

##################################################################################################
#RGB and colormap functions
##################################################################################################
def save_rgb_image(filepath,image):
    minval = np.min(image)
    maxval = np.max(image)
    image = (image - minval)/ (maxval - minval + 0.01)
    image *= 255
    image = image.astype(np.uint8)
    tifffile.imwrite(filepath, image) 
##################################################################################################
def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray
##################################################################################################
def create_smooth_diverging_cmap(cmap1="viridis", cmap2="inferno", n=256, cutoff=0.8, center_width=0.3):
    """
    Creates a custom map with a SMOOTH fade to white in the center.
    
    Parameters:
    - cutoff: (0.8) stops the original maps before they get too yellow.
    - center_width: (0.03) The fraction of the map dedicated to the FADE to white.
    """
    # 1. Setup Indices
    # We split the map into 4 zones: Left Body, Left Fade, Right Fade, Right Body
    center_idx = n // 2
    fade_len = int(n * center_width / 2)  # Number of pixels for the fade on EACH side
    body_len = center_idx - fade_len      # Remaining pixels for the main colormaps

    # 2. Get Key Colors (The anchors for our interpolation)
    # We need the specific RGBA values at the cutoff points to start the fade
    c1_map = plt.colormaps[cmap1]
    c2_map = plt.colormaps[cmap2]
    
    # Anchor colors
    left_end_color = c1_map(cutoff)   # The Teal color where Viridis stops
    right_start_color = c2_map(cutoff) # The Orange color where Plasma starts
    white = np.array([1.0, 1.0, 1.0, 1.0]) # Pure White

    # 3. Generate the Main Bodies
    # Sample from 0 to cutoff
    left_body = c1_map(np.linspace(0, cutoff, body_len))
    right_body = c2_map(np.linspace(cutoff, 0, body_len))

    # 4. Generate the Smooth Bridges (The Fix)
    # Linearly interpolate from Left_Anchor -> White -> Right_Anchor
    left_bridge = np.linspace(left_end_color, white, fade_len)
    right_bridge = np.linspace(white, right_start_color, fade_len)

    # 5. Stack Everything: Left Body -> Fade In -> Fade Out -> Right Body
    combined_colors = np.vstack((left_body, left_bridge, right_bridge, right_body))

    return mcolors.LinearSegmentedColormap.from_list(
        f"{cmap1}_smooth_white_{cmap2}", 
        combined_colors
    )
##################################################################################################
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """
    Function to truncate a colormap.

    Args:
        cmap (str or mcolors.Colormap): Colormap to be altered.
        minval (float, optional): Minimum value (0.0 to 1.0) of the original 
                                   colormap to include. Defaults to 0.0.
        maxval (float, optional): Maximum value (0.0 to 1.0) of the original 
                                   colormap to include. Defaults to 1.0.
        n (int, optional): Number of intervals in the new colormap gradient. 
                           Defaults to 100.

    Returns:
        mcolors.LinearSegmentedColormap: The new, truncated colormap instance.
    """
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    
    # Create the new colormap by sampling the original
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap
##################################################################################################
def generate_cmap(color,ncolors,cmap_name='cmap'):
    tempRGB = []
    tempBGR = []
    if 'tuple' in str(type(color)):
        RGB=list(color)
        RGB = tuple(RGB[0:3])
        BRG=[RGB[2],RGB[1],RGB[0]]
        
        colors = [(0,0,0),RGB]
        cmap1 = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=ncolors)
        try:
            cmap = colormaps.get_cmap(cmap1)
        except:
            cmap = plt.get_cmap(cmap1).copy()
        tempRGB=[(0,0,0),RGB]
        tempBGR=[(0,0,0),BRG]
    elif 'list' in str(type(color)) and len(color)==2:
        tempRGB = list()
        tempBGR = list()
        for c in range(len(color)):
            tempRGB.append(list(color[c]))
            tempBGR.append([tempRGB[c][2],tempRGB[c][1],tempRGB[c][0]])

        cmap1 = mcolors.LinearSegmentedColormap.from_list(cmap_name, color, N=ncolors)
        try:
            cmap = colormaps.get_cmap(cmap1)
        except:
            cmap = plt.get_cmap(cmap1).copy()
    elif 'str' in str(type(color)):
        if 'custom' in color or 'Custom' in color:
            if 'OWB' in color:
                colors = [(1,0.5,0), (1,1,1), (0,0,1)]
            elif 'OKB' in color:
                colors = [(1,0.5,0), (0,0,0), (0,0,1)]
            elif 'BWO' in color:
                colors = [(0,0,1), (1,1,1), (1,0.5,0)]
            elif 'BKO' in color:
                colors = [(0,0,1), (0,0,0), (1,0.5,0)]
            elif 'MWG' in color:
                colors = [(1,0,1), (1,1,1), (0,1,0)]
            elif 'MKG' in color:
                colors = [(1,0,1), (0,0,0), (0,1,0)]
            elif 'GWM' in color:
                colors = [(0,1,0), (1,1,1), (1,0,1)]
            elif 'GKM' in color:
                colors = [(0,1,0), (0,0,0), (1,0,1)]
            elif 'GWB' in color:
                colors = [(0,1,0), (1,1,1), (0,0,1)]
            elif 'GKB' in color:
                colors = [(0,1,0), (0,0,0), (0,0,1)]
            elif 'BWG' in color:
                colors = [(0,0,1), (1,1,1), (0,1,0)]
            elif 'BKG' in color:
                colors = [(0,0,1), (0,0,0), (0,1,0)]
            elif 'RWB' in color:
                colors = [(1,0,0), (1,1,1), (0,0,1)]
            elif 'RKB' in color:
                colors = [(1,0,0), (0,0,0), (0,0,1)]
            elif 'BWR' in color:
                colors = [(0,0,1), (1,1,1), (1,0,0)]
            elif 'BKR' in color:
                colors = [(0,0,1), (0,0,0), (1,0,0)]

            cmap1 = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=ncolors)
            try:
                cmap = colormaps.get_cmap(cmap1)
            except:
                cmap = plt.get_cmap(cmap1).copy()


        elif len(color)==1:
            tempRGB,tempBGR=ConvertColorLetter(color)
            colors = [(0, (0, 0, 0)),(1, tempRGB)]
            cmap1 = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=ncolors)
            try:
                cmap = colormaps.get_cmap(cmap1)
            except:
                cmap = plt.get_cmap(cmap1).copy()
        else:
            try:
                cmap = plt.get_cmap(color, ncolors)
            except:
                try:
                    cmap = colormaps.get_cmap(cmap1)
                except:
                    cmap = plt.get_cmap(cmap1).copy()
        tempRGB=[]
        tempBGR=[]
        for c in range(ncolors):
            tempRGB.append([cmap(c)[0],cmap(c)[1],cmap(c)[2]])
            tempBGR.append([cmap(c)[2],cmap(c)[1],cmap(c)[0]])
    elif 'dict' in str(type(color)):
        if color['type'] == 'smooth_diverging':
            cmap = create_smooth_diverging_cmap(cmap1=color['cmap1'],cmap2= color['cmap2'], cutoff=color['cutoff'], center_width=color['center_width'], n=ncolors) # Increased width slightly for visibility
            cmap = truncate_colormap(cmap, color['minVal'], color['maxVal'] , n=ncolors)

            tempRGB = cmap(np.linspace(0, 1, ncolors))

            # tempRGB=[]
            # tempBGR=[]
            # for c in range(ncolors):
            #     tempRGB.append([cmap(c)[0],cmap(c)[1],cmap(c)[2]])
            #     tempBGR.append([cmap(c)[2],cmap(c)[1],cmap(c)[0]])


        #     # print("AARON CMAP LOADED")
        #     # import pickle
        #     # with open('custom_cmap (4).p', 'rb') as file:
        #     #     cmap = pickle.load(file)
        #     cmap = create_smooth_diverging_cmap(cmap1="inferno",cmap2="viridis", cutoff=0.85, center_width=0.3, n=ncolors) # Increased width slightly for visibility
        #     cmap = truncate_colormap(cmap, 0.5, 1, n=ncolors)

        #     # cmap = create_smooth_diverging_cmap(cmap1="viridis", cmap2="inferno",cutoff=0.85, center_width=0.3, n=ncolors) # Increased width slightly for visibility
        #     # cmap = truncate_colormap(cmap, 0.1, 0.9, n=ncolors)

        #     # # --- Verification Plot ---
        #     # data = np.random.randn(20, 20)
        #     # plt.figure(figsize=(8, 6))
        #     # plt.imshow(data, cmap=truncated_cmap,vmin=-2,vmax=2)
        #     # plt.colorbar(label="Z Score")
        #     # plt.title(f"Smooth Diverging Map (Cutoff=0.8)")
        #     # plt.show()

        # elif 'aaron_diverging' == color:
        #     # print("AARON CMAP LOADED")
        #     # import pickle
        #     # with open('custom_cmap_diverging (2).p', 'rb') as file:
        #     #     cmap = pickle.load(file)
        #     cmap = create_smooth_diverging_cmap(cmap1="inferno",cmap2="viridis", cutoff=0.85, center_width=0.3, n=ncolors) # Increased width slightly for visibility
        #     cmap = truncate_colormap(cmap, 0.1, 0.9, n=ncolors)
    else:
        raise Exception("Unknown Color provided")

    return cmap, tempRGB, tempBGR
##################################################################################################
def export_colorbar(cmap_name,cmap,ncolors,cmin,cmax,figname,figdir):

    plt.rcParams['pdf.fonttype'] = 42
    # plt.rcParams['ps.fonttype'] = 'truetype'


    gradient = np.linspace(0, 1, ncolors)
    gradient = np.vstack((gradient, gradient))
    fig, axs = plt.subplots(nrows=1, figsize=(6, 1))
    fig.subplots_adjust(top=1 - 0.35 / 1, bottom=0.15 / 1,
                        left=0.2, right=0.99)
    axs.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap),interpolation='none')
    # axs.set_title(cmap_name, fontsize=20,fontname='arial')
    # axs.text(-0.01, 0.5, str(cmin), va='center', ha='right', fontsize=20,
    #         transform=axs.transAxes,fontname='arial')
    # axs.text(1.01, 0.5, str(cmax), va='center', ha='left', fontsize=20,
    #         transform=axs.transAxes,fontname='arial')
    axs.set_title(cmap_name, fontsize=20)
    axs.text(-0.01, 0.5, str(cmin), va='center', ha='right', fontsize=20,
            transform=axs.transAxes)
    axs.text(1.01, 0.5, str(cmax), va='center', ha='left', fontsize=20,
            transform=axs.transAxes)
    axs.set_axis_off()
    plt.savefig(os.path.join(figdir,figname+".pdf"),format="pdf",bbox_inches='tight',pad_inches=0.05,dpi=600)
    # plt.savefig(os.path.join(figdir,figname+".eps"),format="eps",bbox_inches='tight',pad_inches=0.05,dpi=600)
##################################################################################################
def export_colorbar_new(cmap,style,figName,figSaveDir, figsize = (10,2), fontsize=10, lw = 0.5, label = 'ΔF',minVal=np.nan,midVal=np.nan,maxVal=np.nan,minLabel='',midLabel='',maxLabel='',save_PDF=True):
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams['xtick.major.width'] = lw
    plt.rcParams['ytick.major.width'] = lw
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ["Arial"]
    plt.rcParams['pdf.fonttype'] = 42
    labelColor = (0,0,0)
    overlayColor = (1,1,1)
    if 'custom' in cmap and 'W' in cmap:
        overlayColor = (0,0,0)
    if style == 'centered':
        if np.isnan(minVal):
            minVal = -1
        if np.isnan(midVal):
            midVal = 0
        if np.isnan(maxVal):
            maxVal = 1
    elif style == 'end':
        if np.isnan(minVal):
            minVal = 0
        if np.isnan(midVal):
            midVal = 0.5
        if np.isnan(maxVal):
            maxVal = 1
    minVal = str(minVal)
    midVal = str(midVal)
    maxVal = str(maxVal)

    if style == 'circular' and cmap == 'hsv':
        size = 500
        # Create coordinate grid
        x = np.linspace(-1, 1, size)
        y = np.linspace(-1, 1, size)
        X, Y = np.meshgrid(x, y)
    
        # Convert to polar coordinates
        R = np.sqrt(X**2 + Y**2)
        Theta = np.arctan2(Y, X)
    
        # Normalize angle to [0, 1] for hue and reverse direction
        H = (1 - (Theta % (2 * np.pi)) / (2 * np.pi)) % 1.0
    
        # Normalize radius to [0, 1] for brightness
        V = np.clip(R, 0, 1)
    
        # Saturation is set to 1
        S = np.ones_like(H)
    
        # Convert HSV to RGB using matplotlib's hsv colormap
        HSV = np.stack((H, S, V), axis=-1)
        RGB = plt.cm.hsv(HSV[..., 0])[:, :, :3] * HSV[..., 2][..., np.newaxis]
    
        # Apply circular mask and set outside region to white
        mask = R <= 1
        RGB_masked = np.ones_like(RGB)  # Initialize with white
        for i in range(3):
            RGB_masked[..., i][mask] = RGB[..., i][mask]
    
        fig,ax = clean_subplots(1,1,figsize=figsize)
        ax.imshow(RGB_masked, extent=(-1, 1, -1, 1))
        ax.axis('off')
    
        # Add angle labels
        label_radius = 1.1
        angles_deg = [0, 90, 180, 270]
        for angle in angles_deg:
            rad = np.deg2rad(angle)
            x_pos = label_radius * np.cos(rad)
            y_pos = label_radius * np.sin(rad)
            ax.text(x_pos, y_pos, f'{angle}°', ha='center', va='center', fontsize=20, color='black')

        ax.plot([0,np.cos(45)], [0,np.sin(45)], color=(0,0,0), linewidth=lw)  # Vertical line
    else:
        xsize = 151
        ysize = 1001
        
        # if 'tuple' in str(type(cmap)):
        colormap,rgb,bgr=generate_cmap(cmap,ysize)

        # print(rgb)
        maxColor = rgb[-1][0:3]
        if len(rgb) > 5:
            midColor = rgb[int(ysize/2)][0:3]
        elif len(rgb) == 5:
            midColor = rgb[2][0:3]
        elif len(rgb) == 4:
            midColor = [(rgb[1][0]+rgb[2][0])/2,(rgb[1][1]+rgb[2][1])/2,(rgb[1][2]+rgb[2][2])/2]
        elif len(rgb) == 3:
            midColor = rgb[1][0:3]
        else:
            midColor = [(rgb[0][0]+rgb[1][0])/2,(rgb[0][1]+rgb[1][1])/2,(rgb[0][2]+rgb[1][2])/2]
        minColor = rgb[0][0:3]
        # print(f'maxColor = {maxColor}')
        # print(f'midColor = {midColor}')
        # print(f'minColor = {minColor}')
        if np.sum(maxColor) < 1.5:
            maxOverlayColor = (1,1,1)
        else:
            maxOverlayColor = (0,0,0)
        if np.sum(midColor) < 1.5:
            midOverlayColor = (1,1,1)
        else:
            midOverlayColor = (0,0,0)
        if np.sum(minColor) < 1.5:
            minOverlayColor = (1,1,1)
        else:
            minOverlayColor = (0,0,0)
        # print(f'maxOverlayColor = {maxOverlayColor}')
        # print(f'midOverlayColor = {midOverlayColor}')
        # print(f'minOverlayColor = {minOverlayColor}')
        x = np.linspace(-1, 1, xsize)
        y = np.linspace(-1, 1, ysize)
        X, Y = np.meshgrid(x, y)
        fig,ax = clean_subplots(1,6,figsize=figsize)

        col = 0
        ax[col].imshow(np.transpose(Y),interpolation='none',cmap=colormap)
        ax[col].set_facecolor('none')
        ax[col].set_yticklabels([])
        ax[col].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            left=False,      # ticks along the bottom edge are off
            right=False) # labels along the bottom edge are off

        ax[col].set_xlabel(label,fontsize=fontsize)
        if style == 'centered':
            xticks = np.array([0,ysize/2,ysize])
            ax[col].set_xticks(xticks)
            ax[col].set_xticklabels([minVal,midVal,maxVal],fontsize=fontsize,color=labelColor)
        else:
            xticks = np.array([0,ysize])
            ax[col].set_xticks(xticks)
            ax[col].set_xticklabels([minVal,maxVal],fontsize=fontsize,color=labelColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)

        col = 1
        ax[col].imshow(np.transpose(Y),interpolation='none',cmap=colormap)
        ax[col].set_facecolor('none')
        ax[col].set_yticklabels([])
        ax[col].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            left=False,      # ticks along the bottom edge are off
            right=False) # labels along the bottom edge are off
        ax[col].set_xlabel(label,fontsize=fontsize)
        if style == 'centered':
            xticks = np.array([0,ysize/2,ysize])
            ax[col].set_xticks([])
            # ax[col].set_xticklabels([minVal,midVal,maxVal],fontsize=fontsize)
            ax[col].text(0,xsize,minLabel,ha='left',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(0,xsize/2,minVal,ha='right',va='center',fontsize=fontsize,color=labelColor)
            ax[col].text(ysize/2,xsize/2,midVal,ha='center',va='center',fontsize=fontsize,color=midOverlayColor)
            ax[col].text(ysize,xsize,maxLabel,ha='right',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(ysize,xsize/2,maxVal,ha='left',va='center',fontsize=fontsize,color=labelColor)
        else:
            xticks = np.array([0,ysize])
            ax[col].set_xticks([])
            # ax[col].set_xticklabels([minVal,maxVal],fontsize=fontsize)
            ax[col].text(0,xsize/2,minVal,ha='right',va='center',fontsize=fontsize,color=labelColor)
            ax[col].text(ysize,xsize/2,maxVal,ha='left',va='center',fontsize=fontsize,color=labelColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)

        col = 2
        ax[col].imshow(np.transpose(Y),interpolation='none',cmap=colormap)
        ax[col].set_facecolor('none')
        ax[col].set_yticklabels([])
        ax[col].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            left=False,      # ticks along the bottom edge are off
            right=False) # labels along the bottom edge are off
        ax[col].set_xlabel(label,fontsize=fontsize)
        if style == 'centered':
            xticks = np.array([0,ysize/2,ysize])
            ax[col].set_xticks([])
            # ax[col].set_xticklabels([minVal,midVal,maxVal],fontsize=fontsize)
            ax[col].text(0,xsize,minLabel,ha='left',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(0,xsize/2,minVal,ha='left',va='center',fontsize=fontsize,color=minOverlayColor)
            ax[col].text(ysize/2,xsize/2,midVal,ha='center',va='center',fontsize=fontsize,color=midOverlayColor)
            ax[col].text(ysize,xsize,maxLabel,ha='right',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(ysize,xsize/2,maxVal,ha='right',va='center',fontsize=fontsize,color=maxOverlayColor)
        else:
            xticks = np.array([0,ysize])
            ax[col].set_xticks([])
            # ax[col].set_xticklabels([minVal,maxVal],fontsize=fontsize)
            ax[col].text(0,xsize/2,minVal,ha='left',va='center',fontsize=fontsize,color=minOverlayColor)
            ax[col].text(ysize,xsize/2,maxVal,ha='right',va='center',fontsize=fontsize,color=maxOverlayColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)


        col = 3
        ax[col].imshow(Y,interpolation='none',cmap=colormap)
        ax[col].invert_yaxis()
        ax[col].yaxis.tick_right()
        ax[col].set_facecolor('none')
        ax[col].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
        ax[col].set_ylabel(label,fontsize=fontsize)
        ax[col].yaxis.set_label_position("right")
        if style == 'centered':
            ax[col].text(0, ysize/2,label,ha='right',va='center',rotation=90,fontsize=fontsize,color=midOverlayColor)
            yticks = np.array([0,ysize/2,ysize])
            ax[col].set_yticks(yticks)
            ax[col].set_yticklabels([minVal,midVal,maxVal],fontsize=fontsize,color=labelColor)
        else:
            ax[col].text(xsize, ysize/2,label,ha='left',va='center',rotation=-90,fontsize=fontsize,color=midOverlayColor)
            yticks = np.array([0,ysize])
            ax[col].set_yticks(yticks)
            ax[col].set_yticklabels([minVal,maxVal],fontsize=fontsize,color=labelColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)



        col = 4
        ax[col].imshow(Y,interpolation='none',cmap=colormap)
        ax[col].invert_yaxis()
        ax[col].yaxis.tick_right()
        ax[col].set_facecolor('none')
        ax[col].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
        # ax[col].set_ylabel(label,fontsize=fontsize)
        if style == 'centered':
            yticks = np.array([0,ysize/2,ysize])
            ax[col].set_yticks([])
            # ax[col].set_yticklabels([minVal,midVal,maxVal],fontsize=fontsize)
            ax[col].text(xsize, 0,minLabel,ha='left',va='bottom',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize/2, 0,minVal,ha='center',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(xsize, ysize/2,label,ha='left',va='center',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize/2, ysize/2,midVal,ha='center',va='center',fontsize=fontsize,color=labelColor)
            ax[col].text(xsize, ysize,maxLabel,ha='left',va='top',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize/2, ysize,maxVal,ha='center',va='bottom',fontsize=fontsize,color=labelColor)
        else:
            yticks = np.array([0,ysize])
            ax[col].set_yticks([])
            # ax[col].set_yticklabels([minVal,maxVal],fontsize=fontsize)
            ax[col].text(xsize, 0,minLabel,ha='left',va='bottom',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize/2, 0,minVal,ha='center',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text(xsize, ysize/2,label,ha='left',va='center',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize, ysize,maxLabel,ha='left',va='top',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text(xsize/2, ysize,maxVal,ha='center',va='bottom',fontsize=fontsize,color=labelColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)

        col = 5
        ax[col].imshow(np.concatenate((Y,Y),axis = 1),interpolation='none',cmap=colormap)
        ax[col].invert_yaxis()
        ax[col].yaxis.tick_right()
        ax[col].set_facecolor('none')
        ax[col].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            width=lw,
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
        # ax[col].set_ylabel(label,fontsize=fontsize)
        if style == 'centered':
            yticks = np.array([0,ysize/2,ysize])
            ax[col].set_yticks([])
            # ax[col].set_yticklabels([minVal,midVal,maxVal],fontsize=fontsize)
            ax[col].text((2*xsize)/2, 0,minLabel,ha='center',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, 0,minVal,ha='center',va='bottom',fontsize=fontsize,color=minOverlayColor)
            ax[col].text((2*xsize), ysize/2,label,ha='left',va='center',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, ysize/2,midVal,ha='center',va='center',fontsize=fontsize,color=midOverlayColor)
            ax[col].text((2*xsize)/2, ysize,maxLabel,ha='center',va='bottom',fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, ysize,maxVal,ha='center',va='top',fontsize=fontsize,color=maxOverlayColor)
        else:
            yticks = np.array([0,ysize])
            ax[col].set_yticks([])
            # ax[col].set_yticklabels([minVal,maxVal],fontsize=fontsize)
            ax[col].text((2*xsize)/2, 0,minLabel,ha='center',va='top',fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, 0,minVal,ha='center',va='bottom',fontsize=fontsize,color=minOverlayColor)
            ax[col].text((2*xsize), ysize/2,label,ha='left',va='center',rotation=-90,fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, ysize,maxLabel,ha='center',va='bottom',fontsize=fontsize,color=labelColor)
            ax[col].text((2*xsize)/2, ysize,maxVal,ha='center',va='top',fontsize=fontsize,color=maxOverlayColor)
        ax[col].spines['left'].set_linewidth(lw)
        ax[col].spines['bottom'].set_linewidth(lw)
        ax[col].spines['right'].set_linewidth(lw)
        ax[col].spines['top'].set_linewidth(lw)

    if save_PDF:
        with PdfPages(os.path.join(figSaveDir,figName)) as pdf:
            pdf.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a pdf page
##################################################################################################
def ConvertColorLetter(tempColor):
    if len(tempColor)==1:
        tempRGB=(1, 1, 1)
        tempBGR=(255, 255, 255)
        if tempColor=='k':
            tempRGB=(0, 0, 0)
            tempBGR=(0, 0, 0)
        elif tempColor=='w':
            tempRGB=(1, 1, 1)
            tempBGR=(255, 255, 255)
        elif tempColor=='r':
            tempRGB=(1, 0, 0)                   
            tempBGR=(0, 0, 255)
        elif tempColor=='g':
            tempRGB=(0, 1, 0)
            tempBGR=(0, 255, 0)
        elif tempColor=='b':
            tempRGB=(0, 0, 1)                           
            tempBGR=(255, 0, 0)
        elif tempColor=='c':
            tempRGB=(0, 1, 1)
            tempBGR=(255, 255, 0)
        elif tempColor=='m':
            tempRGB=(1, 0, 1)                        
            tempBGR=(255, 0, 255)
        elif tempColor=='y':
            tempRGB=(1, 1, 0)
            tempBGR=(0, 255, 255)
        else:
            tempRGB=tempColor
            tempBGR=tempColor
    return tempRGB,tempBGR
##################################################################################################
def convert2RGB(input_frame,cmin,cmax,cmap,ColorScalar,nan_color=[0.5,0.5,0.5]):
    frame=copy.deepcopy(input_frame)
    frame_height,frame_width=np.shape(frame)
    frame_RGB=np.zeros((frame_height,frame_width,3),dtype='float32')
    try:
        frame[frame<float(cmin)]=float(cmin)
        frame[frame>float(cmax)]=float(cmax)
        frame=frame-float(cmin)
        if (float(cmax)-float(cmin))>0:
            frame=frame/(float(cmax)-float(cmin))
        frame_RGB=cmap(frame)
        # frame=frame.astype('float32')
        # frame[frame<float(cmin)]=float(cmin)
        # frame[frame>float(cmax)]=float(cmax)
        # frame=frame-np.min(frame)
        # frame=frame/np.max(frame)
        # # frame=frame*ColorScalar
        # # frame=np.round(frame)
        # # frame=frame/np.max(frame)
        # frame_RGB=cmap(frame)
        if np.any(np.isnan(frame)):
            frame_mask=np.isnan(frame)
            frame_RGB_R=frame_RGB[:,:,0]
            frame_RGB_R[frame_mask]=float(nan_color[0])
            frame_RGB[:,:,0]=frame_RGB_R
            frame_RGB_G=frame_RGB[:,:,1]
            frame_RGB_G[frame_mask]=float(nan_color[1])
            frame_RGB[:,:,1]=frame_RGB_G
            frame_RGB_B=frame_RGB[:,:,2]
            frame_RGB_B[frame_mask]=float(nan_color[2])
            frame_RGB[:,:,2]=frame_RGB_B
            frame_RGB[:,:,3]=1
    except:
        pass
    return frame_RGB
##################################################################################################
def convert2RGB_px_colormap(input_frame,cmin,cmax,color_scalar,px_cmap,nan_color=[0.5,0.5,0.5]):
    frame=copy.deepcopy(input_frame)
    frame_height,frame_width=np.shape(frame)
    frame_RGB=np.zeros((frame_height,frame_width,4),dtype='float32')
    # try:
    frame[frame<float(cmin)]=float(cmin)
    frame[frame>float(cmax)]=float(cmax)
    frame=frame-float(cmin)
    frame=frame/(float(cmax)-float(cmin))
    frame=np.round(frame*color_scalar).astype('int')
    if np.any(frame>0):
        # for y in range(frame_height):
        #     for x in range(frame_width):
        #         if np.any(frame[y,x]>0):
        #             frame_RGB[y,x,:]=px_cmap[y,x,frame[y,x],:]
        pxs=np.array(frame.nonzero())
        npx=pxs.shape[1]
        for px in range(npx):
            frame_RGB[pxs[0,px],pxs[1,px],:]=px_cmap[pxs[0,px],pxs[1,px],frame[pxs[0,px],pxs[1,px]]-1,:]
    if np.any(np.isnan(frame)):
        frame_mask=np.isnan(frame)
        frame_RGB_R=frame_RGB[:,:,0]
        frame_RGB_R[frame_mask]=float(nan_color[0])
        frame_RGB[:,:,0]=frame_RGB_R
        frame_RGB_G=frame_RGB[:,:,1]
        frame_RGB_G[frame_mask]=float(nan_color[1])
        frame_RGB[:,:,1]=frame_RGB_G
        frame_RGB_B=frame_RGB[:,:,2]
        frame_RGB_B[frame_mask]=float(nan_color[2])
        frame_RGB[:,:,2]=frame_RGB_B
        frame_RGB[:,:,3]=1
    # except:
    #     pass
    return frame_RGB
##################################################################################################
def add_plot_scaleBar(ax, plot_scaleBar, color = (0,0,0),includeLabel = False, vert_unit = 'AU',horz_unit = 's'):
    if not 'includeVert' in plot_scaleBar:
        plot_scaleBar['includeVert'] = True
    if not 'includeHorz' in plot_scaleBar:
        plot_scaleBar['includeHorz'] = True
    if not 'horzLabel_vertAdjust' in plot_scaleBar:
        plot_scaleBar['horzLabel_vertAdjust'] = 0.03
    if not 'vertLabel_horzAdjust' in plot_scaleBar:
        plot_scaleBar['vertLabel_horzAdjust'] = 0.02
    if not 'includeLabel' in plot_scaleBar:
        plot_scaleBar['includeLabel'] = False
    if not 'fontsize' in plot_scaleBar:
        plot_scaleBar['fontsize'] = 10

    xlim = ax.get_xlim()    
    ylim = ax.get_ylim()
    sizeX = xlim[1] - xlim[0]
    sizeY = ylim[1] - ylim[0]
    vertRatio=sizeY*plot_scaleBar['vertAdjust']
    horzRatio=sizeX*plot_scaleBar['horzAdjust']
    if 'manual_height' in plot_scaleBar:
        height = copy.deepcopy(plot_scaleBar['manual_height'])
        if sizeY < 0:
            height = -height
    else:
        height = sizeY*plot_scaleBar['height_per']
        height = smart_floor(height)
    if 'length_fr' in plot_scaleBar and not 'length' in plot_scaleBar:
        horz_length = copy.deepcopy(plot_scaleBar['length_fr'])
    elif 'length_s' in plot_scaleBar and not 'length_fr' in plot_scaleBar and not 'length' in plot_scaleBar:
        horz_length = copy.deepcopy(plot_scaleBar['length_s'])
    elif 'length' in plot_scaleBar:
        horz_length = copy.deepcopy(plot_scaleBar['length'])
    else:
        raise Exception("PROBLEM WITH HORZONTAL SCALE BAR LENGTH SPECIFICATION")

    if plot_scaleBar['corner'] == 'BL':
        plot_scaleBar['horz_xcoords']=[xlim[0]+horzRatio,xlim[0]+horzRatio+horz_length]
        plot_scaleBar['horz_ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio]
        plot_scaleBar['horz_labelcoord']=[np.mean(plot_scaleBar['horz_xcoords']),plot_scaleBar['horz_ycoords'][0]-vertRatio*plot_scaleBar['horzLabel_vertAdjust']]
        plot_scaleBar['vert_xcoords']=[xlim[0]+horzRatio,xlim[0]+horzRatio]
        plot_scaleBar['vert_ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio+height]
        plot_scaleBar['vert_labelcoord']=[plot_scaleBar['vert_xcoords'][0]+horzRatio*plot_scaleBar['vertLabel_horzAdjust'], np.mean(np.mean(plot_scaleBar['vert_ycoords']))]
    elif plot_scaleBar['corner'] == 'BR':
        plot_scaleBar['horz_xcoords']=[xlim[1]-horzRatio,xlim[1]-horzRatio-horz_length]
        plot_scaleBar['horz_ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio]
        plot_scaleBar['horz_labelcoord']=[np.mean(plot_scaleBar['horz_xcoords']),plot_scaleBar['horz_ycoords'][0]-vertRatio*plot_scaleBar['horzLabel_vertAdjust']]
        plot_scaleBar['vert_xcoords']=[xlim[1]-horzRatio,xlim[1]-horzRatio]
        plot_scaleBar['vert_ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio+height]
        plot_scaleBar['vert_labelcoord']=[plot_scaleBar['vert_xcoords'][0]-horzRatio*plot_scaleBar['vertLabel_horzAdjust'], np.mean(np.mean(plot_scaleBar['vert_ycoords']))]
    elif plot_scaleBar['corner'] == 'TL':
        plot_scaleBar['horz_xcoords']=[xlim[0]+horzRatio,xlim[0]+horzRatio+horz_length]
        plot_scaleBar['horz_ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio]
        plot_scaleBar['horz_labelcoord']=[np.mean(plot_scaleBar['horz_xcoords']),plot_scaleBar['horz_ycoords'][0]+vertRatio*plot_scaleBar['horzLabel_vertAdjust']]
        plot_scaleBar['vert_xcoords']=[xlim[0]+horzRatio,xlim[0]+horzRatio]
        plot_scaleBar['vert_ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio-height]
        plot_scaleBar['vert_labelcoord']=[plot_scaleBar['vert_xcoords'][0]+horzRatio*plot_scaleBar['vertLabel_horzAdjust'], np.mean(np.mean(plot_scaleBar['vert_ycoords']))]
    elif plot_scaleBar['corner'] == 'TR':
        plot_scaleBar['horz_xcoords']=[xlim[1]-horzRatio,xlim[1]-horzRatio-horz_length]
        plot_scaleBar['horz_ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio]
        plot_scaleBar['horz_labelcoord']=[np.mean(plot_scaleBar['horz_xcoords']),plot_scaleBar['horz_ycoords'][0]+vertRatio*plot_scaleBar['horzLabel_vertAdjust']]
        plot_scaleBar['vert_xcoords']=[xlim[1]-horzRatio,xlim[1]-horzRatio]
        plot_scaleBar['vert_ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio-height]
        plot_scaleBar['vert_labelcoord']=[plot_scaleBar['vert_xcoords'][0]-horzRatio*plot_scaleBar['vertLabel_horzAdjust'], np.mean(np.mean(plot_scaleBar['vert_ycoords']))]

    if plot_scaleBar['includeHorz']:
        ax.plot(plot_scaleBar['horz_xcoords'],plot_scaleBar['horz_ycoords'],linestyle = '-',color=color,linewidth=plot_scaleBar['lw'],solid_capstyle='butt')
    if plot_scaleBar['includeVert']:
        ax.plot(plot_scaleBar['vert_xcoords'],plot_scaleBar['vert_ycoords'],linestyle = '-',color=color,linewidth=plot_scaleBar['lw'],solid_capstyle='butt')



    if 'length_fr' in plot_scaleBar and 'fr' in horz_unit:
        horz_label = str(copy.deepcopy(plot_scaleBar['length_fr']))+' '+horz_unit
    elif 'length_s' in plot_scaleBar and 's' in horz_unit:
        horz_label = str(copy.deepcopy(plot_scaleBar['length_s']))+' '+horz_unit
    elif len(horz_unit) == 0:
        horz_label = str(copy.deepcopy(plot_scaleBar['length']))
    else:
        horz_label = horz_unit

    height = np.abs(height)
    if height>10 and np.mod(height,10) == 0:
        height = int(height)
    vert_label = str(height)+' '+vert_unit

    if plot_scaleBar['includeLabel'] or includeLabel:
        if plot_scaleBar['includeHorz']:
            ax.text(plot_scaleBar['horz_labelcoord'][0],plot_scaleBar['horz_labelcoord'][1],\
                    horz_label,color=color,fontsize=plot_scaleBar['fontsize'],ha='center',va='top')
        if plot_scaleBar['includeVert']:
            if plot_scaleBar['corner'] == 'BL':
                ax.text(plot_scaleBar['vert_labelcoord'][0],plot_scaleBar['vert_labelcoord'][1],\
                        vert_label,color=color,fontsize=plot_scaleBar['fontsize'],ha='left',va='center')
            elif plot_scaleBar['corner'] == 'BR':
                ax.text(plot_scaleBar['vert_labelcoord'][0],plot_scaleBar['vert_labelcoord'][1],\
                        vert_label,color=color,fontsize=plot_scaleBar['fontsize'],ha='right',va='center')
            elif plot_scaleBar['corner'] == 'TL':
                ax.text(plot_scaleBar['vert_labelcoord'][0],plot_scaleBar['vert_labelcoord'][1],\
                        vert_label,color=color,fontsize=plot_scaleBar['fontsize'],ha='left',va='center')
            elif plot_scaleBar['corner'] == 'TR':
                ax.text(plot_scaleBar['vert_labelcoord'][0],plot_scaleBar['vert_labelcoord'][1],\
                        vert_label,color=color,fontsize=plot_scaleBar['fontsize'],ha='right',va='center')
    return ax,plot_scaleBar
##################################################################################################
##################################################################################################
def pickler(savePath, saveName, data, extension='.pkl', verbose=True):
    """Save an object to a .pkl/.p file, optionally printing file details.

    This is a lightweight counterpart to safe_pickler: it has no scratch-directory
    "safe save" behavior. You give it a filename and object and it writes the
    file directly, while (optionally) reporting the file size and save time.

    Parameters
    ----------
    savePath : str
        The path where the file will be saved.
    saveName : str
        The name of the file to be saved. May include a '.pkl' or '.p'
        extension; if no extension is present, `extension` is appended.
    data : object
        The object to be saved.
    extension : str
        The file extension to append to `saveName` when it has none. Defaults to '.pkl'.
    verbose : bool
        If True, print the file size and save time.

    Returns
    -------
    None

    Example
    -------
    savePath = '/some/path'
    saveName = 'my_data'
    data = my_data                 # the object to pickle
    extension = '.pkl'             # appended to saveName only if it has no extension
    verbose = True                 # True | False
    pickler(
        savePath, saveName, data, extension=extension, verbose=verbose)
    """

    import pickle
    import timeit
    save_timer = timeit.default_timer()
    savePath_str = copy.deepcopy(savePath)
    savePath_str = savePath_str.replace("\\", "/")

    saveName_lower = saveName.lower()
    if saveName_lower.endswith('.pkl'):
        final_save_name = saveName
    elif saveName_lower.endswith('.p'):
        final_save_name = saveName
    else:
        final_save_name = saveName + extension

    full_path = os.path.join(savePath, final_save_name)
    if verbose:
        print("  Path: " + savePath_str)
        print("  File: " + final_save_name)

    if not os.path.exists(savePath):
        print("<<WARNING>> PATH DOESNT EXIST")
        print("<<WARNING>> PATH DOESNT EXIST")
        print("<<WARNING>> PATH DOESNT EXIST")

    if verbose:
        print("  Saving " + final_save_name + " to " + savePath_str, end='...', flush=True)
    f = open(full_path, 'wb')
    pickle.dump(data, f)
    f.close()

    if not os.path.exists(full_path):
        print("<<WARNING>> FILE DOESNT EXIST")
        print("<<WARNING>> FILE DOESNT EXIST")
        print("<<WARNING>> FILE DOESNT EXIST")
    else:
        if verbose:
            print('Finished!', flush=True)
            fileSize_GB = os.path.getsize(full_path) / 1e9
            fileSize_GiB = fileSize_GB * 0.931323
            print("  FileSize:      " + str(np.round(fileSize_GB, decimals=2)) + " GB (" + str(np.round(fileSize_GiB, decimals=2)) + " GiB)")
            save_time = timeit.default_timer() - save_timer
            print("  Saving Took a total of : " + str((np.round((save_time / 60) * 1000) / 1000)) + " min", flush=True)
def unpickler(savePath, saveName, extension='.pkl', verbose=True):
    """Load an object from a .pkl/.p file, optionally printing file details.

    This is a lightweight counterpart to safe_unpickler: it has no scratch-directory
    "safe load" behavior. You give it a filename and it returns the loaded data,
    while (optionally) reporting the file size, last-modified time, and load time.

    Parameters
    ----------
    savePath : str
        The path where the file is saved.
    saveName : str
        The name of the file to be loaded. May include a '.pkl' or '.p' extension;
        if no extension is present, `extension` is appended.
    extension : str
        The file extension to append to `saveName` when it has none. Defaults to '.pkl'.
    verbose : bool
        If True, print the file size, last-modified time, and load time.

    Returns
    -------
    data : object
        The loaded object (an empty dict {} if the file does not exist).

    Example
    -------
    savePath = '/some/path'
    saveName = 'my_data'
    extension = '.pkl'
    verbose = True
    data = unpickler(
        savePath, saveName, extension=extension, verbose=verbose
    )
    """

    import pickle
    import timeit
    load_timer = timeit.default_timer()
    savePath_str = copy.deepcopy(savePath)
    savePath_str = savePath_str.replace("\\", "/")

    saveName_lower = saveName.lower()
    if saveName_lower.endswith('.pkl'):
        final_save_name = saveName
    elif saveName_lower.endswith('.p'):
        final_save_name = saveName
    else:
        final_save_name = saveName + extension

    full_path = os.path.join(savePath, final_save_name)
    if verbose:
        print("  Path: " + savePath_str)
        print("  File: " + final_save_name)

    if not os.path.exists(full_path):
        print("<<WARNING>> FILE DOESNT EXIST")
        print("<<WARNING>> FILE DOESNT EXIST")
        print("<<WARNING>> FILE DOESNT EXIST")
        data = {}
    else:
        if verbose:
            fileSize_GB = os.path.getsize(full_path) / 1e9
            fileSize_GiB = fileSize_GB * 0.931323
            print("  FileSize:      " + str(np.round(fileSize_GB, decimals=2)) + " GB (" + str(np.round(fileSize_GiB, decimals=2)) + " GiB)")
            t = os.path.getmtime(full_path)
            print("  Last Modified: " + str(datetime.fromtimestamp(t)))
            print("  Loading " + final_save_name + " from " + savePath_str, end='...', flush=True)
        f = open(full_path, 'rb')
        data = pickle.load(f)
        f.close()
        if verbose:
            print('Finished!', flush=True)
            load_time = timeit.default_timer() - load_timer
            print("  Loading Took a total of : " + str((np.round((load_time / 60) * 1000) / 1000)) + " min", flush=True)
    return data
##################################################################################################

def zoom_lookup(zoom_s,figParams,verbose=True):

    if not 'plot_scaleBar' in figParams:
        figParams['plot_scaleBar'] = {}


    figParams['xticks'] = {}
    figParams['xticks']['trialStart'] = [-6,-4,-2,0,2,4,6,8,10]
    figParams['xticks']['goCue'] = [-8,-6,-4,-2,0,2,4,6,8]
    figParams['xticks']['cTimes'] = [-8,-6,-4,-2,0,2,4,6,8]
    figParams['xticks']['postContLick1'] = [-8,-6,-4,-2,0,2,4,6,8]
    figParams['xdata_correction'] = {}
    figParams['xdata_correction']['trialStart'] = False
    figParams['xdata_correction']['goCue'] = False
    figParams['xdata_correction']['cTimes'] = False
    figParams['xdata_correction']['postContLick1'] = False
    figParams['plot_scaleBar']['length_s'] = 2
    figParams['xdata_correction'] = {}
    figParams['xticks_zoom'] = {}
    figParams['xlim_zoom'] = {}
    figParams['zoom_label'] = ''
        
    if zoom_s == 6:
        zoom = True
        figParams['zoom_label'] = '6sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4,5,6]
        figParams['xlim_zoom']['trialStart'] = [-0.4,6.2]
        figParams['xticks_zoom']['goCue'] = [-3,-2,-1,0,1,2,3]
        figParams['xlim_zoom']['goCue'] = [-3.1,3.1]
        figParams['xticks_zoom']['cTimes'] = [-3,-2,-1,0,1,2,3]
        figParams['xlim_zoom']['cTimes'] = [-3.1,3.1]
        figParams['xticks_zoom']['postContLick1'] = [-3,-2,-1,0,1,2,3]
        figParams['xlim_zoom']['postContLick1'] = [-3.1,3.1]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == -6:
        zoom = True
        figParams['zoom_label'] = '6sZOOMOffset'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4,5,6]
        figParams['xlim_zoom']['trialStart'] = [-0.4,6.2]
        figParams['xticks_zoom']['goCue'] = [-1,0,1,2,3,4,5]
        figParams['xlim_zoom']['goCue'] = [-1.1,5.1]
        figParams['xticks_zoom']['cTimes'] = [-1,0,1,2,3,4,5]
        figParams['xlim_zoom']['cTimes'] = [-1.1,5.1]
        figParams['xticks_zoom']['postContLick1'] = [-1,0,1,2,3,4,5]
        figParams['xlim_zoom']['postContLick1'] = [-1.1,5.1]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == 5:
        zoom = True
        figParams['zoom_label'] = '4sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4,5]
        figParams['xlim_zoom']['trialStart'] = [-0.4,5.2]
        figParams['xticks_zoom']['goCue'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['goCue'] = [-2.6,2.6]
        figParams['xticks_zoom']['cTimes'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['cTimes'] = [-2.6,2.6]
        figParams['xticks_zoom']['postContLick1'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['postContLick1'] = [-2.6,2.6]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == -5:
        zoom = True
        figParams['zoom_label'] = '4sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4,5]
        figParams['xlim_zoom']['trialStart'] = [-0.4,5.2]
        figParams['xticks_zoom']['goCue'] = [-1,0,1,2,3,4]
        figParams['xlim_zoom']['goCue'] = [-1.1,4.1]
        figParams['xticks_zoom']['cTimes'] = [-1,0,1,2,3,4]
        figParams['xlim_zoom']['cTimes'] = [-1.1,4.1]
        figParams['xticks_zoom']['postContLick1'] = [-1,0,1,2,3,4]
        figParams['xlim_zoom']['postContLick1'] = [-1.1,4.1]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == 4:
        zoom = True
        figParams['zoom_label'] = '4sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4]
        figParams['xlim_zoom']['trialStart'] = [-0.4,4.2]
        figParams['xticks_zoom']['goCue'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['goCue'] = [-2.1,2.1]
        figParams['xticks_zoom']['cTimes'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['cTimes'] = [-2.1,2.1]
        figParams['xticks_zoom']['postContLick1'] = [-2,-1,0,1,2]
        figParams['xlim_zoom']['postContLick1'] = [-2.1,2.1]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == -4:
        zoom = True
        figParams['zoom_label'] = '4sZOOMOffset'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3,4]
        figParams['xlim_zoom']['trialStart'] = [-0.4,4.2]
        figParams['xticks_zoom']['goCue'] = [0,1,2,3]
        figParams['xlim_zoom']['goCue'] = [-0.5,3.5]
        figParams['xticks_zoom']['cTimes'] = [0,1,2,3]
        figParams['xlim_zoom']['cTimes'] = [-0.5,3.5]
        figParams['xticks_zoom']['postContLick1'] = [0,1,2,3]
        figParams['xlim_zoom']['postContLick1'] = [-0.5,3.5]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 1
    elif zoom_s == 3:
        zoom = True
        figParams['zoom_label'] = '3sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3]
        figParams['xlim_zoom']['trialStart'] = [-0.4,3.2]
        figParams['xticks_zoom']['goCue'] = [-1.5,-1,-0.5,0,0.5,1,1.5]
        figParams['xlim_zoom']['goCue'] = [-1.6,1.6]
        figParams['xticks_zoom']['cTimes'] = [-1.5,-1,-0.5,0,0.5,1,1.5]
        figParams['xlim_zoom']['cTimes'] = [-1.6,1.6]
        figParams['xticks_zoom']['postContLick1'] = [-1.5,-1,-0.5,0,0.5,1,1.5]
        figParams['xlim_zoom']['postContLick1'] = [-1.6,1.6]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 0.5
    elif zoom_s == -3:
        zoom = True
        figParams['zoom_label'] = '3sZOOM'
        figParams['xticks_zoom']['trialStart'] = [0,1,2,3]
        figParams['xlim_zoom']['trialStart'] = [-0.4,3.2]
        figParams['xticks_zoom']['goCue'] = [0,0.5,1,1.5,2,2.5,3]
        figParams['xlim_zoom']['goCue'] = [-0.2,3.2]
        figParams['xticks_zoom']['cTimes'] =  [0,0.5,1,1.5,2,2.5,3]
        figParams['xlim_zoom']['cTimes'] = [-0.2,3.2]
        figParams['xticks_zoom']['postContLick1'] =  [0,0.5,1,1.5,2,2.5,3]
        figParams['xlim_zoom']['postContLick1'] = [-0.2,3.2]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 0.5
    elif zoom_s == 2.5:
        zoom = True
        figParams['zoom_label'] = '2_5sZOOM'
        figParams['xticks_zoom']['trialStart'] = [2,3]
        figParams['xlim_zoom']['trialStart'] = [1.35,3.85]
        figParams['xticks_zoom']['goCue'] = [-1,-0.5,0,0.5,1]
        figParams['xlim_zoom']['goCue'] = [-1.25,1.25]
        figParams['xticks_zoom']['cTimes'] = [-1,-0.5,0,0.5,1]
        figParams['xlim_zoom']['cTimes'] = [-1.25,1.25]
        figParams['xticks_zoom']['postContLick1'] = [-1,-0.5,0,0.5,1]
        figParams['xlim_zoom']['postContLick1'] = [-1.25,1.25]
        figParams['xdata_correction']['trialStart'] = False
        figParams['xdata_correction']['goCue'] = False
        figParams['xdata_correction']['cTimes'] = False
        figParams['xdata_correction']['postContLick1'] = False
        figParams['plot_scaleBar']['length_s'] = 0.25
    else:
        zoom = False

    if zoom and verbose:
        print("Zoom mode activated: "+figParams['zoom_label'])
        print("          xlim_zoom: "+str(figParams['xlim_zoom']))
        print("        xticks_zoom: "+str(figParams['xticks_zoom']))

    return figParams,zoom
#########################################################################################################


def aligned_ROI_cluster_trace_raster_fig(ROI_clustering,figParams,ROI_type_key,ROI_score,align_data,zoom):
    if not 'clean' in figParams:
        figParams['clean'] = False
    if not 'clusterList' in figParams:
        figParams['clusterList'] = []
    if len(figParams['clusterList']) == 0:
        figParams['clusterList'] = list(range(ROI_clustering[ROI_type_key][ROI_score]['n_clusters']))
    if not 'ylim_adjust' in figParams:
        figParams['ylim_adjust'] = [-5,5]
    if figParams['group'] == 'trialStart':
        groupLabel = "Trial Start Aligned"
    elif figParams['group'] == 'cTimes':
        groupLabel = "First Lick Aligned"
    elif figParams['group'] == 'goCue':
        groupLabel = "Go Cue Aligned"
    elif figParams['group'] == 'postContLick1':
        groupLabel = "Second Lick Aligned"


    a = figParams['export_align_data'].index(align_data)
    align_data_short = figParams['export_align_data_label'][a]
    k = figParams['export_ROI_type_keys'].index(ROI_type_key)

    total = 0
    anmList = []
    for c0,c in enumerate(figParams['clusterList']):
        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
            total = total + ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']
            for r in range(len(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['ROI_info'])):
                if not ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['ROI_info'][r][3] in anmList:
                    anmList.append(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['ROI_info'][r][3])

    cmap, _, _ = generate_cmap(figParams['cmap'],1000,cmap_name=str(figParams['cmap']))
    cmap.set_bad(color=figParams['nan_color'])
    
    if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 0 or \
        ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
        binLabel = ""
    else:
        binLabel = "bin"+str(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key])

    trace_label = figParams['group']+" "+align_data_short+" "+binLabel+" "+figParams['plot_data_label_full']
    trace_label1 = figParams['group']+"_"+align_data_short+"_"+binLabel+"_"+figParams['plot_data_label_full']+"_"+figParams['plot_data_label']
    
    trace_label_byROI = figParams['group']+" "+align_data_short+" "+binLabel+" "+figParams['plot_data_byROI_label_full']
    trace_label_byROI1 = figParams['group']+"_+"+align_data_short+"_"+binLabel+"_"+figParams['plot_data_byROI_label_full']+"_"+figParams['plot_data_byROI_label']
    ##################################
    nRows=1
    nCols=len(figParams['export_ttypes'])
    figsize=(nCols*figParams['page2_hscalar'], nRows*figParams['page2_vscalar'])
    if nCols == 1:
        nCols+=1
    if figParams['clean']:
        print(f'figsize: {figsize} nRows: {nRows} nCols: {nCols}')
    
    fig,ax=clean_subplots(nRows,nCols,figsize=figsize,constrained_layout = False)
    maxVal=0
    minVal=0
    for c0,c in enumerate(figParams['clusterList']):
        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
            if 'allSess' in figParams['trialGrouping']:
                for t1,ttype in enumerate(list(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][figParams['trialGrouping']].keys())):
                    if ttype in figParams['scaling_ttypes']:
                        for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                            maxVal=np.nanmax([maxVal,np.nanpercentile(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][\
                                figParams['trialGrouping']][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:],95)])
                            minVal=np.nanmin([minVal,np.nanpercentile(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][\
                                figParams['trialGrouping']][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:],0)])
            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                for b,be in enumerate(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][figParams['trialGrouping']].keys()):
                    if be in figParams['scaling_be']:
                        for t1,ttype in enumerate(list(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][figParams['trialGrouping']][be].keys())):
                            if ttype in figParams['scaling_ttypes']:
                                for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                    maxVal=np.nanmax([maxVal,np.nanpercentile(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][\
                                        figParams['trialGrouping']][be][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:],95)])
                                    minVal=np.nanmin([minVal,np.nanpercentile(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['group']][align_data][\
                                        figParams['trialGrouping']][be][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:],0)])

    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
        load_alignment_xdata(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams'],ROI_type_key,figParams['group'],zoom,figParams)

    for t1,ttype in enumerate(figParams['export_ttypes']):
        emptyRow = np.ones((1,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['group']][align_data]),dtype='float32')*np.nan
        edges = [0]
        yticks = []
        yticklabels = []
        ytickcolors = []
        if figParams['addSpacers']:
            allClusterData = copy.deepcopy(emptyRow)
        else:
            allClusterData = np.zeros((0,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['group']][align_data]),dtype='float32')
        for c0,c in enumerate(figParams['clusterList']):
            if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                yticks.append(edges[-1]+ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']/2)
                yticklabels.append(ROI_clustering[ROI_type_key][ROI_score]['cluster_labels'][c]+"\n("+\
                    str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])+", "+str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['perROIs'])+"%)")
                for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                    roiTrace =  extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['group'],figParams['trialGrouping'],align_data,figParams['trialGrouping_behavior_epoch'],ttype,r,figParams['fillNaNs'])
                    allClusterData = np.concatenate((allClusterData,\
                                    np.expand_dims(roiTrace,axis=0)),axis=0)
                if figParams['addSpacers']:
                    allClusterData = np.concatenate((allClusterData,emptyRow),axis=0)
                    edges.append(allClusterData.shape[0]-1)
                else:
                    edges.append(allClusterData.shape[0])
                ytickcolors.append(ROI_clustering[ROI_type_key][ROI_score]['cluster_colors'][c])   # First label
        
        if 'manualContrastLimits' in figParams and figParams['manualContrastLimits'] is not None:
            vmin = figParams['manualContrastLimits'][0]
            vmax = figParams['manualContrastLimits'][1]
        else:
            if 'delta' in figParams['trialGrouping']:
            
                vmin=-1*np.nanmax([np.abs(maxVal),np.abs(minVal)])*figParams['imScalar_byROI'][a]
                vmax=np.nanmax([np.abs(maxVal),np.abs(minVal)])*figParams['imScalar_byROI'][a]
                # print("Delta scaling: "+str(vmin)+" to "+str(vmax))
            else:
                vmin=0
                vmax=maxVal*figParams['imScalar_byROI'][a]
                print("maxVal: "+str(maxVal))
                print("Positive scaling: "+str(vmin)+" to "+str(vmax))
            
        
        allClusterData_RGB = convert2RGB(allClusterData,\
            vmin,vmax,\
            cmap,figParams['colorScalar'],figParams['nan_color'])
        allClusterData_RGB = allClusterData_RGB[:,:,0:3]

        if figParams['ROI_aspects'][k] == 'auto':
            im=ax[t1].imshow(allClusterData_RGB[:,plotIdxs,:],extent=[xlim_s[0],xlim_s[1],allClusterData.shape[0],0],\
                interpolation='none',aspect='auto',clip_on=False)
        else:
            im=ax[t1].imshow(allClusterData_RGB[:,plotIdxs,:],extent=[xlim_s[0],xlim_s[1],allClusterData.shape[0],0],\
                interpolation='none',aspect=(figParams['ROI_aspects'][k]*allClusterData.shape[1])/allClusterData.shape[0],clip_on=False)
        ax[t1].spines['top'].set_visible(False)
        ax[t1].spines['right'].set_visible(False)
        ax[t1].spines['bottom'].set_visible(False)
        ax[t1].spines['left'].set_visible(False)
        ax[t1].set_facecolor('none')
        if figParams['clusterGrid']:
            ax[t1].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
            ax[t1].set_yticks([])
            for e in edges:
                if figParams['addSpacers']:
                    ax[t1].axhline(e,color=figParams['imshow_lineColor'],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                else:
                    ax[t1].axhline(e-0.5,color=figParams['imshow_lineColor'],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
        else:
            ax[t1].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
            ax[t1].set_yticks([])
            for c,yticklabel in enumerate(yticklabels):
                if figParams['addSpacers']:
                    ax[t1].plot([xlim_s[0],xlim_s[0]],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                else:
                    ax[t1].plot([xlim_s[0],xlim_s[0]],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
        if figParams['addClusterLabels'] and t1 == 0:
            ax[t1].set_yticks(yticks)
            ax[t1].tick_params(axis='y', which='both',length=0)
            ax[t1].set_yticklabels(yticklabels,fontdict = {'fontsize': figParams['page2_fontSize']-2,\
                'verticalalignment': figParams['imshow_ytick_verticalalignment'],\
                    'horizontalalignment': figParams['imshow_ytick_horizontalalignment'],\
                        'rotation':figParams['imshow_ytick_rotation']})
            ytick_labels = ax[t1].get_yticklabels()
            for c in range(len(ytickcolors)):
                ytick_labels[c].set_color(ytickcolors[c])   # First label

        ax[t1].set_xlim(xlim_s)
        if not figParams['clean']:
            ax[t1].set_xticks(xticks_s)
            ax[t1].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['page2_fontSize']-2,'verticalalignment': 'top','horizontalalignment': 'center'})
        else:
            ax[t1].set_xticks([])
        if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 0 or \
            ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
            ax[t1].set_xlabel(groupLabel+"(s)",fontsize=figParams['page2_fontSize'],color=(0,0,0))
        else:
            ax[t1].set_xlabel(groupLabel+"(s; B"+\
                str(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key])+")",\
                    fontsize=figParams['page2_fontSize'],color=(0,0,0))

        figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
        ax[t1],figParams['im_scaleBar'] = add_plot_scaleBar(ax[t1],figParams['im_scaleBar'],figParams['imshow_lineColor'],True,' '+figParams['im_scaleBar']['vertUnit'],str(figParams['plot_scaleBar']['length_s'])+' '+figParams['im_scaleBar']['horzUnit'])
        trial_structure_times = load_clustering_trial_structure_times(ROI_clustering,ROI_type_key,ROI_score,figParams['group'],align_data,fix_overlaps=True,verbose = False)
        ax[t1] = add_trial_structure_features('time',ax[t1],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['group'],align_data,ttype,[],\
            True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['page2_fontSize'],0,0-figParams['ylim_adjust'][0],\
                figParams['horzCueLineOn'],figParams['horzCueImScalar'])
        
        ax[t1].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])

        if figParams['group'] == ROI_clustering[ROI_type_key][ROI_score]['clustering_group'] and align_data == ROI_clustering[ROI_type_key][ROI_score]['clustering_align_data']:
            if ROI_clustering[ROI_type_key][ROI_score]['clustering_method'] == 'window_peak' or ROI_clustering[ROI_type_key][ROI_score]['clustering_method'] == 'window_stat':
                ax[t1].axvline(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][0] - \
                    ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                    color=figParams['window_lineColor'],linewidth=figParams['imshow_lineWidth']+0.5,alpha=figParams['alpha'])
                ax[t1].axvline(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][1] - \
                    ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                    color=figParams['window_lineColor'],linewidth=figParams['imshow_lineWidth']+0.5,alpha=figParams['alpha'])
                if ttype in ROI_clustering[ROI_type_key][ROI_score]['clustering_ttypes']:
                    ax[t1].text(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][0] - \
                        ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                            0-figParams['ylim_adjust'][0], ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_short_label'],\
                        ha='left',va='bottom',fontsize=figParams['page2_fontSize'],color=figParams['window_lineColor'],alpha=figParams['alpha'])
            elif ROI_clustering[ROI_type_key][ROI_score]['clustering_method'] == 'multi_window_peak':
                for w in range(len(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'])):
                    if ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][0] - \
                        ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data] >= xlim_s[0] and \
                        ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][0] - \
                            ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data] <= xlim_s[1]:
                        ax[t1].axvline(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][0] - \
                            ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                            color=figParams['window_lineColor'],linewidth=figParams['imshow_lineWidth']+0.5,alpha=figParams['alpha'])
                        if ttype in ROI_clustering[ROI_type_key][ROI_score]['clustering_ttypes']:
                            ax[t1].text(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][0] - \
                                ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                                    0-figParams['ylim_adjust'][0], ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_short_label'][w],\
                                ha='left',va='bottom',fontsize=figParams['page2_fontSize'],color=figParams['window_lineColor'],alpha=figParams['alpha'])
                    if ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][1] - \
                        ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data] >= xlim_s[0] and \
                        ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][1] - \
                            ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data] <= xlim_s[1]:
                        ax[t1].axvline(ROI_clustering[ROI_type_key][ROI_score]['clustering_target_window_s'][w][1] - \
                            ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['group']]['traceAlignParams']['all_align_shift_s'][align_data],\
                                color=figParams['window_lineColor'],linewidth=figParams['imshow_lineWidth']+0.5,alpha=figParams['alpha'])
        if 'allSess' in figParams['trialGrouping']:
            ax[t1].text(xlim_s[0], -5, ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+ttype,\
                ha='left',va='bottom',fontsize=figParams['page2_fontSize'],color=(0,0,0))
        elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
            be = figParams['trialGrouping_behavior_epoch']
            ax[t1].text(xlim_s[0], -5, ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+be+" "+ttype,\
                ha='left',va='bottom',fontsize=figParams['page2_fontSize'],color=(0,0,0))
        
        if figParams['filt_trace_byROI']:
            if align_data in figParams['filt_trace_byROI_window'].keys():
                if ROI_type_key in figParams['filt_trace_byROI_window'][align_data].keys():
                    filtLabel = "Filt"+str(figParams['filt_trace_byROI_window'][align_data][ROI_type_key])+" "
                else:
                    filtLabel = ""
            else:
                filtLabel = ""
        else:
            filtLabel = ""
        if t1+1 == len(figParams['export_ttypes']):
            if not figParams['clean']:
                ax[t1].text(xlim_s[1],allClusterData.shape[0]/2,\
                    trace_label_byROI+"\n(n = "+str(total)+" ROIs "+str(len(anmList))+" anm; "+filtLabel+"\n"+figParams['cmap']+" "+\
                        str(np.round(vmin,decimals=2))+"<->"+str(np.round(vmax,decimals=2))+")"+"\n\n\n\n\n",fontsize=figParams['page2_fontSize']-4,rotation=-90,ha='center',va='center')
            else:
                ax[t1].text(xlim_s[1],allClusterData.shape[0]/2,"(n = "+str(total)+" "+str(len(anmList))+" anm)\n"+str(np.round(vmin,decimals=2))+"<->"+str(np.round(vmax,decimals=2))+")"+"\n\n",fontsize=figParams['page2_fontSize']-4,rotation=-90,ha='center',va='center')

        if not figParams['clean']:
            plt.subplots_adjust(wspace=0.01,hspace=0.01)
            # set_dynamic_suptitle(fig,figParams['batchID']+"\n"+"Clust: "+ROI_clustering[ROI_type_key][ROI_score]['detailed_clustering_label_title']+"(n = "+str(total)+" ROIs "+str(len(anmList))+" anm)",\
            #                 margin = 0.01,fontsize=figParams['page2_fontSize'],color=(0,0,0),va='bottom')
            #                 # y=0.95+(1/(ROI_clustering[ROI_type_key][ROI_score]['clustering_nROI']*figParams['ROI_vscalar'][k])),fontsize=figParams['page2_fontSize'],color=(0,0,0))
    
    
    
    if figParams['clean']:
        for t1 in range(len(figParams['export_ttypes'])):
            try:
                pos = ax[t1].get_position()
                # print(f'pos for col {t1} relwidth {pos.width:.3f} width {figsize[0]*pos.width:.3f}in {2.2*figsize[0]*pos.width:.3f}cm relheight {pos.height:.3f} height {figsize[1]*pos.height:.3f}in {2.2*figsize[1]*pos.height:.3f}cm')
            except:
                pass



    return fig,ax
#########################################################################################################
def load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False):

    curr_sessions = list(summaryInfo['anm_details'][anmIdx][ROI_type_key]['sessions'].keys())
    sess = curr_sessions[0]
    anm = summaryInfo['anms'][anmIdx]

    if traceAlignParams:
        orig_align_framerate = traceAlignParams['all_alignment_frame_params'][ROI_type_key][anm][sess]['orig_align_framerate']
        align_framerate = traceAlignParams['all_alignment_frame_params'][ROI_type_key][anm][sess]['align_framerate']
        orig_align_shift_fr = traceAlignParams['all_alignment_frame_params'][ROI_type_key][anm][sess]['orig_align_shift_fr']
        align_reduce_factor = traceAlignParams['align_reduce_factor'][ROI_type_key]
        align_shift_s = traceAlignParams['align_shift_s']
        align_shift_fr = traceAlignParams['all_alignment_frame_params'][ROI_type_key][anm][sess]['align_shift_fr']
    else:
        orig_align_framerate = float(copy.deepcopy(summaryInfo['anm_details'][anmIdx][ROI_type_key]['sessions'][sess]['fr']))
        align_framerate = float(copy.deepcopy(summaryInfo['anm_details'][anmIdx][ROI_type_key]['sessions'][sess]['fr']))
        orig_align_shift_fr = 0
        align_reduce_factor = 1
        align_shift_s = 0
        align_shift_fr = 0

    fr_adjust = 1
    trial=0
    finding_clean=True
    orig_StartFrame=0
    orig_SampleFrame=np.nan
    orig_SampleFrames=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['GeneralParams']['CuePeriod_frames'])
    orig_DelayFrame=np.nan
    orig_DelayFrames=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['GeneralParams']['DelayPeriod_frames'])
    orig_ResponseFrame=np.nan
    orig_ResponseFrames=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['GeneralParams']['ResponsePeriod_frames'])
    while finding_clean:
        if int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['ImagingMatch']):
            if int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['cleanTrial']) in \
                summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['CleanTrialClass']:
                finding_clean=False
                for e in range(len(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'])):
                    if "TrigTrialStart" in summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['Label'][e]:
                        orig_StartFrame=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'][e] - fr_adjust)
                    if "SamplePeriod" in summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['Label'][e]:
                        orig_SampleFrame=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'][e] - fr_adjust)
                    if "DelayPeriod" in summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['Label'][e]:
                        orig_DelayFrame=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'][e] - fr_adjust)
                    if "ResponseCue" in summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['Label'][e]:
                        orig_ResponseFrame=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'][e] - fr_adjust)
                    if "AnswerPeriod" in summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['Label'][e]:
                        orig_AnswerFrame=int(summaryInfo[ROI_type_key+'_pooledData'][anmIdx]['sessions'][sess]['Behavior']['Trials'][trial]['EventStruct']['event_frame'][e] - fr_adjust)
            else:
                trial+=1
        else:
            trial+=1
    
    # orig_SampleFrames = orig_DelayFrame-orig_SampleFrame-fr_adjust
    # orig_DelayFrames = orig_ResponseFrame-orig_DelayFrame-fr_adjust
    precision = 3

    StartTime=np.round(orig_StartFrame/orig_align_framerate,decimals=precision).astype('float32')
    SampleTime=np.round(orig_SampleFrame/orig_align_framerate,decimals=precision).astype('float32')
    SampleTimes=np.round(orig_SampleFrames/orig_align_framerate,decimals=precision).astype('float32')
    DelayTime=np.round(orig_DelayFrame/orig_align_framerate,decimals=precision).astype('float32')
    DelayTimes=np.round(orig_DelayFrames/orig_align_framerate,decimals=precision).astype('float32')
    ResponseTime=np.round(orig_ResponseFrame/orig_align_framerate,decimals=precision).astype('float32')
    ResponseTimes=np.round(orig_ResponseFrames/orig_align_framerate,decimals=precision).astype('float32')
    AnswerTime=np.round(orig_AnswerFrame/orig_align_framerate,decimals=precision).astype('float32')

    if fix_overlaps:
        if SampleTime - StartTime + SampleTimes > DelayTime - StartTime:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing SampleTimes overlap "+str(SampleTimes))
        elif  (DelayTime - StartTime) - (SampleTime - StartTime + SampleTimes) > 0.02:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing SampleTimes Gap "+str(SampleTimes))
        if DelayTime - StartTime + DelayTimes > ResponseTime - StartTime:
            DelayTimes = ResponseTime - DelayTime
            if verbose:
                print("Fixing DelayTimes overlap "+str(DelayTimes))
        elif  (ResponseTime - StartTime) - (DelayTime - StartTime + DelayTimes) > 0.02:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing DelayTimes Gap "+str(SampleTimes))

    if verbose:
        print("Sample:   "+str([SampleTime-StartTime,SampleTime-StartTime+SampleTimes]))
        print("Delay:    "+str([DelayTime-StartTime,DelayTime-StartTime+DelayTimes]))
        print("Response: "+str([ResponseTime-StartTime,ResponseTime-StartTime+ResponseTimes]))


    if align_reduce_factor != 1:
        StartFrame = int(np.floor(orig_StartFrame/align_reduce_factor))
        if not np.isnan(orig_SampleFrame):
            SampleFrame = int(np.floor(orig_SampleFrame/align_reduce_factor))
        SampleFrames = int(np.floor(orig_SampleFrames/align_reduce_factor))
        if not np.isnan(orig_DelayFrame):
            DelayFrame = int(np.floor(orig_DelayFrame/align_reduce_factor))
        DelayFrames = int(np.floor(orig_DelayFrames/align_reduce_factor))
        if not np.isnan(orig_ResponseFrame):
            ResponseFrame = int(np.floor(orig_ResponseFrame/align_reduce_factor))
        ResponseFrames = int(np.floor(orig_ResponseFrames/align_reduce_factor))
        if not np.isnan(orig_AnswerFrame):
            AnswerFrame = int(np.floor(orig_AnswerFrame/align_reduce_factor))
    else:
        # StartFrame = copy.deepcopy(orig_StartFrame) + fr_adjust
        StartFrame = copy.deepcopy(orig_StartFrame)
        SampleFrame = copy.deepcopy(orig_SampleFrame)
        SampleFrames = copy.deepcopy(orig_SampleFrames)
        DelayFrame = copy.deepcopy(orig_DelayFrame)
        DelayFrames = copy.deepcopy(orig_DelayFrames)
        ResponseFrame = copy.deepcopy(orig_ResponseFrame)
        ResponseFrames = copy.deepcopy(orig_ResponseFrames)
        AnswerFrame = copy.deepcopy(orig_AnswerFrame)

    trial_structure_times = {}
    trial_structure_times['StartFrame'] = StartFrame
    trial_structure_times['SampleFrame'] = SampleFrame
    trial_structure_times['SampleFrames'] = SampleFrames
    trial_structure_times['DelayFrame'] = DelayFrame
    trial_structure_times['DelayFrames'] = DelayFrames
    trial_structure_times['ResponseFrame'] = ResponseFrame
    trial_structure_times['ResponseFrames'] = ResponseFrames
    trial_structure_times['AnswerFrame'] = AnswerFrame
    trial_structure_times['StartTime'] = StartTime
    trial_structure_times['SampleTime'] = SampleTime
    trial_structure_times['SampleTimes'] = SampleTimes
    trial_structure_times['DelayTime'] = DelayTime
    trial_structure_times['DelayTimes'] = DelayTimes
    trial_structure_times['ResponseTime'] = ResponseTime
    trial_structure_times['ResponseTimes'] = ResponseTimes
    trial_structure_times['AnswerTime'] = AnswerTime
    trial_structure_times['orig_align_framerate'] = orig_align_framerate
    trial_structure_times['orig_align_shift_fr'] = orig_align_shift_fr
    trial_structure_times['align_reduce_factor'] = align_reduce_factor
    trial_structure_times['align_shift_s'] = align_shift_s
    trial_structure_times['align_shift_fr'] = align_shift_fr
    trial_structure_times['align_framerate'] = align_framerate
    trial_structure_times['fr_adjust'] = fr_adjust
    
    return trial_structure_times

def add_trial_structure_features(xMode,ax,figParams,trial_structure_times,ROI_type_key,ROI_score,group,align_data,ttype,ttypes,\
    vertLinesOn,horzLinesOn,boxesOn,textOn,lineColor,vertLineWidth,vertLineAlpha,horzLineWidth,horzLineAlpha,fontSize,plotminVal,plotmaxVal,\
        horzCueLineOn=False,horzCueScalar=0.02,startMark = True, delayMark = False, goMark = True, horzLoc = 'top'):        
    
    StartFrame = trial_structure_times['StartFrame']
    SampleFrame = trial_structure_times['SampleFrame']
    SampleFrames = trial_structure_times['SampleFrames']
    DelayFrame = trial_structure_times['DelayFrame']
    DelayFrames = trial_structure_times['DelayFrames']
    ResponseFrame = trial_structure_times['ResponseFrame']
    ResponseFrames = trial_structure_times['ResponseFrames']
    AnswerFrame = trial_structure_times['AnswerFrame']
    StartTime = trial_structure_times['StartTime']
    SampleTime = trial_structure_times['SampleTime']
    SampleTimes = trial_structure_times['SampleTimes']
    DelayTime = trial_structure_times['DelayTime']
    DelayTimes = trial_structure_times['DelayTimes']
    ResponseTime = trial_structure_times['ResponseTime']
    ResponseTimes = trial_structure_times['ResponseTimes']
    AnswerTime = trial_structure_times['AnswerTime']
    align_shift_s = trial_structure_times['align_shift_s']
    align_shift_fr = trial_structure_times['align_shift_fr']

    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    if plotmaxVal < plotminVal:
        vertLine = ylim
    else:
        vertLine = [plotminVal,plotmaxVal]

    # print(f'plotminVal,plotmaxVal: {plotminVal},{plotmaxVal}')
    # print(f'ylim: {ylim}')
    # print(f'vertLine: {vertLine}')
    trialType = ''
    if "R" in ttype and not "L" in ttype:
        trialType = 'R'
        for t in ttypes:
            if "L" in t:
                trialType = ''
    elif "L" in ttype and not "R" in ttype:
        trialType = 'L'
        for t in ttypes:
            if "R" in t:
                trialType = ''
    #####################################################################################################
    if boxesOn:
        if 'trialStart' in group:
            if 'time' in xMode:
                SampleRectangle=Rectangle(((SampleTime-StartTime),plotminVal),\
                                            SampleTimes,np.abs(plotmaxVal-plotminVal))
                DelayRectangle=Rectangle(((DelayTime-StartTime),plotminVal),\
                                            DelayTimes,np.abs(plotmaxVal-plotminVal))
                ResponseRectangle=Rectangle(((ResponseTime-StartTime),plotminVal),\
                                            ResponseTimes,np.abs(plotmaxVal-plotminVal))
            else:
                SampleRectangle=Rectangle(((SampleFrame-StartFrame),plotminVal),\
                                            SampleFrames,np.abs(plotmaxVal-plotminVal))
                DelayRectangle=Rectangle(((DelayFrame-StartFrame),plotminVal),\
                                            DelayFrames,np.abs(plotmaxVal-plotminVal))
                ResponseRectangle=Rectangle(((ResponseFrame-StartFrame),plotminVal),\
                                            ResponseFrames,np.abs(plotmaxVal-plotminVal))

            patches = []
            colors = []
            if startMark:
                patches.append(SampleRectangle)
                if trialType == 'R':
                    colors.append(figParams['facecolors_right'][0])
                elif trialType == 'L':
                    colors.append(figParams['facecolors_left'][0])
                else:
                    colors.append(figParams['facecolors_generic'][0])
            if delayMark:
                patches.append(DelayRectangle)
                if trialType == 'R':
                    colors.append(figParams['facecolors_right'][1])
                elif trialType == 'L':
                    colors.append(figParams['facecolors_left'][1])
                else:
                    colors.append(figParams['facecolors_generic'][1])
            if goMark:
                patches.append(ResponseRectangle)
                if trialType == 'R':
                    colors.append(figParams['facecolors_right'][2])
                elif trialType == 'L':
                    colors.append(figParams['facecolors_left'][2])
                else:
                    colors.append(figParams['facecolors_generic'][2])

            # if trialType == 'R':
            #     GeneralPatches = PatchCollection(patches, facecolor=figParams['facecolors_right'][colorIdx], \
            #         alpha=figParams['patch_alpha'], edgecolor = figParams['edgecolor'])
            # elif trialType == 'L':
            #     GeneralPatches = PatchCollection(patches, facecolor=figParams['facecolors_left'][colorIdx], \
            #         alpha=figParams['patch_alpha'], edgecolor = figParams['edgecolor'])
            # else:
            #     GeneralPatches = PatchCollection(patches, facecolor=figParams['facecolors_generic'][colorIdx], \
            #         alpha=figParams['patch_alpha'], edgecolor = figParams['edgecolor'])

            GeneralPatches = PatchCollection(patches, facecolor=colors, \
                alpha=figParams['patch_alpha'], edgecolor = figParams['edgecolor'])

            ax.add_collection(GeneralPatches);
    if vertLinesOn:
        if 'time' in xMode:
            if 'trialStart' in group:
                if startMark:
                    # ax.axvline(SampleTime-StartTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(SampleTime-StartTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)

                if delayMark:
                    # ax.axvline(DelayTime-StartTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(DelayTime-StartTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                if goMark:
                    # ax.axvline(ResponseTime-StartTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(ResponseTime-StartTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    # ax.axvline(AnswerTime-StartTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
            elif 'goCue' in group:
                if startMark:
                    # ax.axvline(SampleTime-ResponseTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(SampleTime-ResponseTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                        
                if delayMark:
                    # ax.axvline(DelayTime - ResponseTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(DelayTime - ResponseTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                if goMark:
                    # ax.axvline(ResponseTime - ResponseTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(ResponseTime - ResponseTime),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    # ax.axvline(AnswerTime-StartTime,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
        else:
            if 'trialStart' in group or 'goCue' in group:
                if startMark:
                    # ax.axvline(SampleFrame-StartFrame - align_shift_fr,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(SampleFrame-StartFrame - align_shift_fr),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                if delayMark:
                    # ax.axvline(DelayFrame-StartFrame - align_shift_fr,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(DelayFrame-StartFrame - align_shift_fr),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                if goMark:
                    # ax.axvline(ResponseFrame-StartFrame - align_shift_fr,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    ax.plot(np.ones(2)*(ResponseFrame-StartFrame - align_shift_fr),vertLine,\
                        color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                    # ax.axvline(AnswerFrame-StartFrame - align_shift_fr,\
                    #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
    if boxesOn or vertLinesOn:
        if not 'goCue' in group:
            if 'time' in xMode:
                # ax.axvline(0, color=lineColor,linewidth=vertLineWidth, alpha=vertLineAlpha)
                ax.plot(np.ones(2)*(0),vertLine,\
                    color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                # ax.axvline(-align_shift_fr,\
                #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
            else:
                # ax.axvline(-align_shift_fr,\
                #     color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
                ax.plot(np.ones(2)*(- align_shift_fr),vertLine,\
                    color=lineColor, linewidth=vertLineWidth,alpha=vertLineAlpha)
    
    
    if horzLoc == 'top':
        if plotmaxVal < plotminVal:
            horzCueSize = np.abs(ylim[1]-ylim[0]) * horzCueScalar
            horzCueShift = horzCueSize
            # print(horzCueSize)
        else:
            horzCueSize = np.abs(plotmaxVal-plotminVal) * horzCueScalar
            horzCueShift = -horzCueSize
    else:
        if plotmaxVal < plotminVal:
            horzCueSize = np.abs(ylim[1]-ylim[0]) * horzCueScalar
            horzCueShift = -horzCueSize
        else:
            horzCueSize = np.abs(plotmaxVal-plotminVal) * horzCueScalar
            horzCueShift = horzCueSize

    if horzCueLineOn and ('trialStart' in group or 'goCue' in group):
        mixedCue = False
        if trialType == 'L':
            cue = 'left'
        elif trialType == 'R':
            cue = 'right'
        else:
            cue = 'left'
            mixedCue = True
        sampleTrace_ydata,sampleTrace_xdata,sampleTrace_env = generate_pulsed_tone(
            carrier_hz=figParams['cueInfo'][cue]['sim']['carrier_hz'] ,       # set your tone frequency here
            pulse_rate_hz=figParams['cueInfo'][cue]['sim']['pulse_rate_hz'] ,       # 4 pulses per second
            on_duration_s=figParams['cueInfo'][cue]['sim']['on_duration_s'] ,   # ON duration per pulse
            duration_s=float(SampleTimes) ,          # total length
            sr=figParams['cueInfo'][cue]['sim']['sr'] ,
            amplitude=figParams['cueInfo'][cue]['sim']['amplitude'] ,
            fade_ms=figParams['cueInfo'][cue]['sim']['fade_ms'] ,             # tweak for softer/harder edges
            return_envelope=True)
        cue = 'go'
        responseTrace_ydata,responseTrace_xdata,responseTrace_env = generate_pulsed_tone(
            carrier_hz=figParams['cueInfo'][cue]['sim']['carrier_hz'] ,       # set your tone frequency here
            pulse_rate_hz=figParams['cueInfo'][cue]['sim']['pulse_rate_hz'] ,       # 4 pulses per second
            on_duration_s=figParams['cueInfo'][cue]['sim']['on_duration_s'] ,      # ON duration per pulse
            duration_s=float(ResponseTimes) ,          # total length
            sr=figParams['cueInfo'][cue]['sim']['sr'] ,
            amplitude=figParams['cueInfo'][cue]['sim']['amplitude'] ,
            fade_ms=figParams['cueInfo'][cue]['sim']['fade_ms'] ,             # tweak for softer/harder edges
            return_envelope=True)
        responseTrace_ydata[responseTrace_env<=0] = np.nan
        horzCueLineWidth = 0.5
        # horzCueAdjust = 0.0
        if horzLoc == 'top':
            if plotmaxVal < plotminVal:
                horzCueLineMax = np.nanmin([np.nanmin(sampleTrace_ydata*horzCueSize+plotmaxVal),np.nanmin(responseTrace_ydata*horzCueSize+plotmaxVal)])
                ax.set_ylim([ylim[0],horzCueLineMax*1.01])
            else:
                horzCueLineMax = np.nanmax([np.nanmax(sampleTrace_ydata*horzCueSize+plotmaxVal),np.nanmax(responseTrace_ydata*horzCueSize+plotmaxVal)])
                plotmaxVal = plotmaxVal - (horzCueLineMax - plotmaxVal)
            horzPlotLoc = plotmaxVal
            va = 'top'
        else:
            if plotmaxVal < plotminVal:
                horzCueLineMin = np.nanmax([np.nanmax(sampleTrace_ydata*horzCueSize+plotminVal),np.nanmax(responseTrace_ydata*horzCueSize+plotminVal)])
                ax.set_ylim([horzCueLineMin*0.99,ylim[1]])
            else:
                horzCueLineMin = np.nanmin([np.nanmin(sampleTrace_ydata*horzCueSize+plotminVal),np.nanmin(responseTrace_ydata*horzCueSize+plotminVal)])
                plotminVal = plotminVal + (plotminVal - horzCueLineMin)
            horzPlotLoc = plotminVal
            va = 'bottom'
        # print(f'horzLoc {horzLoc} horzPlotLoc {horzPlotLoc}')
    if horzLinesOn:
        if 'time' in xMode:
            if 'trialStart' in group:
                if trialType == 'R':
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(sampleTrace_xdata+(SampleTime-StartTime),sampleTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_right'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-StartTime,SampleTime-StartTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_right'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if delayMark:
                        ax.plot([DelayTime-StartTime,DelayTime-StartTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_right'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata+(ResponseTime-StartTime),responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_right'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-StartTime,ResponseTime-StartTime+ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_right'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                elif trialType == 'L':
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(sampleTrace_xdata+(SampleTime-StartTime),sampleTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_left'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-StartTime,SampleTime-StartTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_left'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if delayMark:
                        ax.plot([DelayTime-StartTime,DelayTime-StartTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_left'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata+(ResponseTime-StartTime),responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_left'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-StartTime,ResponseTime-StartTime+ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_left'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                else:
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(np.insert(sampleTrace_xdata, 0, -1*sampleTrace_xdata[1])+(SampleTime-StartTime),np.insert(sampleTrace_env, 0, 0)*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_generic'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-StartTime,SampleTime-StartTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_generic'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if delayMark:
                        ax.plot([DelayTime-StartTime,DelayTime-StartTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_generic'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata+(ResponseTime-StartTime),responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_generic'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-StartTime,ResponseTime-StartTime+ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_generic'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                            
            elif 'goCue' in group:
                if trialType == 'R':
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(sampleTrace_xdata+(SampleTime-ResponseTime),sampleTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_right'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-ResponseTime,SampleTime-ResponseTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_right'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if delayMark:
                        ax.plot([DelayTime-ResponseTime,DelayTime-ResponseTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_right'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata,responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_right'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-ResponseTime,ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_right'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')

                elif trialType == 'L':
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(sampleTrace_xdata+(SampleTime-ResponseTime),sampleTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_left'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-ResponseTime,SampleTime-ResponseTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_left'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')

                    if delayMark:
                        ax.plot([DelayTime-ResponseTime,DelayTime-ResponseTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_left'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata,responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_left'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-ResponseTime,ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_left'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')

                else:
                    if startMark:
                        if horzCueLineOn:
                            ax.plot(np.insert(sampleTrace_xdata, 0, -1*sampleTrace_xdata[1])+(SampleTime-ResponseTime),np.insert(sampleTrace_env, 0, 0)*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_generic'][0],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([SampleTime-ResponseTime,SampleTime-ResponseTime+SampleTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_generic'][0],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')

                    if delayMark:
                        ax.plot([DelayTime-ResponseTime,DelayTime-ResponseTime+DelayTimes],[horzPlotLoc,horzPlotLoc],\
                            linestyle='-',color=figParams['facecolors_generic'][1],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    if goMark:
                        if horzCueLineOn:
                            ax.plot(responseTrace_xdata,responseTrace_ydata*horzCueSize+horzPlotLoc+horzCueShift,\
                                linestyle='-',color=figParams['facecolors_generic'][2],linewidth=horzCueLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                        else:
                            ax.plot([ResponseTime-ResponseTime,ResponseTimes],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=figParams['facecolors_generic'][2],linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
        else:
            if 'trialStart' in group:
                if trialType == 'R':
                    ax.plot(np.array([SampleFrame-StartFrame,SampleFrame-StartFrame+SampleFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_right'][0],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([DelayFrame-StartFrame,DelayFrame-StartFrame+DelayFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_right'][1],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([ResponseFrame-StartFrame,ResponseFrame-StartFrame+ResponseFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_right'][2],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                elif trialType == 'L':
                    ax.plot(np.array([SampleFrame-StartFrame,SampleFrame-StartFrame+SampleFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_left'][0],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([DelayFrame-StartFrame,DelayFrame-StartFrame+DelayFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_left'][1],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([ResponseFrame-StartFrame,ResponseFrame-StartFrame+ResponseFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_left'][2],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                else:
                    ax.plot(np.array([SampleFrame-StartFrame,SampleFrame-StartFrame+SampleFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_generic'][0],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([DelayFrame-StartFrame,DelayFrame-StartFrame+DelayFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_generic'][1],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')
                    ax.plot(np.array([ResponseFrame-StartFrame,ResponseFrame-StartFrame+ResponseFrames]) - align_shift_fr,\
                        [horzPlotLoc,horzPlotLoc],linestyle='-',color=figParams['facecolors_generic'][2],\
                            linewidth=horzLineWidth,alpha=horzLineAlpha,solid_capstyle='butt')            
    if textOn:
        if 'time' in xMode:
            if ('cTimes' in group or 'maxExtCL' in group):
                ax.text(0,horzPlotLoc,"Contact",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
            elif 'goCue' in group:
                ax.text(0,horzPlotLoc,"Go Cue",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
            elif 'postContLick1' in group:
                ax.text(0,horzPlotLoc,"2nd Lick",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
        else:
            if ('cTimes' in group or 'maxExtCL' in group):
                ax.text(-align_shift_fr,horzPlotLoc,"Contact",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
            elif 'goCue' in group:
                ax.text(-align_shift_fr,horzPlotLoc,"Go Cue",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
            elif 'postContLick1' in group:
                ax.text(-align_shift_fr,horzPlotLoc,"2nd Lick",ha='left',va=va,fontsize=fontSize,color=(0,0,0))
    
    return ax


def simple_load_alignment_xdata(traceAlignParams,ROI_type_key):
    nFr = copy.deepcopy(traceAlignParams['align_length_fr'][ROI_type_key])
    fr  = copy.deepcopy(traceAlignParams['align_framerate'][ROI_type_key])
    align_shift_s = copy.deepcopy(traceAlignParams['align_shift_s'])
    align_shift_fr = copy.deepcopy(traceAlignParams['align_shift_fr'][ROI_type_key])
    align_reduce_factor = copy.deepcopy(traceAlignParams['align_reduce_factor'][ROI_type_key])
    xdata_s = copy.deepcopy(traceAlignParams['align_xdata_s'][ROI_type_key])
    xdata_fr = copy.deepcopy(traceAlignParams['align_xdata_fr'][ROI_type_key])
    return xdata_s,xdata_fr,nFr,fr,align_shift_s,align_shift_fr,align_reduce_factor

def load_alignment_xdata(traceAlignParams,ROI_type_key,group,zoom,figParams):


    nFr = copy.deepcopy(traceAlignParams['align_length_fr'][ROI_type_key])
    framerate  = copy.deepcopy(traceAlignParams['align_framerate'][ROI_type_key])
    align_shift_s = copy.deepcopy(traceAlignParams['align_shift_s'])
    align_shift_fr = copy.deepcopy(traceAlignParams['align_shift_fr'][ROI_type_key])
    align_reduce_factor = copy.deepcopy(traceAlignParams['align_reduce_factor'][ROI_type_key])
    xdata_s = copy.deepcopy(traceAlignParams['align_xdata_s'][ROI_type_key])
    xdata_fr = copy.deepcopy(traceAlignParams['align_xdata_fr'][ROI_type_key])
    dataCorrection = (1/framerate)/2

    if zoom:
        if 'xticks_zoom' in figParams:
            xticks_s = copy.deepcopy(figParams['xticks_zoom'][group])
        else:
            print("<<WARNING>> Please provide 'xticks_zoom' in figParams for zoomed plots.")
            xticks_s = []
    else:
        if 'xticks' in figParams:
            xticks_s = copy.deepcopy(figParams['xticks'][group])
        else:
            print("<<WARNING>> Please provide 'xticks' in figParams for non-zoomed plots.")
            xticks_s = []



    # print(dataCorrection)

    # xticks_fr = []
    # xticks_idx = []
    # xticks_s = []
    # xticks = []
    if zoom:
        # xlim_s = copy.deepcopy(figParams['xlim_zoom'][group])
        xlim_idx = [np.argmin(np.abs(xdata_s-figParams['xlim_zoom'][group][0])),\
                np.argmin(np.abs(xdata_s-figParams['xlim_zoom'][group][1]))]
        # xlim_fr = [np.argmin(np.abs(xdata_s-figParams['xlim_zoom'][group][0])),\
        #         np.argmin(np.abs(xdata_s-figParams['xlim_zoom'][group][1]))]          
        
        for x in figParams['xticks_zoom'][group]:
            # xticks.append(str(x))
            tempIdx = np.argmin(np.abs(xdata_s-x))
            # xticks_idx.append(tempIdx)
            # xticks_s.append(xdata_s[tempIdx])
            # xticks_fr.append(xdata_fr[tempIdx])
    else:
        # xlim_s = [np.min(xdata_s),np.max(xdata_s)]
        # xlim_fr = [-0.5,nFr+0.5]
        for x in figParams['xticks'][group]:
            # xticks.append(str(x))
            tempIdx = np.argmin(np.abs(xdata_s-x))
            # xticks_idx.append(tempIdx)
            # xticks_s.append(xdata_s[tempIdx])
            # xticks_fr.append(xdata_fr[tempIdx])
        xlim_idx = [0,len(xdata_s)-1]
    plotIdxs = np.arange(xlim_idx[0],xlim_idx[1]+1).astype(int)

    # xlim_s = [np.round(xdata_s[xlim_idx[0]],decimals=2),np.round(xdata_s[xlim_idx[1]],decimals=2)]
    # xlim_s[0] = xlim_s[0] - np.round(dataCorrection,decimals=2)
    # xlim_s[1] = xlim_s[1] + np.round(dataCorrection,decimals=2)
    xlim_s = [xdata_s[plotIdxs[0]]-dataCorrection,xdata_s[plotIdxs[-1]]+dataCorrection]
    # if xlim_s[1] - xlim_s[0] <= 4:
    #     for i in range(len(xticks_s)):
    #         xticks_s[i] = np.round(xticks_s[i],decimals=1)
    # else:
    #     for i in range(len(xticks_s)):
    #         xticks_s[i] = int(np.round(xticks_s[i],decimals=0))
    # for i in range(len(xticks_fr)):
    #     xticks_fr[i] = int(np.round(xticks_fr[i],decimals=0))
    
    # print('align_shift_s '+str(align_shift_s)+' xlim_s '+str(xlim_s)+' align_shift_fr '+str(align_shift_fr)+' xlim_fr '+str(xlim_fr))
    # return xdata_s,xlim_s,xticks_s,xdata_fr,xlim_fr,xticks_fr,xticks_idx,plotIdxs
    return xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate

def xdata_correction(xdata,figParams,group):
    if not 'xdata_correction' in figParams.keys():
        figParams['xdata_correction'] = {}
    if not group in figParams['xdata_correction']:
        figParams['xdata_correction'][group] = False
    if figParams['xdata_correction'][group]:
        xdata1 = copy.deepcopy(xdata)
        xdata1[xdata1>0] = np.nan
        xadjust = xdata[np.nanargmin(np.abs(xdata1))]
        xdata = xdata - xadjust
    else:
        pass
    return xdata

def extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,group,trialGrouping,align_data,be,ttype,r,fillNaNs = False):
    if 'allSess' in trialGrouping:
        roiTrace = copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:])
    if 'byBehavEpoch' in trialGrouping or 'byPrePost' in trialGrouping or 'deltaPrePost' in trialGrouping:
        roiTrace = copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][figParams['plot_data_byROI']][figParams['plot_stat_byROI']][r,:])
    if np.any(np.isfinite(roiTrace)):
        pass
    else:
        if fillNaNs:
            roiTrace[np.isnan(roiTrace)] = 0
    
    if 'filt_trace_byROI' in figParams and 'filt_trace_byROI_window' in figParams:
        if align_data in figParams['filt_trace_byROI_window']:
            if figParams['filt_trace_byROI_window'][align_data]:
                if ROI_type_key in figParams['filt_trace_byROI_window'][align_data]:
                    # print(figParams['filt_trace_byROI_window'][align_data][ROI_type_key])
                    roiTrace = running_mean_convolve(roiTrace,figParams['filt_trace_byROI_window'][align_data][ROI_type_key])
    return roiTrace

def extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,group,align_data,trialGrouping,clustered_data,ttype,be,plot_stat,plot_data,plot_data_errors,plot_data_error_pos,plot_data_error_neg,shift,invertData=False):
    mainTrace = np.ones((1,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][group][align_data]),dtype='float32')*np.nan
    upperTrace = np.ones((1,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][group][align_data]),dtype='float32')*np.nan
    lowerTrace = np.ones((1,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][group][align_data]),dtype='float32')*np.nan
    if invertData:
        adjust =-1
    else:
        adjust=1

    if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
        if 'allSess' in trialGrouping:
            mainTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift
            if plot_data_errors:
                if 'boot_CI' in plot_data_error_pos:
                    upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_pos])*adjust-shift
                    lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_neg])*adjust-shift
                else:
                    upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data]*adjust+\
                                ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_pos])-shift
                    lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data]*adjust-\
                                ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_neg])-shift
            else:
                upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift
                lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift
        elif 'byBehavEpoch' in trialGrouping or 'byPrePost' in trialGrouping or 'deltaPrePost' in trialGrouping:
            mainTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift
            if plot_data_errors:
                if 'boot_CI' in plot_data_error_pos:
                    upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_pos])*adjust-shift
                    lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_neg])*adjust-shift
                    # print(str(np.min(lowerTrace))+" - "+str(np.min(mainTrace))+" - "+str(np.min(upperTrace)))
                    # print(str(np.max(lowerTrace))+" - "+str(np.max(mainTrace))+" - "+str(np.max(upperTrace)))
                else:
                    upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data]*adjust+\
                                ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_pos])-shift
                    lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data]*adjust-\
                                ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data_error_neg])-shift
            else:
                upperTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift
                lowerTrace = (ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][group][align_data][trialGrouping][be][ttype]['cluster_summaryStats'][clustered_data][plot_stat][plot_data])*adjust-shift


    return mainTrace,upperTrace,lowerTrace

def running_mean_convolve(x, N):
    return np.convolve(x, np.ones(N) / float(N), 'same')

def load_clustering_trial_structure_times(ROI_clustering,ROI_type_key,ROI_score,group,align_data,fix_overlaps=True,verbose = False):

    StartFrame = ROI_clustering[ROI_type_key][ROI_score]['all_StartFrame'][group][align_data]
    SampleFrame = ROI_clustering[ROI_type_key][ROI_score]['all_SampleFrame'][group][align_data]
    SampleFrames = ROI_clustering[ROI_type_key][ROI_score]['all_SampleFrames'][group][align_data]
    DelayFrame = ROI_clustering[ROI_type_key][ROI_score]['all_DelayFrame'][group][align_data]
    DelayFrames = ROI_clustering[ROI_type_key][ROI_score]['all_DelayFrames'][group][align_data]
    ResponseFrame = ROI_clustering[ROI_type_key][ROI_score]['all_ResponseFrame'][group][align_data]
    ResponseFrames = ROI_clustering[ROI_type_key][ROI_score]['all_ResponseFrames'][group][align_data]
    AnswerFrame = ROI_clustering[ROI_type_key][ROI_score]['all_AnswerFrame'][group][align_data]
    StartTime = ROI_clustering[ROI_type_key][ROI_score]['all_StartTime'][group][align_data]
    SampleTime = ROI_clustering[ROI_type_key][ROI_score]['all_SampleTime'][group][align_data]
    SampleTimes = ROI_clustering[ROI_type_key][ROI_score]['all_SampleTimes'][group][align_data]
    DelayTime = ROI_clustering[ROI_type_key][ROI_score]['all_DelayTime'][group][align_data]
    DelayTimes = ROI_clustering[ROI_type_key][ROI_score]['all_DelayTimes'][group][align_data]
    ResponseTime = ROI_clustering[ROI_type_key][ROI_score]['all_ResponseTime'][group][align_data]
    ResponseTimes = ROI_clustering[ROI_type_key][ROI_score]['all_ResponseTimes'][group][align_data]
    AnswerTime = ROI_clustering[ROI_type_key][ROI_score]['all_AnswerTime'][group][align_data]
    
    align_shift_s = ROI_clustering[ROI_type_key][ROI_score]['allParams'][group]['traceAlignParams']['all_align_shift_s'][align_data]
    align_shift_fr = ROI_clustering[ROI_type_key][ROI_score]['allParams'][group]['traceAlignParams']['all_align_shift_fr'][align_data]

    trace_nFr = ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][group][align_data]
    trace_framerate =  ROI_clustering[ROI_type_key][ROI_score]['all_trace_framerate'][group][align_data]
    orig_trace_nFr = ROI_clustering[ROI_type_key][ROI_score]['all_orig_trace_nFr'][group][align_data]
    orig_trace_framerate = ROI_clustering[ROI_type_key][ROI_score]['all_orig_trace_framerate'][group][align_data]



    if fix_overlaps:
        if SampleTime - StartTime + SampleTimes > DelayTime - StartTime:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing SampleTimes overlap "+str(SampleTimes))
        elif  (DelayTime - StartTime) - (SampleTime - StartTime + SampleTimes) > 0.02:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing SampleTimes Gap "+str(SampleTimes))
        if DelayTime - StartTime + DelayTimes > ResponseTime - StartTime:
            DelayTimes = ResponseTime - DelayTime
            if verbose:
                print("Fixing DelayTimes overlap "+str(DelayTimes))
        elif  (ResponseTime - StartTime) - (DelayTime - StartTime + DelayTimes) > 0.02:
            SampleTimes = DelayTime - SampleTime
            if verbose:
                print("Fixing DelayTimes Gap "+str(SampleTimes))

    if verbose:
        print("Sample:   "+str([SampleTime-StartTime,SampleTime-StartTime+SampleTimes]))
        print("Delay:    "+str([DelayTime-StartTime,DelayTime-StartTime+DelayTimes]))
        print("Response: "+str([ResponseTime-StartTime,ResponseTime-StartTime+ResponseTimes]))
        
    trial_structure_times = {}
    trial_structure_times['StartFrame'] = StartFrame
    trial_structure_times['SampleFrame'] = SampleFrame
    trial_structure_times['SampleFrames'] = SampleFrames
    trial_structure_times['DelayFrame'] = DelayFrame
    trial_structure_times['DelayFrames'] = DelayFrames
    trial_structure_times['ResponseFrame'] = ResponseFrame
    trial_structure_times['ResponseFrames'] = ResponseFrames
    trial_structure_times['AnswerFrame'] = AnswerFrame
    trial_structure_times['StartTime'] = StartTime
    trial_structure_times['SampleTime'] = SampleTime
    trial_structure_times['SampleTimes'] = SampleTimes
    trial_structure_times['DelayTime'] = DelayTime
    trial_structure_times['DelayTimes'] = DelayTimes
    trial_structure_times['ResponseTime'] = ResponseTime
    trial_structure_times['ResponseTimes'] = ResponseTimes
    trial_structure_times['AnswerTime'] = AnswerTime
    trial_structure_times['align_shift_s'] = align_shift_s
    trial_structure_times['align_shift_fr'] = align_shift_fr
    trial_structure_times['trace_nFr'] = trace_nFr
    trial_structure_times['trace_framerate'] = trace_framerate
    trial_structure_times['orig_trace_nFr'] = orig_trace_nFr
    trial_structure_times['orig_trace_framerate'] = orig_trace_framerate

    return trial_structure_times


#################################################################################################
def generate_pulsed_tone(
    carrier_hz: float = 1000.0,
    pulse_rate_hz: float = 4.0,   # pulses per second
    on_duration_s: float = 0.15,  # ON time per pulse
    duration_s: float = 5.0,      # total signal duration
    sr: int = 44100,              # sample rate
    amplitude: float = 0.8,       # 0..1 (will be clipped to avoid >1)
    fade_ms: float = 5.0,         # attack/release crossfade to avoid clicks
    start_phase: float = 0.0,     # radians
    return_envelope: bool = False
):
    

    """
    Generate a pulsed sine tone:
    - Pulses repeat at pulse_rate_hz, each ON for on_duration_s (OFF is whatever remains of the period).
    - The sine runs only during the ON portions (via amplitude envelope).
    - Adds short attack/release fades to reduce clicks.

    Returns
    -------
    y : np.ndarray (float32)
        The audio signal in [-1, 1].
    env : np.ndarray (float32, optional)
        The envelope used (if return_envelope=True).
    """
    assert carrier_hz > 0, "carrier_hz must be positive."
    assert pulse_rate_hz > 0, "pulse_rate_hz must be positive."
    assert duration_s > 0, "duration_s must be positive."
    assert sr > 0, "sr must be positive."
    assert 0 <= amplitude <= 1.0, "amplitude must be in [0, 1]."

    t = np.arange(int(duration_s * sr)) / sr
    period_s = 1.0 / pulse_rate_hz
    if on_duration_s > period_s:
        raise ValueError(
            f"on_duration_s ({on_duration_s}) must be <= pulse period ({period_s:.4f}s) for {pulse_rate_hz} Hz pulses."
        )

    # Base envelope: 1 during ON windows, 0 otherwise (rectangular)
    # ON when (t % period) < on_duration_s
    env = ((np.mod(t, period_s) < on_duration_s)).astype(np.float32)

    # Smooth edges with short attack/release to avoid clicks
    fade_s = max(0.0, fade_ms / 1000.0)
    if fade_s > 0:
        n_fade = int(round(fade_s * sr))
        if n_fade > 0:
            # Build a 1D fade kernel that rises from 0->1 then falls 1->0
            # We'll apply as convolution on the edges by constructing a smoothing window.
            # Simpler approach: convolve with a Hann window then renormalize.
            win_len = 2 * n_fade + 1
            # Hann window (smooth)
            win = np.hanning(win_len)
            win = win / win.max()  # normalize peak to 1
            # Convolve envelope with window (smoothing the transitions)
            smoothed = np.convolve(env, win, mode="same")
            # Rescale so that plateau remains approx 1
            # The plateau region of a box convolved with this window ~ sum(win) -> normalize by sum(win)
            smoothed = smoothed / np.sum(win)
            # Re-normalize to [0, 1] (safety due to edge effects)
            smoothed = np.clip(smoothed * (on_duration_s / (on_duration_s + 2 * fade_s)), 0.0, 1.0)
            env = smoothed.astype(np.float32)

    # Carrier
    phase = 2 * np.pi * carrier_hz * t + start_phase
    carrier = np.sin(phase).astype(np.float32)

    # Apply envelope and amplitude
    y = (amplitude * env * carrier).astype(np.float32)
    t = np.arange(len(y)) / sr
    if return_envelope:
        return y, t, env
    return y, t

def ROI_cluster_traceSets(figSaveDir,figSubDir,batchID,figName,ROI_clustering,figParams,zoom,close_all_figs,PDF_export_active=False):

    class _NoPdf:
        """Null stand-in for PdfPages: accepts savefig calls and writes nothing."""
        def __enter__(self):        return self
        def __exit__(self, *exc):   return False
        def savefig(self, *args, **kwargs): pass

    def maybe_pdf(path, active=None):
        active = PDF_export_active if active is None else active
        return PdfPages(path) if active else _NoPdf()

    #######################################################################################################
    if not 'clean' in figParams.keys():
        figParams['clean'] = False
    if not 'constrained_layout' in figParams.keys():
        figParams['constrained_layout'] = False
    if not 'legendLoc' in figParams.keys():
        figParams['legendLoc'] = 'upper right'
    if not 'legendBboxToAnchor' in figParams.keys():
        figParams['legendBboxToAnchor'] = (1, 1)
    if not 'filt_trace_byROI' in figParams.keys():
        figParams['filt_trace_byROI'] = False
    if not 'traceSet_split_ttypes' in figParams.keys():
        figParams['traceSet_split_ttypes'] = False
    if not 'traceSet_includeTrialCounts' in figParams.keys():
        figParams['traceSet_includeTrialCounts'] = True
    if not 'traceSet_overlay' in figParams.keys():
        figParams['traceSet_overlay'] = False
    if not 'trace_heightAdjust' in figParams.keys():
        figParams['trace_heightAdjust'] = 1
    if not 'manColWidth' in figParams.keys():
        figParams['manColWidth'] = 0
    if not 'anatColor' in figParams.keys():
        figParams['anatColor'] = False
    if not 'widthAdjust' in figParams.keys():
        if figParams['traceSet_includeTrialCounts']:
            figParams['widthAdjust'] = 0.1
        else:
            figParams['widthAdjust'] = 0.2
    if not 'gridspec_kw' in figParams.keys():
        figParams['gridspec_kw'] = None
    ROI_type_key = figParams['export_ROI_type_keys'][0]
    ROI_score = figParams['export_ROI_scores'][0]
    if len(figSubDir) > 0:
        figDir=os.path.join(figSaveDir,figSubDir)
        try:
            os.mkdir(figDir)
        except:
            pass
        figDir=os.path.join(figSaveDir,figSubDir,'Sets')
        try:
            os.mkdir(figDir)
        except:
            pass
        figDir=os.path.join(figSaveDir,figSubDir,'Sets',batchID)
        try:
            os.mkdir(figDir)
        except:
            pass
        figDir=os.path.join(figSaveDir,figSubDir,'Sets',batchID,figParams['trialGrouping'])
        try:
            os.mkdir(figDir)
        except:
            pass
        if zoom:
            figDir=os.path.join(figSaveDir,figSubDir,'Sets',batchID,figParams['trialGrouping'],figParams['zoom_label'])
        else:
            figDir=os.path.join(figSaveDir,figSubDir,'Sets',batchID,figParams['trialGrouping'],'Full')
            
    else:
        figDir=figSaveDir
        try:
            os.mkdir(figDir)
        except:
            pass
    if zoom:
        figName = figName + "_"+figParams['zoom_label']+".pdf"
    else:
        figName = figName + ".pdf"
    try:
        os.mkdir(figDir)
    except:
        pass
    print(f'figDir {figSaveDir}')
    print("Exporting: "+figName,end='')
    with maybe_pdf(os.path.join(figDir,figName)) as pdf:
        #######################################################################################################
        overlayColors = []
        overlayLabels = []
        overlayCmapLabel = ''
        if figParams['traceSet_overlay']:
            print("Collecting Overlay Info: ")
            for a in range(len(figParams['traceSet_trace_plotScalar'])):
                if figParams['traceSet_trace_plotScalar'][a] == 0:
                    print("<<WARNING>> ADJUSTING figParams['traceSet_trace_plotScalar'] from 0 -> 1")
                    figParams['traceSet_trace_plotScalar'][a] = 1
            
            for cg in range(len(figParams['export_traceSets'].keys())):
                for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                    if not figParams['export_traceSets'][cg][traceSet]['overlayLabel'] in overlayLabels:
                        print("Adding cg"+str(cg)+" traseSet"+str(traceSet))
                        overlayLabels.append(figParams['export_traceSets'][cg][traceSet]['overlayLabel'])
                        overlayColors.append(figParams['export_traceSets'][cg][traceSet]['overlayColor'])
                        overlayCmapLabel = overlayCmapLabel+","+str(figParams['export_traceSets'][cg][traceSet]['overlayLabel'])+":"+str(figParams['export_traceSets'][cg][traceSet]['overlayColor'])
            # print("overlayLabels = "+str(overlayLabels))
            # print("overlayColors = "+str(overlayColors))
            # print("overlayCmapLabel = "+str(overlayCmapLabel))
            
        foundTtypeScaling = False
        for cg in range(len(figParams['export_traceSets'].keys())):
            for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                    if 'allSess' in figParams['trialGrouping']:
                        if figParams['export_traceSets'][cg][traceSet]['ttypes'][s] in figParams['scaling_ttypes']:
                            foundTtypeScaling = True
                    elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                        if figParams['export_traceSets'][cg][traceSet]['ttypes'][s] in figParams['scaling_ttypes'] and figParams['export_traceSets'][cg][traceSet]['be'][s] in figParams['scaling_be']:
                            foundTtypeScaling = True
        if not foundTtypeScaling:
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
            print("<<WARNING>> No ttypes found in traceSets match the scaling_ttypes specified in figParams!")
        #######################################################################################################
        #Generate blank manual im contrast if needed
        if figParams['manual_im_contrast']:
            print("NOTE - Manual Image Contrasts: ")
            for group in figParams['manual_im_contrasts'].keys():
                for trialGrouping in figParams['manual_im_contrasts'][group].keys():
                    for align_data in figParams['manual_im_contrasts'][group][trialGrouping].keys():
                        for ROI_type_key in figParams['export_ROI_type_keys']:
                            if not ROI_type_key in figParams['manual_im_contrasts'][group][trialGrouping][align_data].keys():
                                figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key] = {}
                            for ROI_score in figParams['export_ROI_scores']:
                                if not ROI_score in figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key].keys():
                                    figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key][ROI_score] = {}

        if figParams['manual_plot_contrast']:
            print("NOTE - Manual Plot Contrasts: ")
            for group in figParams['manual_plot_contrasts'].keys():
                for trialGrouping in figParams['manual_plot_contrasts'][group].keys():
                    for align_data in figParams['manual_plot_contrasts'][group][trialGrouping].keys():
                        for ROI_type_key in figParams['export_ROI_type_keys']:
                            if not ROI_type_key in figParams['manual_plot_contrasts'][group][trialGrouping][align_data].keys():
                                figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key] = {}
                            for ROI_score in figParams['export_ROI_scores']:
                                if not ROI_score in figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key].keys():
                                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key][ROI_score] = {}
        ################################################################################################################################################################################################
        nTS = 0
        for cg in figParams['export_traceSets'].keys():
            nTS = np.max([nTS,len(figParams['export_traceSets'][cg].keys())])
        nCG=len(figParams['export_traceSets'].keys())
        nTtypes = 0
        for cg in range(len(figParams['export_traceSets'].keys())):
            for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                nTtypes = np.max([nTtypes,len(figParams['export_traceSets'][cg][traceSet]['ttypes'])])
        if figParams['traceSet_split_ttypes']:
            subRows=nTtypes+1
        else:
            subRows=2
        if figParams['traceSet_overlay']:
            nRealRows = 1
            nRows = subRows
            nTSreal = 1
            nCGreal = nCG
        else:
            nCGreal = 1
            nRealRows = nCG
            nRows = subRows * nCG
            nTSreal = copy.deepcopy(nTS)
        
        nrealCols = nTSreal * nCGreal* len((figParams['traceSet_export_align_data']))*len(figParams['export_ROI_type_keys'])*len(figParams['export_ROI_scores'])
        nCols = nTSreal * nCGreal * len((figParams['traceSet_export_align_data']))*len(figParams['export_ROI_type_keys'])*len(figParams['export_ROI_scores']) * 2
        
        if nRows == 1:
            nRows+=1
            delRow= True
        else:
            delRow = False
        if nCols == 1:
            nCols+=1
            delCol = True
        else:
            delCol = False
        # print("nCG = "+str(nCG))
        # print("nTS = "+str(nTS))
        # print("nTtypes = "+str(nTtypes))
        # print("nrealCols = "+str(nrealCols))
        # print("nCols = "+str(nCols))
        # print("nRealRows = "+str(nRealRows))
        # print("nRows = "+str(nRows))
        overlayMerges = {}
        horzLines = {}
        for row in range(nRows):
            overlayMerges[row] = {}
            for col in range(nCols):
                if figParams['traceSet_split_ttypes']:
                    overlayMerges[row][col] = {}
                    for s in range(nTtypes):
                        overlayMerges[row][col][s] = np.array([])
                else:
                    overlayMerges[row][col] = np.array([])
                horzLines[row,col] = []
        if not 'ROI_type_im_contrast' in figParams.keys():
            figParams['ROI_type_im_contrast'] = False
        if figParams['manual_im_contrast']:
            print("IMPORTANT Transferring manual_im_contrasts from figParams to im_contrasts...")
            im_contrasts = copy.deepcopy(figParams['manual_im_contrasts'][figParams['traceSet_group']][figParams['trialGrouping']])
        else:
            im_contrasts = {}
            for a,align_data in enumerate(figParams['traceSet_export_align_data']):
                im_contrasts[align_data] = {}
                im_contrasts[align_data]['im_maxVal'] = 0
                im_contrasts[align_data]['im_minVal'] = 0
                im_contrasts[align_data]['im_maxVal_cont'] = 0.1
                im_contrasts[align_data]['im_minVal_cont'] = 0
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    im_contrasts[align_data][ROI_type_key] = {}
                    for ROI_score in figParams['export_ROI_scores']:
                        im_contrasts[align_data][ROI_type_key][ROI_score] = {}
                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal'] = copy.deepcopy(im_contrasts[align_data]['im_maxVal'])
                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'] = copy.deepcopy(im_contrasts[align_data]['im_minVal'])
                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = copy.deepcopy(im_contrasts[align_data]['im_maxVal_cont'])
                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = copy.deepcopy(im_contrasts[align_data]['im_minVal_cont'])
        
        if figParams['clean']:
            if figParams['traceSet_group'] == 'trialStart':
                groupLabel = "Trial Start "
            elif figParams['traceSet_group'] == 'cTimes':
                groupLabel = "First Lick "
            elif figParams['traceSet_group'] == 'goCue':
                groupLabel = "Go Cue "
            elif figParams['traceSet_group'] == 'postContLick1':
                groupLabel = "2nd Lick  "
        else:
            if figParams['traceSet_group'] == 'trialStart':
                groupLabel = "Trial Start Aligned Time "
            elif figParams['traceSet_group'] == 'cTimes':
                groupLabel = "First Lick Aligned Time "
            elif figParams['traceSet_group'] == 'goCue':
                groupLabel = "Go Cue Aligned Time "
            elif figParams['traceSet_group'] == 'postContLick1':
                groupLabel = "Second Lick Aligned Time "


        if not 'ROI_type_plot_contrast' in figParams.keys():
            figParams['ROI_type_plot_contrast'] = False
        if figParams['manual_plot_contrast']:
            print("IMPORTANT Transferring manual_plot_contrasts from figParams to plot_contrasts...")

            plot_contrasts = copy.deepcopy(figParams['manual_plot_contrasts'][figParams['traceSet_group']][figParams['trialGrouping']])

        else:
            plot_contrasts = {}
            for a,align_data in enumerate(figParams['traceSet_export_align_data']):
                plot_contrasts[align_data] = {}
                plot_contrasts[align_data]['plot_shift_init'] = 0
                plot_contrasts[align_data]['plot_maxVal']=-1e6
                plot_contrasts[align_data]['plot_minVal']=1e6
                plot_contrasts[align_data]['plot_shift_final'] = 0
                plot_contrasts[align_data]['plot_maxVal_final']=-1e6
                plot_contrasts[align_data]['plot_minVal_final']=1e6
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    plot_contrasts[align_data][ROI_type_key] = {}
                    for ROI_score in figParams['export_ROI_scores']:
                        plot_contrasts[align_data][ROI_type_key][ROI_score] = {}
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_init'] = copy.deepcopy(plot_contrasts[align_data]['plot_shift_init'])
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal'] = copy.deepcopy(plot_contrasts[align_data]['plot_maxVal'])
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal'] = copy.deepcopy(plot_contrasts[align_data]['plot_minVal'])
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_final'] = copy.deepcopy(plot_contrasts[align_data]['plot_shift_final'])
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'] = copy.deepcopy(plot_contrasts[align_data]['plot_maxVal_final'])
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] = copy.deepcopy(plot_contrasts[align_data]['plot_minVal_final'])
        

        ####################################################################################
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            for cg in range(len(figParams['export_traceSets'].keys())):
                for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                    for ROI_type_key in figParams['export_ROI_type_keys']:
                        for ROI_score in figParams['export_ROI_scores']:
                            if not figParams['manual_im_contrast']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    invertData = False
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                        if ttype in figParams['scaling_ttypes']:
                                            if 'allSess' in figParams['trialGrouping']:
                                                for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                                    roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,False)
                                                    if invertData:
                                                        roiTrace = -roiTrace

                                                    im_contrasts[align_data]['im_maxVal']=np.nanmax([im_contrasts[align_data]['im_maxVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                    im_contrasts[align_data]['im_minVal']=np.nanmin([im_contrasts[align_data]['im_minVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                                                    im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']=np.nanmax([im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                    im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']=np.nanmin([im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if be in figParams['scaling_be']:
                                                    for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                                        # print("r = "+str(r)+", be = "+str(be),'ttype = '+str(ttype),'c = '+str(c),'align_data = '+str(align_data),'ROI_type_key = '+str(ROI_type_key),'ROI_score = '+str(ROI_score),'traceSet_group = '+str(figParams['traceSet_group']),'trialGrouping = '+str(figParams['trialGrouping']))
                                                        roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,False)
                                                        if invertData:
                                                            roiTrace = -roiTrace
                                                        # print(roiTrace)
                                                        im_contrasts[align_data]['im_maxVal']=np.nanmax([im_contrasts[align_data]['im_maxVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                        im_contrasts[align_data]['im_minVal']=np.nanmin([im_contrasts[align_data]['im_minVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                                                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']=np.nanmax([im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                        im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']=np.nanmin([im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                            if not figParams['manual_plot_contrast']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    invertData = False
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],\
                                            figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],plot_contrasts[align_data]['plot_shift_init']*s,invertData)
                                    if ttype in figParams['scaling_ttypes']:
                                        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                            if 'allSess' in figParams['trialGrouping']:
                                                plot_contrasts[align_data]['plot_maxVal']=np.nanmax([plot_contrasts[align_data]['plot_maxVal'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data]['plot_minVal']=np.nanmin([plot_contrasts[align_data]['plot_minVal'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal']=np.nanmax([plot_contrasts[align_data]['plot_maxVal'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal']=np.nanmin([plot_contrasts[align_data]['plot_minVal'],np.nanmin(lowerTrace)])
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],\
                                            figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_init']*s,invertData)
                                    if ttype in figParams['scaling_ttypes']:
                                        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                            if 'allSess' in figParams['trialGrouping']:
                                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if be in figParams['scaling_be']:
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal'],np.nanmin(lowerTrace)])
                                    
            if not figParams['manual_im_contrast']:
                if 'delta' in figParams['trialGrouping']:
                    if im_contrasts[align_data]['im_minVal']>=im_contrasts[align_data]['im_maxVal']:
                        print("Warning: im_minVal >= im_maxVal for align_data "+align_data+". Setting to -1 and 1. im_minVal = "+str(im_contrasts[align_data]['im_minVal'])+", im_maxVal = "+str(im_contrasts[align_data]['im_maxVal']))
                        im_contrasts[align_data]['im_minVal']=-1
                        im_contrasts[align_data]['im_maxVal']=1
                    im_contrasts[align_data]['im_minVal_cont'] = -1*np.nanmax([np.absolute(im_contrasts[align_data]['im_maxVal']),np.absolute(im_contrasts[align_data]['im_minVal'])])*figParams['traceSet_minCont'][a]
                    im_contrasts[align_data]['im_maxVal_cont'] = np.nanmax([np.absolute(im_contrasts[align_data]['im_maxVal']),np.absolute(im_contrasts[align_data]['im_minVal'])])*figParams['traceSet_maxCont'][a]
                else:
                    if im_contrasts[align_data]['im_minVal']>=im_contrasts[align_data]['im_maxVal']:
                        print("Warning: im_minVal >= im_maxVal for align_data "+align_data+". Setting to 0 and 1. im_minVal = "+str(im_contrasts[align_data]['im_minVal'])+", im_maxVal = "+str(im_contrasts[align_data]['im_maxVal']))
                        im_contrasts[align_data]['im_minVal']=0
                        im_contrasts[align_data]['im_maxVal']=1
                    im_contrasts[align_data]['im_minVal_cont'] = im_contrasts[align_data]['im_minVal']*figParams['traceSet_minCont'][a]
                    im_contrasts[align_data]['im_maxVal_cont'] = im_contrasts[align_data]['im_maxVal']*figParams['traceSet_maxCont'][a]
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    for ROI_score in figParams['export_ROI_scores']:
                        if 'delta' in figParams['trialGrouping']:
                            if im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']>=im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']:
                                print("Warning: im_minVal >= im_maxVal for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to -1 and 1. im_minVal = "+str(im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'])+", im_maxVal = "+str(im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']))
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']=-1
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']=1
                            
                            # figParams['traceSet_byROIType_minCont'] = {}
                            # figParams['traceSet_byROIType_maxCont'] = {}
                            # figParams['traceSet_byROIType_minCont'][a] = {}
                            # figParams['traceSet_byROIType_maxCont'][a] = {}
                            
                            # figParams['traceSet_byROIType_minCont'][a]['dendrites'] = 0.3
                            # figParams['traceSet_byROIType_maxCont'][a]['dendrites'] = 0.3
                            # figParams['traceSet_byROIType_minCont'][a]['somas'] = 0.6
                            # figParams['traceSet_byROIType_maxCont'][a]['somas'] = 0.6

                            if 'traceSet_byROIType_minCont' in figParams.keys():
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = -1*np.nanmax([np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']),\
                                    np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'])])*figParams['traceSet_byROIType_minCont'][a][ROI_type_key]
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = np.nanmax([np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']),\
                                    np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'])])*figParams['traceSet_byROIType_maxCont'][a][ROI_type_key]
                            else:
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = -1*np.nanmax([np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']),\
                                    np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'])])*figParams['traceSet_minCont'][a]
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = np.nanmax([np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']),\
                                    np.absolute(im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal'])])*figParams['traceSet_maxCont'][a]
                        else:
                            if im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']>=im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']:
                                print("Warning: im_minVal >= im_maxVal for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to 0 and 1.")
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']=0
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']=1
                            im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']*figParams['traceSet_minCont'][a]
                            im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']*figParams['traceSet_maxCont'][a]

                            if 'traceSet_byROIType_minCont' in figParams.keys():
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']*figParams['traceSet_byROIType_minCont'][a][ROI_type_key]
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']*figParams['traceSet_byROIType_maxCont'][a][ROI_type_key]

                            else:
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal']*figParams['traceSet_minCont'][a]
                                im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont'] = im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal']*figParams['traceSet_maxCont'][a]
                                
                        print(f'{align_data} {ROI_type_key} {ROI_score} im_minVal = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_minVal"]:3f}, im_maxVal = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_maxVal"]:3f}')
                        print(f'                                   im_minVal_cont = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_minVal_cont"]:3f}, im_maxVal_cont = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_maxVal_cont"]:3f}')

            if not figParams['manual_plot_contrast']:
                plot_contrasts[align_data]['plot_maxVal'] = plot_contrasts[align_data]['plot_maxVal']*1.05
                if 'delta' in figParams['trialGrouping']:
                    if plot_contrasts[align_data]['plot_minVal']>=plot_contrasts[align_data]['plot_maxVal']:
                        print("Warning: plot_minVal >= plot_maxVal for align_data "+align_data+". Setting to -1 and 1.")
                        plot_contrasts[align_data]['plot_minVal']=-1
                        plot_contrasts[align_data]['plot_maxVal']=1
                else:
                    if plot_contrasts[align_data]['plot_minVal']>=plot_contrasts[align_data]['plot_maxVal']:
                        print("Warning: plot_minVal >= plot_maxVal for align_data "+align_data+". Setting to 0 and 1.")
                        plot_contrasts[align_data]['plot_minVal']=0
                        plot_contrasts[align_data]['plot_maxVal']=1
                plot_contrasts[align_data]['plot_shift_final'] = figParams['traceSet_trace_plotScalar'][a]*(plot_contrasts[align_data]['plot_maxVal']-plot_contrasts[align_data]['plot_minVal'])
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    for ROI_score in figParams['export_ROI_scores']:
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal'] = plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']*1.05
                        if 'delta' in figParams['trialGrouping']:
                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']>=plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']:
                                print("Warning: plot_minVal >= plot_maxVal for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to -1 and 1.")
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']=-1
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']=1
                        else:
                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']>=plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']:
                                print("Warning: plot_minVal >= plot_maxVal for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to 0 and 1.")
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal']=0
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']=1
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_final'] = figParams['traceSet_trace_plotScalar'][a]\
                        *(plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal']-plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal'])

                for cg in range(len(figParams['export_traceSets'].keys())):
                    for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                        for ROI_type_key in figParams['export_ROI_type_keys']:
                            for ROI_score in figParams['export_ROI_scores']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    invertData = False
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],\
                                                                                             plot_contrasts[align_data]['plot_shift_final']*s,invertData)
                                    if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                        if 'delta' in figParams['trialGrouping']:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                        else:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                # print(ROI_type_key+" "+str(ROI_score)+" "+str(c)+" "+ttype+" "+be)
                                                # print(lowerTrace)
                                                plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],\
                                                                                             plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_final']*s,invertData)
                                    if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                        if 'delta' in figParams['trialGrouping']:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                        else:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] > 0:
                                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] = 0

                plot_contrasts[align_data]['plot_maxVal_final'] = \
                    plot_contrasts[align_data]['plot_maxVal_final']+\
                        (plot_contrasts[align_data]['plot_maxVal_final']-plot_contrasts[align_data]['plot_minVal_final'])*0.01
                plot_contrasts[align_data]['plot_minVal_final'] = \
                    plot_contrasts[align_data]['plot_minVal_final']-\
                        (plot_contrasts[align_data]['plot_maxVal_final']-plot_contrasts[align_data]['plot_minVal_final'])*0.02
                
                if 'delta' in figParams['trialGrouping']:
                    if plot_contrasts[align_data]['plot_minVal_final']>=plot_contrasts[align_data]['plot_maxVal_final']:
                        print("Warning: plot_minVal_final >= plot_maxVal_final for align_data "+align_data+". Setting to -1 and 1.")
                        plot_contrasts[align_data]['plot_minVal_final']=-1
                        plot_contrasts[align_data]['plot_maxVal_final']=1
                else:
                    if plot_contrasts[align_data]['plot_minVal_final']>=plot_contrasts[align_data]['plot_maxVal_final']:
                        print("Warning: plot_minVal_final >= plot_maxVal_final for align_data "+align_data+". Setting to 0 and 1.")
                        plot_contrasts[align_data]['plot_minVal_final']=0
                        plot_contrasts[align_data]['plot_maxVal_final']=1

                for ROI_type_key in figParams['export_ROI_type_keys']:
                    for ROI_score in figParams['export_ROI_scores']:
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final'] = \
                            plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']+\
                                (plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']-plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'])*0.01
                        plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] = \
                            plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']-\
                                (plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']-plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'])*0.02
                        if 'delta' in figParams['trialGrouping']:
                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']>=plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']:
                                print("Warning: plot_minVal_final >= plot_maxVal_final for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to -1 and 1.")
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=-1
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=1
                        else:
                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']>=plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']:
                                print("Warning: plot_minVal_final >= plot_maxVal_final for align_data "+align_data+", ROI type "+ROI_type_key+", ROI score "+str(ROI_score)+". Setting to 0 and 1.")
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']=0
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']=1
                            if plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] > 0:
                                plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final'] = 0

        max_nTrials = 0
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            for cg in range(len(figParams['export_traceSets'].keys())):
                for ts,traceSet in enumerate(figParams['export_traceSets'][cg].keys()):
                    for ROI_type_key in figParams['export_ROI_type_keys']:
                        for ROI_score in figParams['export_ROI_scores']:
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                    ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                else:
                                    ttypeLabel = ttype
                                c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                    if 'delta' in figParams['trialGrouping']:
                                        if 'allSess' in figParams['trialGrouping']:
                                            max_nTrials = np.nanmax([max_nTrials,np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'])])
                                        elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                            max_nTrials = np.nanmax([max_nTrials,np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'])])
                                    else:
                                        if 'allSess' in figParams['trialGrouping']:
                                            max_nTrials = np.nanmax([max_nTrials,np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials']),\
                                                                                np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['nTrials'])])
                                        elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                            max_nTrials = np.nanmax([max_nTrials,np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials']),\
                                                                                np.nanmax(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['nTrials'])])
        traceSets = {}
        for cg in range(len(figParams['export_traceSets'].keys())):
            traceSets[cg] = {}
            for ts in range(nTS):
                traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                traceSets[cg][traceSet] = {}
                for a,align_data in enumerate(figParams['traceSet_export_align_data']):
                    traceSets[cg][traceSet][align_data] = {}
                    for ROI_type_key in figParams['export_ROI_type_keys']:
                        traceSets[cg][traceSet][align_data][ROI_type_key] = {}
                        for ROI_score in figParams['export_ROI_scores']:
                            traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score] = {}
        allClusterTrialCounts = {}
        for cg in range(len(figParams['export_traceSets'].keys())):
            allClusterTrialCounts[cg] = {}
            for ts in range(nTS):
                traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                allClusterTrialCounts[cg][traceSet] = {}
                for a,align_data in enumerate(figParams['traceSet_export_align_data']):
                    allClusterTrialCounts[cg][traceSet][align_data] = {}
                    for ROI_type_key in figParams['export_ROI_type_keys']:
                        allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key] = {}
                        for ROI_score in figParams['export_ROI_scores']:
                            allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.array([])
        allClusterTrialCounts_byTtype = {}
        for cg in range(len(figParams['export_traceSets'].keys())):
            allClusterTrialCounts_byTtype[cg] = {}
            for ts in range(nTS):
                traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                allClusterTrialCounts_byTtype[cg][traceSet] = {}
                for a,align_data in enumerate(figParams['traceSet_export_align_data']):
                    allClusterTrialCounts_byTtype[cg][traceSet][align_data] = {}
                    for ROI_type_key in figParams['export_ROI_type_keys']:
                        allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key] = {}
                        for ROI_score in figParams['export_ROI_scores']:
                            allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score] = {}
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.array([])

        ################################################################################################################################################################################################
        figsize=(nCols*figParams['traceSet_hscalar'],nRows*figParams['traceSet_vscalar'])
        print(f'Generating figure with size {figsize} nRows {nRows} nCols {nCols}')
        fig,ax=clean_subplots(nRows,nCols,figsize=figsize,constrained_layout = figParams['constrained_layout'],gridspec_kw=figParams['gridspec_kw'])
        # plt.subplots_adjust(wspace=0.1, hspace=0.2)
        ################################################################################################################################################################################################
        #Images
        row = -1
        col = -2
        realCol = -1
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            ##################################
            align_data_short = copy.deepcopy(align_data)
            align_data_short = align_data_short.replace("allTC_","")
            ##################################
            for rt,ROI_type_key in enumerate(figParams['export_ROI_type_keys']):
                for ROI_score in figParams['export_ROI_scores']:
                    for cg in range(len(figParams['export_traceSets'].keys())):
                        realCol+=1
                        col+=2
                        realRow = -1
                        if figParams['traceSet_overlay']:
                            realRow+=1

                        for ts in range(nTS):
                            # print(f'align_data {align_data} realCol {realCol} col {col} realRow {realRow} nRows {nRows} ts {ts} ROI type {ROI_type_key} ROI score {ROI_score} trace group {cg}')
                            if not figParams['traceSet_overlay']:
                                realRow+=1
                            traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                            ####################################################################
                            #Prep Images
                            cmap_labels = ''
                            cmaps = {}
                            if 'cmap' in figParams['export_traceSets'][cg][traceSet]:
                                if len(figParams['export_traceSets'][cg][traceSet]['clusters']) == len(figParams['export_traceSets'][cg][traceSet]['cmap']):
                                    # print('using custom cmaps by ttype')
                                    tempCmap = copy.deepcopy(figParams['traceSet_cmap'])
                                    cmap,_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                    cmap.set_bad(color=figParams['nan_color'])
                                    cmap_label=str(tempCmap)

                                    for s in range(len(figParams['export_traceSets'][cg][traceSet]['clusters'])):
                                        tempCmap = copy.deepcopy(figParams['export_traceSets'][cg][traceSet]['cmap'][s])
                                        cmaps[s],_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                        cmaps[s].set_bad(color=figParams['nan_color'])
                                        cmap_labels=cmap_labels+str(tempCmap)+','
                                elif isinstance(figParams['export_traceSets'][cg][traceSet]['cmap'],tuple):
                                    tempCmap = copy.deepcopy(figParams['export_traceSets'][cg][traceSet]['cmap'])
                                    cmap,_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                    cmap.set_bad(color=figParams['nan_color'])
                                    cmap_label=str(tempCmap)
                                    tempCmap = copy.deepcopy(figParams['export_traceSets'][cg][traceSet]['cmap'])
                                    cmap_labels=cmap_labels+str(tempCmap)
                                    for s in range(len(figParams['export_traceSets'][cg][traceSet]['clusters'])):
                                        cmaps[s],_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                        cmaps[s].set_bad(color=figParams['nan_color'])
                                else:
                                    raise Exception("Error: Length of 'cmap' list does not match number of clusters in trace set "+traceSet+".")
                            else:
                                tempCmap = copy.deepcopy(figParams['traceSet_cmap'])
                                cmap,_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                cmap.set_bad(color=figParams['nan_color'])
                                cmap_label=str(tempCmap)

                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['clusters'])):
                                    tempCmap = copy.deepcopy(figParams['traceSet_cmap'])
                                    cmaps[s],_,_=generate_cmap(tempCmap,figParams['colorScalar'],cmap_name='new cmap')
                                    cmaps[s].set_bad(color=figParams['nan_color'])
                                    cmap_labels=cmap_labels+str(tempCmap)+','
                            if figParams['traceSet_overlay']:
                                cmap_labels = overlayCmapLabel
                            else:
                                if not figParams['traceSet_split_ttypes']:
                                    cmap_labels = cmap_label
                            cmap_labels = cmap_labels.rstrip(',')

                            ################################################
                            if not figParams['global_im_contrast'] and not figParams['ROI_type_im_contrast']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    invertData = False
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    if ttype in figParams['scaling_ttypes']:
                                        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                            if 'allSess' in figParams['trialGrouping']:
                                                for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                                    roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,False)
                                                    if invertData:
                                                        roiTrace = -roiTrace
                                                    im_contrasts[align_data]['im_maxVal']=np.nanmax([im_contrasts[align_data]['im_maxVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                    im_contrasts[align_data]['im_minVal']=np.nanmin([im_contrasts[align_data]['im_minVal'],\
                                                        np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if be in figParams['scaling_be']:
                                                    for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                                        roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,False)
                                                        if invertData:
                                                            roiTrace = -roiTrace
                                                        im_contrasts[align_data]['im_maxVal']=np.nanmax([im_contrasts[align_data]['im_maxVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_maxPer'])])
                                                        im_contrasts[align_data]['im_minVal']=np.nanmin([im_contrasts[align_data]['im_minVal'],\
                                                            np.nanpercentile(roiTrace,figParams['traceSet_image_minPer'])])
                                if 'delta' in figParams['trialGrouping']:
                                    im_contrasts[align_data]['im_minVal_cont'] = -1*np.nanmax([np.absolute(im_contrasts[align_data]['im_maxVal']),np.absolute(im_contrasts[align_data]['im_minVal'])])*figParams['traceSet_minCont'][a]
                                    im_contrasts[align_data]['im_maxVal_cont'] = np.nanmax([np.absolute(im_contrasts[align_data]['im_maxVal']),np.absolute(im_contrasts[align_data]['im_minVal'])])*figParams['traceSet_maxCont'][a]
                                else:
                                    im_contrasts[align_data]['im_minVal_cont'] = im_contrasts[align_data]['im_minVal']*figParams['traceSet_minCont'][a]
                                    im_contrasts[align_data]['im_maxVal_cont'] = im_contrasts[align_data]['im_maxVal']*figParams['traceSet_maxCont'][a]
                                print(f'{align_data} {ROI_type_key} {ROI_score} im_minVal_cont = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_minVal_cont"]:3f}, im_maxVal_cont = {im_contrasts[align_data][ROI_type_key][ROI_score]["im_maxVal_cont"]:3f}')

                            ################################################
                            xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                                load_alignment_xdata(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams'],ROI_type_key,figParams['traceSet_group'],zoom,figParams)
                            ################################################
                            # Collect traces
                            traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['nROIs'] = np.zeros(len(figParams['export_traceSets'][cg][traceSet]['ttypes']))
                            traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['orig_nTrials'] = {} 
                            traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['nTrials'] = {} 
                            emptyRow = np.ones((1,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['traceSet_group']][align_data]),dtype='float32')*np.nan
                            emptyCounts = np.ones((1,2),dtype='float32')*np.nan
                            edges = [0]
                            yticks = []
                            yticklabels = []
                            ytickcolors = []
                            if figParams['traceSet_addSpacers']:
                                allClusterData = copy.deepcopy(emptyRow)
                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = copy.deepcopy(emptyCounts)
                            else:
                                allClusterData = np.zeros((0,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['traceSet_group']][align_data]),dtype='float32')
                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.ones((0,2),dtype='float32')*np.nan
                            allClusterData_byTtype = {}
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.ones((0,2),dtype='float32')*np.nan
                                allClusterData_byTtype[s] = {}
                                allClusterData_byTtype[s]['allClusterData'] = np.zeros((0,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['traceSet_group']][align_data]),dtype='float32')
                                allClusterData_byTtype[s]['yticks'] = []
                                allClusterData_byTtype[s]['yticklabels'] = []
                                allClusterData_byTtype[s]['ytickcolors'] = []
                            # print(figParams['export_traceSets'][cg][traceSet]['ttypes'])
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                    ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                else:
                                    ttypeLabel = ttype
                                c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                invertData = False
                                if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                    if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                        invertData=True

                                # print(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])
                                if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                    yticks.append(edges[-1]+ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']/2)
                                    allClusterData_byTtype[s]['yticks'].append(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']/2)
                                    if figParams['traceSet_overlay']:
                                        clustLabel = 'OVERLAY'
                                    else:
                                        if 'cluster_simple_labels' in ROI_clustering[ROI_type_key][ROI_score].keys():
                                            clustLabel = ROI_clustering[ROI_type_key][ROI_score]['cluster_simple_labels'][c]
                                        else:
                                            clustLabel = 'Cluster '+str(c)
                                    if 'allSess' in figParams['trialGrouping']:
                                        nMasks = ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['cluster_summaryStats'][figParams['clustered_data']][figParams['plot_stat']]['nRepeats']
                                        traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['orig_nTrials'][s] = \
                                            copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'])
                                        traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['nTrials'][s] = \
                                            copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['nTrials'])
                                        if not figParams['clean']:
                                            label=clustLabel+"\n"+\
                                                    " All "+ttypeLabel+" Tr\n"+\
                                                    "(n="+str(nMasks)+"/"+str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])+","+\
                                                    str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['perROIs'])+"%)"
                                        else:
                                            label=ttype
                                    elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                        nMasks = ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['cluster_summaryStats'][figParams['clustered_data']][figParams['plot_stat']]['nRepeats']
                                        traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['orig_nTrials'][s] = \
                                            copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'])
                                        traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['nTrials'][s] = \
                                            copy.deepcopy(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['nTrials'])
                                        if not figParams['clean']:
                                            label= clustLabel+"\n"\
                                                    +be+" "+ttypeLabel+" Tr\n"+\
                                                    "(n="+str(nMasks)+"/"+str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])+","+\
                                                    str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['perROIs'])+"%)"
                                        else:
                                            label= be+"\n"+ttypeLabel

                                    yticklabels.append(label)
                                    allClusterData_byTtype[s]['yticklabels'].append(label)
                                    for r in range(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs']):
                                        traceSets[cg][traceSet][align_data][ROI_type_key][ROI_score]['nROIs'][s]+=1
                                        if 'delta' in figParams['trialGrouping']:
                                            if 'allSess' in figParams['trialGrouping']:
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData = np.concatenate((allClusterData,\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.concatenate((allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score],\
                                                                np.expand_dims(np.array(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'][r]),axis=0)),axis=0)
                                                
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData_byTtype[s]['allClusterData'] = np.concatenate((allClusterData_byTtype[s]['allClusterData'],\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.concatenate((allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s],\
                                                                np.expand_dims(np.array(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'][r]),axis=0)),axis=0)
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData = np.concatenate((allClusterData,\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.concatenate((allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score],\
                                                                np.expand_dims(np.array(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'][r]),axis=0)),axis=0)
                                                
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData_byTtype[s]['allClusterData'] = np.concatenate((allClusterData_byTtype[s]['allClusterData'],\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.concatenate((allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s],\
                                                                np.expand_dims(np.array(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'][r]),axis=0)),axis=0)
                                        else:
                                            if 'allSess' in figParams['trialGrouping']:
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData = np.concatenate((allClusterData,\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.concatenate((allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score],\
                                                                np.expand_dims(np.array([ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'][r],\
                                                                                        ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['nTrials'][r]]),axis=0)),axis=0)
                                                
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData_byTtype[s]['allClusterData'] = np.concatenate((allClusterData_byTtype[s]['allClusterData'],\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.concatenate((allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s],\
                                                                np.expand_dims(np.array([ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['orig_nTrials'][r],\
                                                                                        ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['nTrials'][r]]),axis=0)),axis=0)
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData = np.concatenate((allClusterData,\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.concatenate((allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score],\
                                                                np.expand_dims(np.array([ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'][r],\
                                                                                        ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['nTrials'][r]]),axis=0)),axis=0)
                                                
                                                roiTrace = extract_ROI_trace(ROI_clustering,ROI_type_key,ROI_score,c,figParams,figParams['traceSet_group'],figParams['trialGrouping'],align_data,be,ttype,r,figParams['traceSet_fillNaNs'])
                                                if invertData:
                                                    roiTrace = -roiTrace
                                                allClusterData_byTtype[s]['allClusterData'] = np.concatenate((allClusterData_byTtype[s]['allClusterData'],\
                                                                np.expand_dims(roiTrace,axis=0)),axis=0)
                                                allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s] = np.concatenate((allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s],\
                                                                np.expand_dims(np.array([ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['orig_nTrials'][r],\
                                                                                        ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['nTrials'][r]]),axis=0)),axis=0)
                                    
                                    ###########################
                                    if figParams['anatColor']:
                                        tempColor = figParams['anatColors'][rt]
                                    else:
                                        if 'colors' in figParams['export_traceSets'][cg][traceSet].keys():
                                            tempColor = figParams['export_traceSets'][cg][traceSet]['colors'][s]
                                        else:
                                            tempColor = ROI_clustering[ROI_type_key][ROI_score]['cluster_colors'][c]
                                        
                                    if figParams['traceSet_addSpacers']:
                                        allClusterData = np.concatenate((allClusterData,emptyRow),axis=0)
                                        allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score] = np.concatenate((allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score],emptyCounts),axis=0)
                                        edges.append(allClusterData.shape[0]-1)
                                    else:
                                        edges.append(allClusterData.shape[0])
            
                                    ytickcolors.append(tempColor)
                                    allClusterData_byTtype[s]['ytickcolors'].append(tempColor)
                            # print(ytickcolors)
                            #######################################################################################################
                            if not figParams['global_im_contrast'] and not figParams['ROI_type_im_contrast']:
                                maxVal = im_contrasts[align_data]['im_maxVal_cont']
                                minVal = im_contrasts[align_data]['im_minVal_cont']
                            elif figParams['global_im_contrast'] and not figParams['ROI_type_im_contrast']:
                                maxVal = im_contrasts[align_data]['im_maxVal_cont']
                                minVal = im_contrasts[align_data]['im_minVal_cont']
                            elif not figParams['global_im_contrast'] and figParams['ROI_type_im_contrast']:
                                maxVal = im_contrasts[align_data][ROI_type_key][ROI_score]['im_maxVal_cont']
                                minVal = im_contrasts[align_data][ROI_type_key][ROI_score]['im_minVal_cont']
                            else:
                                raise Exception("Cannot have both global_im_contrast and ROI_type_im_contrast set to True")
                            # print("figParams['global_im_contrast'] = "+str(figParams['global_im_contrast'])+", figParams['ROI_type_im_contrast'] = "+str(figParams['ROI_type_im_contrast']))
                            # print(minVal,maxVal)
                            # print(f'minVal {minVal} maxVal {maxVal} for align_data {align_data} ROI type {ROI_type_key} ROI score {ROI_score}')
                            if not figParams['clean']:
                                if figParams['filt_trace_byROI']:
                                    if align_data in figParams['filt_trace_byROI_window'].keys():
                                        if ROI_type_key in figParams['filt_trace_byROI_window'][align_data].keys():
                                            filtLabel = "Filt"+str(figParams['filt_trace_byROI_window'][align_data][ROI_type_key])+" "
                                        else:
                                            filtLabel = ""
                                    else:
                                        filtLabel = ""
                                else:
                                    filtLabel = ""
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    binLabel = ""
                                else:
                                    binLabel = "bin"+str(int(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key]))
                                if figParams['traceSet_overlay']:
                                    clustLabel = ' OVERLAY'+"\n"+\
                                        ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['overlayLabel1']+" | "
                                else:
                                    if 'cluster_detail_labels' in ROI_clustering[ROI_type_key][ROI_score].keys():
                                        clustLabel = "\n"+ROI_clustering[ROI_type_key][ROI_score]['cluster_detail_labels'][c]+"\n"+\
                                            ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['label']+" | "
                                    else:
                                        clustLabel = "\nCluster "+str(c)+"\n"+\
                                            ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['label']+" | "  
                                label1 = "("+filtLabel+"Img.Lims = "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+" "+str(cmap_labels)+")"
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    label = ROI_clustering[ROI_type_key][ROI_score]['detailed_clustering_label_title_oneLine']+\
                                        clustLabel+\
                                        figParams['traceSet_export_align_data_label'][a]+" | "+figParams['traceSet_group']+" "+figParams['plot_data_label_byROI']+\
                                        "\n("+filtLabel+"Img.Lims = "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+" "+str(cmap_labels)+")"
                                else:
                                    label = ROI_clustering[ROI_type_key][ROI_score]['detailed_clustering_label_title_oneLine']+\
                                        clustLabel+\
                                        figParams['traceSet_export_align_data_label'][a]+" "+\
                                            " ("+binLabel+") | "+\
                                            figParams['traceSet_group']+" "+figParams['plot_data_label_byROI']+\
                                        "\n("+filtLabel+"Img.Lims = "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+" "+str(cmap_labels)+")"
                            else:
                                if figParams['filt_trace_byROI']:
                                    if align_data in figParams['filt_trace_byROI_window'].keys():
                                        if ROI_type_key in figParams['filt_trace_byROI_window'][align_data].keys():
                                            filtLabel = "F"+str(figParams['filt_trace_byROI_window'][align_data][ROI_type_key])+" "
                                        else:
                                            filtLabel = ""
                                    else:
                                        filtLabel = ""
                                else:
                                    filtLabel = ""
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    binLabel = ""
                                else:
                                    binLabel = "B"+str(int(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key]))
                                if figParams['traceSet_overlay']:
                                    clustLabel = ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['overlayLabel1']+" | "
                                else:
                                    if 'cluster_detail_labels' in ROI_clustering[ROI_type_key][ROI_score].keys():
                                        clustLabel = ROI_clustering[ROI_type_key][ROI_score]['cluster_detail_labels'][c]+"\n"+\
                                            ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['label']+" | "
                                    else:
                                        clustLabel = ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+figParams['export_traceSets'][cg][traceSet]['label']+" | "
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    label = " ("+filtLabel+" "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+")"
                                    label1 = "("+filtLabel+" "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+")"
                                else:
                                    label = " ("+binLabel+"; "+filtLabel+""+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+")"
                                    label1 =  "("+binLabel+"; "+filtLabel+" "+str(np.round(minVal,decimals=1))+"<->"+str(np.round(maxVal,decimals=1))+")"

                            # print(f'label = {label}')
                            # print(f'label1 = {label1}')
                            # print(f'minVal = {minVal}, maxVal = {maxVal}')
                            allClusterData_RGB = convert2RGB(allClusterData,\
                                minVal,maxVal,\
                                cmap,figParams['colorScalar'],figParams['nan_color'])
                            allClusterData_RGB = allClusterData_RGB[:,:,0:3]
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                allClusterData_byTtype[s]['allClusterData_RGB'] = convert2RGB(allClusterData_byTtype[s]['allClusterData'],\
                                    minVal,maxVal,\
                                    cmaps[s],figParams['colorScalar'],figParams['nan_color'])
                                allClusterData_byTtype[s]['allClusterData_RGB'] = allClusterData_byTtype[s]['allClusterData_RGB'][:,:,0:3]
                            ################################################################################################################################################################################################
                            if figParams['traceSet_split_ttypes']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    if figParams['traceSet_overlay']:
                                        row = s
                                    else:
                                        row = (ts*subRows)+s
                                    if figParams['clean']:
                                        tempStr = str(yticklabels[s])
                                        tempStr = tempStr.replace('\n',' ')
                                        labela = tempStr+" "+label
                                        label1a = tempStr+" "+label1
                                    else:
                                        labela = label
                                        label1a = label1


                                    ################################################################################################
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    if overlayMerges[row][col][s].size == 0:
                                        overlayMerges[row][col][s] = copy.deepcopy(allClusterData_byTtype[s]['allClusterData_RGB'][:,plotIdxs,:])
                                    else:
                                        overlayMerges[row][col][s] = overlayMerges[row][col][s] + copy.deepcopy(allClusterData_byTtype[s]['allClusterData_RGB'][:,plotIdxs,:])
                                    overlayMerges[row][col][s][overlayMerges[row][col][s]>1] = 1
                                    # fig1,ax1 = plt.subplots(1,1)
                                    # ax1.imshow(overlayMerges[row][col][s])
                                    # plt.show()
                                    if figParams['traceSet_aspect']:
                                        ax[row,col].imshow(overlayMerges[row][col][s],extent=[xlim_s[0],xlim_s[1],overlayMerges[row][col][s],0],interpolation='none',aspect=figParams['traceSet_aspect'],clip_on=False)
                                    else:
                                        ax[row,col].imshow(overlayMerges[row][col][s],extent=[xlim_s[0],xlim_s[1],overlayMerges[row][col][s].shape[0],0],interpolation='none',aspect = 'auto',clip_on=False)
                                    ################################################################################################
                                    # if s == len(figParams['export_traceSets'][cg][traceSet]['ttypes'])-1 and ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                    if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                        # if not figParams['clean']:
                                        if figParams['traceSet_split_ttypes'] and s > 0:
                                            ax[row,col].text(xlim_s[0],0-figParams['ylim_adjust'][0],ROI_type_key[0:4]+" "+label1a,ha='left',va='bottom',color = allClusterData_byTtype[s]['ytickcolors'][c],fontsize=figParams['traceSet_fontSize'])
                                        else:
                                            ax[row,col].text(xlim_s[0],0-figParams['ylim_adjust'][0],ROI_type_key[0:4]+" "+labela,ha='left',va='bottom',color = allClusterData_byTtype[s]['ytickcolors'][c],fontsize=figParams['traceSet_fontSize'])
                                        ax[row,col] = imshow_cleanup(ax[row,col],False,False)
                                        # ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                        trial_structure_times = load_clustering_trial_structure_times(ROI_clustering,ROI_type_key,ROI_score,figParams['traceSet_group'],align_data,fix_overlaps=True,verbose = False)
                                        ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,\
                                                            figParams['traceSet_group'],align_data,ttype,figParams['export_traceSets'][cg][traceSet]['ttypes'],\
                                            True,True,False,False,figParams['traceSet_imshow_lineColor'],figParams['traceSet_imshow_lineWidth'],figParams['traceSet_alpha'],\
                                                2,1,0,0,0-figParams['ylim_adjust'][0],figParams['horzCueLineOn'],figParams['horzCueImScalar'])
                                        ax[row,col].set_yticks([])
                                        ax[row,col].set_ylim([allClusterData_byTtype[s]['allClusterData'].shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]]) 
                                        ax[row,col].set_yticks(allClusterData_byTtype[s]['yticks'])
                                        ax[row,col].tick_params(axis='y', which='both',length=0)
                                        if figParams['traceSet_split_ttypes'] and not figParams['clean']:
                                            ax[row,col].set_yticklabels(allClusterData_byTtype[s]['yticklabels'],fontdict = {'fontsize': figParams['traceSet_fontSize']-2,\
                                                'verticalalignment': figParams['traceSet_imshow_ytick_verticalalignment'],\
                                                    'horizontalalignment': figParams['traceSet_imshow_ytick_horizontalalignment'],\
                                                        'rotation':figParams['traceSet_imshow_ytick_rotation']})
                                            ytick_labels = ax[row,col].get_yticklabels()
                                            for c in range(len(allClusterData_byTtype[s]['ytickcolors'])):
                                                ytick_labels[c].set_color(allClusterData_byTtype[s]['ytickcolors'][c])   # First label
                                        else:
                                            ax[row,col].set_yticklabels([])
                                        ax[row,col].set_xlim([xlim_s[0],xlim_s[1]])
                                        if not figParams['clean']:
                                            ax[row,col].set_xticks(xticks_s)
                                            ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['traceSet_fontSize']-2,'verticalalignment': 'top','horizontalalignment': 'center'})
                                        else:
                                            ax[row,col].set_xticks([])
                                        figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                                        if len(figParams['im_scaleBar']['horzUnit'])>0:
                                            figParams['horz_scale_label']=' '+figParams['im_scaleBar']['horzUnit']
                                        else:
                                            figParams['horz_scale_label']=''
                                        if len(figParams['im_scaleBar']['vertUnit'])>0:
                                            figParams['vert_scale_label']=' '+figParams['im_scaleBar']['vertUnit']
                                        else:
                                            figParams['vert_scale_label']=''
                                        ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],\
                                            figParams['traceSet_imshow_lineColor'],True,figParams['vert_scale_label'],str(figParams['plot_scaleBar']['length_s'])+figParams['horz_scale_label'])
                                        pos1 = ax[row,col].get_position()
                                        ax[row,col].set_position([pos1.x0-(1/nrealCols)*figParams['widthAdjust']*realCol, pos1.y0, pos1.width+(1/nrealCols)*figParams['widthAdjust'], pos1.height])
                                        pos1 = ax[row,col].get_position()
                                ################################################################################################
                            else:
                                if figParams['traceSet_overlay']:
                                    row = 0   
                                else:
                                    row = (ts*subRows)+0   
                                if figParams['clean']:
                                    tempStr = str(yticklabels[s])
                                    tempStr = tempStr.replace('\n',' ')
                                    labela = tempStr+" "+label
                                    label1a = tempStr+" "+label1
                                else:
                                    labela = label
                                    label1a = label1
                                ################################################################################################
                                if overlayMerges[row][col].size == 0:
                                    overlayMerges[row][col] = copy.deepcopy(allClusterData_RGB[:,plotIdxs,:])
                                else:
                                    overlayMerges[row][col] = overlayMerges[row][col] + copy.deepcopy(allClusterData_RGB[:,plotIdxs,:])
                                overlayMerges[row][col][overlayMerges[row][col]>1] = 1
                                if figParams['traceSet_overlay']:
                                    if  ts == nTS-1:
                                        displayOn = True
                                    else:
                                        displayOn = False
                                else:
                                    displayOn = True
                                
                                if displayOn:
                                    if figParams['traceSet_aspect']:
                                        ax[row,col].imshow(overlayMerges[row][col],extent=[xlim_s[0],xlim_s[1],overlayMerges[row][col].shape[0],0],interpolation='none',aspect=figParams['traceSet_aspect'],clip_on=False)
                                    else:
                                        ax[row,col].imshow(overlayMerges[row][col],extent=[xlim_s[0],xlim_s[1],overlayMerges[row][col].shape[0],0],interpolation='none',aspect = 'auto',clip_on=False)
                                ################################################################################################
                                if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                    # if not figParams['clean']:
                                    if figParams['traceSet_split_ttypes'] and s > 0:
                                        ax[row,col].text(xlim_s[0],0-figParams['ylim_adjust'][0],label1a,ha='left',va='bottom',color = ytickcolors[c],fontsize=figParams['traceSet_fontSize'])
                                    else:
                                        ax[row,col].text(xlim_s[0],0-figParams['ylim_adjust'][0],labela,ha='left',va='bottom',color = ytickcolors[c],fontsize=figParams['traceSet_fontSize'])
                                    ax[row,col] = imshow_cleanup(ax[row,col],False,False)
                                    # ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                    trial_structure_times = load_clustering_trial_structure_times(ROI_clustering,ROI_type_key,ROI_score,figParams['traceSet_group'],align_data,fix_overlaps=True,verbose = False)
                                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['traceSet_group'],\
                                        align_data,ttype,figParams['export_traceSets'][cg][traceSet]['ttypes'],\
                                        True,True,False,False,figParams['traceSet_imshow_lineColor'],figParams['traceSet_imshow_lineWidth'],figParams['traceSet_alpha'],2,1,0,0,0-figParams['ylim_adjust'][0],figParams['horzCueLineOn'],figParams['horzCueImScalar'])
                                    if figParams['traceSet_clusterGrid']:
                                        ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                        ax[row,col].set_yticks([])
                                        for e in edges:
                                            if figParams['traceSet_addSpacers']:
                                                ax[row,col].axhline(e,color=figParams['traceSet_imshow_lineColor'],linewidth=figParams['traceSet_imshow_lineWidth'],alpha=figParams['traceSet_alpha'])
                                            else:
                                                ax[row,col].axhline(e-0.5,color=figParams['traceSet_imshow_lineColor'],linewidth=figParams['traceSet_imshow_lineWidth'],alpha=figParams['traceSet_alpha'])
                                    else:
                                        ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                        ax[row,col].set_yticks([])
                                        for c,yticklabel in enumerate(yticklabels):
                                            if figParams['traceSet_addSpacers']:
                                                ax[row,col].plot([xlim_s[0],xlim_s[0]],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['traceSet_imshow_lineWidth']+1,alpha=figParams['traceSet_alpha'])
                                            else:
                                                ax[row,col].plot([xlim_s[0],xlim_s[0]],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['traceSet_imshow_lineWidth']+1,alpha=figParams['traceSet_alpha'])
                                    if figParams['traceSet_addClusterLabels']:
                                        ax[row,col].set_yticks(yticks)
                                        ax[row,col].tick_params(axis='y', which='both',length=0)
                                        if figParams['traceSet_split_ttypes'] and not figParams['clean']:
                                            ax[row,col].set_yticklabels(yticklabels,fontdict = {'fontsize': figParams['traceSet_fontSize']-2,\
                                                'verticalalignment': figParams['traceSet_imshow_ytick_verticalalignment'],\
                                                    'horizontalalignment': figParams['traceSet_imshow_ytick_horizontalalignment'],\
                                                        'rotation':figParams['traceSet_imshow_ytick_rotation']})
                                            ytick_labels = ax[row,col].get_yticklabels()
                                            for c in range(len(ytickcolors)):
                                                ytick_labels[c].set_color(ytickcolors[c])  
                                        else:
                                            ax[row,col].set_yticklabels([])
                                    ax[row,col].set_xlim([xlim_s[0],xlim_s[1]])
                                    if not figParams['clean']:
                                        ax[row,col].set_xticks(xticks_s)
                                        ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['traceSet_fontSize']-2,'verticalalignment': 'top','horizontalalignment': 'center'})
                                    else:
                                        ax[row,col].set_xticks([])
                                    figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                                    if len(figParams['im_scaleBar']['horzUnit'])>0:
                                        figParams['horz_scale_label']=' '+figParams['im_scaleBar']['horzUnit']
                                    else:
                                        figParams['horz_scale_label']=''
                                    if len(figParams['im_scaleBar']['vertUnit'])>0:
                                        figParams['vert_scale_label']=' '+figParams['im_scaleBar']['vertUnit']
                                    else:
                                        figParams['vert_scale_label']=''
                                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],\
                                            figParams['traceSet_imshow_lineColor'],True,figParams['vert_scale_label'],\
                                                str(figParams['plot_scaleBar']['length_s'])+figParams['horz_scale_label'])
                                    pos1 = ax[row,col].get_position()
                                    ax[row,col].set_position([pos1.x0-(1/nrealCols)*figParams['widthAdjust']*realCol, pos1.y0, pos1.width+(1/nrealCols)*figParams['widthAdjust'], pos1.height])
                                    pos1 = ax[row,col].get_position()    

        ################################################################################################################################################################################################
        #Trial Counts
        row = -1
        col = -1
        realCol = -1
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            ##################################
            align_data_short = copy.deepcopy(align_data)
            align_data_short = align_data_short.replace("allTC_","")
            ##################################
            for rt,ROI_type_key in enumerate(figParams['export_ROI_type_keys']):
                for ROI_score in figParams['export_ROI_scores']:
                    for cg in range(len(figParams['export_traceSets'].keys())):
                        realCol+=1
                        col+=2
                        realRow = -1
                        if figParams['traceSet_overlay']:
                            realRow+=1
                        for ts in range(nTS):
                            if not figParams['traceSet_overlay']:
                                realRow+=1
                            traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                            if figParams['traceSet_split_ttypes']:
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    if figParams['traceSet_overlay']:
                                        row = s
                                    else:
                                        row = (ts*subRows)+s
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    if figParams['traceSet_includeTrialCounts']:
                                        ################################################################################################
                                        if 'delta' in figParams['trialGrouping']:
                                            if figParams['traceSet_overlay']:
                                                tempColors = [overlayColors[ts],tuple(np.array(overlayColors[ts])*0.5)]
                                                tempLabels = [overlayLabels[ts]+' Trials',overlayLabels[ts]+' -Trials']
                                            else:
                                                tempColors = [(0.3,0.3,0.3),(0.6,0.6,0.6)]
                                                tempLabels = ['Trials','-Trials']
                                            figParams['export_traceSets'][cg][traceSet]['overlayColor']
                                            ax[row,col].plot(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s][:,0],\
                                                             np.arange(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s].shape[0])-0.5,'.-',\
                                                linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[0],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[0])
                                            ax[row,col].plot(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s][:,1],\
                                                             np.arange(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s].shape[0])-0.5,'.-',\
                                                linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[1],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[1])
                                        else:
                                            if figParams['traceSet_overlay']:
                                                tempColors = [overlayColors[ts],tuple(np.array(overlayColors[ts])*0.5)]
                                                tempLabels = [overlayLabels[ts]+' Trials',overlayLabels[ts]+' Clean Trials']
                                            else:
                                                tempColors = [(0,0,0),(0.5,0.5,0.5)]
                                                tempLabels = ['Trials','Clean Trials']
                                            ax[row,col].plot(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s][:,0],\
                                                             np.arange(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s].shape[0])-0.5,'.-',\
                                                linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[0],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[0])
                                            ax[row,col].plot(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s][:,1],\
                                                             np.arange(allClusterTrialCounts_byTtype[cg][traceSet][align_data][ROI_type_key][ROI_score][s].shape[0])-0.5,'.:',\
                                                linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[1],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[1])
                                        ##################################
                                        # if s == len(figParams['export_traceSets'][cg][traceSet]['ttypes'])-1 and ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                        if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                            ax[row,col].legend(frameon=False,fontsize=figParams['traceSet_fontSize']-2,
                                                           handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 

                                            ax[row,col].spines['top'].set_visible(False)
                                            ax[row,col].spines['right'].set_visible(False)
                                            ax[row,col].spines['bottom'].set_visible(False)
                                            ax[row,col].spines['left'].set_visible(False)
                                            ax[row,col].set_yticks([])
                                            ax[row,col].set_ylim([allClusterData_byTtype[s]['allClusterData'].shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                            ax[row,col].set_xlim([-10,max_nTrials+10])
                                            ax[row,col].tick_params(axis='x', which='major', labelsize=figParams['traceSet_fontSize']-2)
                                            ax[row,col].set_xlabel('Number of Trials',fontsize=figParams['traceSet_fontSize'])
                                            pos1 = ax[row,col-1].get_position()
                                            pos2 = ax[row,col].get_position()
                                            ax[row,col].set_position([pos1.x0+pos1.width, pos2.y0, pos2.width*0.35, pos2.height])
                                            pos1 = ax[row,col].get_position()
                                        ################################################################################################
                                    else:
                                        if (figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']:
                                            fig.delaxes(ax[row,col])                                        
                            else:
                                if figParams['traceSet_overlay']:
                                    row = 0   
                                else:
                                    row = (ts*subRows)+0   
                                if figParams['traceSet_includeTrialCounts']:
                                    ################################################################################################
                                    if 'delta' in figParams['trialGrouping']:
                                        tempLabels = ['Trials','-Trials']
                                        if figParams['traceSet_overlay']:
                                            tempColors = [overlayColors[ts],tuple(np.array(overlayColors[ts])*0.5)]
                                            tempLabels = [overlayLabels[ts]+' Trials',overlayLabels[ts]+' -Trials']
                                        else:
                                            tempColors = [(0.3,0.3,0.3),(0.6,0.6,0.6)]
                                            tempLabels = ['Trials','-Trials']
                                        ax[row,col].plot(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score][:,0],\
                                                         np.arange(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score].shape[0])-0.5,'.-',\
                                            linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[0],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[0])
                                        ax[row,col].plot(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score][:,1],\
                                                         np.arange(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score].shape[0])-0.5,'.-',\
                                            linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[1],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[1])
                                    else:
                                        if figParams['traceSet_overlay']:
                                            tempColors = [overlayColors[ts],tuple(np.array(overlayColors[ts])*0.5)]
                                            tempLabels = [overlayLabels[ts]+' Trials',overlayLabels[ts]+' Clean Trials']
                                        else:
                                            tempColors = [(0,0,0),(0.5,0.5,0.5)]
                                            tempLabels = ['Trials','Clean Trials']
                                        ax[row,col].plot(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score][:,0],
                                                         np.arange(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score].shape[0])-0.5,'.-',\
                                            linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[0],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[0])
                                        ax[row,col].plot(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score][:,1],
                                                         np.arange(allClusterTrialCounts[cg][traceSet][align_data][ROI_type_key][ROI_score].shape[0])-0.5,'.:',\
                                            linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColors[1],markersize=figParams['traceSet_plot_markerSize'],label=tempLabels[1])
                                    ##################################
                                    if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                        ax[row,col].legend(frameon=False,fontsize=figParams['traceSet_fontSize']-2,\
                                                    handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 

                                        ax[row,col].spines['top'].set_visible(False)
                                        ax[row,col].spines['right'].set_visible(False)
                                        ax[row,col].spines['bottom'].set_visible(False)
                                        ax[row,col].spines['left'].set_visible(False)
                                        if figParams['traceSet_clusterGrid']:
                                            ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                            ax[row,col].set_yticks([])
                                            for e in edges:
                                                if figParams['traceSet_addSpacers']:
                                                    ax[row,col].axhline(e,color=(0,0,0),linewidth=figParams['traceSet_imshow_lineWidth'],alpha=figParams['traceSet_alpha'])
                                                else:
                                                    ax[row,col].axhline(e-0.5,color=(0,0,0),linewidth=figParams['traceSet_imshow_lineWidth'],alpha=figParams['traceSet_alpha'])
                                        else:
                                            ax[row,col].set_ylim([allClusterData.shape[0]+figParams['ylim_adjust'][1],0-figParams['ylim_adjust'][0]])
                                            ax[row,col].set_yticks([])
                                            for c,yticklabel in enumerate(yticklabels):
                                                if figParams['traceSet_addSpacers']:
                                                    ax[row,col].plot([-5,-5],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['traceSet_imshow_lineWidth']+1,alpha=figParams['traceSet_alpha'])
                                                else:
                                                    ax[row,col].plot([-5,-5],[edges[c],edges[c+1]],color=ytickcolors[c],linewidth=figParams['traceSet_imshow_lineWidth']+1,alpha=figParams['traceSet_alpha'])
                                            
                                        ax[row,col].set_xlim([-10,max_nTrials+10])
                                        ax[row,col].tick_params(axis='x', which='major', labelsize=figParams['traceSet_fontSize']-2)
                                        ax[row,col].set_xlabel('Number of Trials',fontsize=figParams['traceSet_fontSize'])
                                        pos1 = ax[row,col-1].get_position()
                                        pos2 = ax[row,col].get_position()
                                        ax[row,col].set_position([pos1.x0+pos1.width, pos2.y0, pos2.width*0.35, pos2.height])
                                        pos1 = ax[row,col].get_position()
                                    ################################################################################################
                                else:
                                    if (figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']:
                                        fig.delaxes(ax[row,col])
        ################################################################################################################################################################################################
        #PLOTS
        row = 0
        col = -2
        realCol = -1
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            ##################################
            align_data_short = copy.deepcopy(align_data)
            align_data_short = align_data_short.replace("allTC_","")
            ##################################
            for rt,ROI_type_key in enumerate(figParams['export_ROI_type_keys']):
                for ROI_score in figParams['export_ROI_scores']:
                    for cg in range(len(figParams['export_traceSets'].keys())):
                        realCol+=1
                        col+=2
                        realRow = -1
                        if figParams['traceSet_overlay']:
                            realRow+=1
                        for ts in range(nTS):
                            if not figParams['traceSet_overlay']:
                                realRow+=1
                            traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                            if figParams['traceSet_split_ttypes']:
                                if figParams['traceSet_overlay']:
                                    row = subRows-1
                                else:
                                    row = (ts*subRows)+len(figParams['export_traceSets'][cg][traceSet]['ttypes'])
                                    # print(row,ts,subRows)
                            else:
                                if figParams['traceSet_overlay']:
                                    row = 1
                                else:
                                    row = (ts*subRows)+1
                            if not close_all_figs:
                                print("Adding plots TraceSet: "+traceSet+" ts: "+str(ts)+" row: "+str(row)+" col: "+str(col)+" ttypes "+str(figParams['export_traceSets'][cg][traceSet]['ttypes']))
                            ########################################################################################################################################
                            max_nFr = 0
                            max_s = 0
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                    ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                else:
                                    ttypeLabel = ttype
                                c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                    max_nFr = np.nanmax([max_nFr,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['traceSet_group']][align_data]])
                                    max_s = np.nanmax([max_s,ROI_clustering[ROI_type_key][ROI_score]['all_trace_nFr'][figParams['traceSet_group']][align_data]/\
                                                       ROI_clustering[ROI_type_key][ROI_score]['all_trace_framerate'][figParams['traceSet_group']][align_data]])
                            if not figParams['global_plot_contrast'] and not figParams['ROI_type_plot_contrast']:
                                plot_contrasts[align_data]['plot_shift_init'] = 0
                                plot_contrasts[align_data]['plot_maxVal']=-1e6
                                plot_contrasts[align_data]['plot_minVal']=1e6
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],\
                                                                                             plot_contrasts[align_data]['plot_shift_init']*s,invertData)

                                    if ttype in figParams['scaling_ttypes']:
                                        if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                            if 'allSess' in figParams['trialGrouping']:
                                                plot_contrasts[align_data]['plot_maxVal']=np.nanmax([plot_contrasts[align_data]['plot_maxVal'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data]['plot_minVal']=np.nanmin([plot_contrasts[align_data]['plot_minVal'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal']=np.nanmax([plot_contrasts[align_data]['plot_maxVal'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal']=np.nanmin([plot_contrasts[align_data]['plot_minVal'],np.nanmin(lowerTrace)])
                                plot_contrasts[align_data]['plot_maxVal'] = plot_contrasts[align_data]['plot_maxVal']*1.05
                                if 'delta' in figParams['trialGrouping']:
                                    if plot_contrasts[align_data]['plot_minVal']>=plot_contrasts[align_data]['plot_maxVal']:
                                        print("Warning: minVal >= maxVal for plot_contrasts calculation, setting to -1 to 1 range")
                                        plot_contrasts[align_data]['plot_minVal']=-1
                                        plot_contrasts[align_data]['plot_maxVal']=1
                                else:
                                    if plot_contrasts[align_data]['plot_minVal']>=plot_contrasts[align_data]['plot_maxVal']:
                                        print("Warning: minVal >= maxVal for plot_contrasts calculation, setting to 0-1 range")
                                        plot_contrasts[align_data]['plot_minVal']=0
                                        plot_contrasts[align_data]['plot_maxVal']=1
                            if not figParams['global_plot_contrast']:
                                plot_contrasts[align_data]['plot_shift_final'] = figParams['traceSet_trace_plotScalar'][a]*(plot_contrasts[align_data]['plot_maxVal']-plot_contrasts[align_data]['plot_minVal'])
                                plot_contrasts[align_data]['plot_maxVal_final']=-1e6
                                plot_contrasts[align_data]['plot_minVal_final']=1e6
                                for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                    be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                    invertData = False
                                    if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                        if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                            invertData=True
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],\
                                            figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],plot_contrasts[align_data]['plot_shift_final']*s,invertData)
                                    if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                        if 'delta' in figParams['trialGrouping']:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                    plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                        else:
                                            if 'allSess' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                            elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                                if ttype in figParams['scaling_ttypes'] and be in figParams['scaling_be']:
                                                    plot_contrasts[align_data]['plot_maxVal_final']=np.nanmax([plot_contrasts[align_data]['plot_maxVal_final'],np.nanmax(upperTrace)])
                                                plot_contrasts[align_data]['plot_minVal_final']=np.nanmin([plot_contrasts[align_data]['plot_minVal_final'],np.nanmin(lowerTrace)])
                                plot_contrasts[align_data]['plot_maxVal_final'] = plot_contrasts[align_data]['plot_maxVal_final']+\
                                    (plot_contrasts[align_data]['plot_maxVal_final']-plot_contrasts[align_data]['plot_minVal_final'])*0.01
                                plot_contrasts[align_data]['plot_minVal_final'] = plot_contrasts[align_data]['plot_minVal_final']-\
                                    (plot_contrasts[align_data]['plot_maxVal_final']-plot_contrasts[align_data]['plot_minVal_final'])*0.02
                                if 'delta' in figParams['trialGrouping']:
                                    if plot_contrasts[align_data]['plot_minVal_final']>=plot_contrasts[align_data]['plot_maxVal_final']:
                                        print("Warning: minVal >= maxVal for plot_contrasts final calculation, setting to -1 to 1 range")
                                        print("plot_contrasts[align_data]['plot_minVal_final'] = "+str(plot_contrasts[align_data]['plot_minVal_final'])+\
                                              ", plot_contrasts[align_data]['plot_maxVal_final'] = "+str(plot_contrasts[align_data]['plot_maxVal_final']))
                                        plot_contrasts[align_data]['plot_minVal_final']=-1
                                        plot_contrasts[align_data]['plot_maxVal_final']=1
                                else:
                                    if plot_contrasts[align_data]['plot_minVal_final']>=plot_contrasts[align_data]['plot_maxVal_final']:
                                        print("Warning: minVal >= maxVal for plot_contrasts final calculation, setting to 0-1 range")
                                        print("plot_contrasts[align_data]['plot_minVal_final'] = "+str(plot_contrasts[align_data]['plot_minVal_final'])+\
                                              ", plot_contrasts[align_data]['plot_maxVal_final'] = "+str(plot_contrasts[align_data]['plot_maxVal_final']))
                                        plot_contrasts[align_data]['plot_minVal_final']=0
                                        plot_contrasts[align_data]['plot_maxVal_final']=1
                                    if plot_contrasts[align_data]['plot_minVal_final'] > 0:
                                        plot_contrasts[align_data]['plot_minVal_final'] = 0
                            ####################################################################
                            xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                                load_alignment_xdata(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams'],ROI_type_key,figParams['traceSet_group'],zoom,figParams)
                            ####################################################################
                            plotLabels = []
                            for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                    ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                else:
                                    ttypeLabel = ttype
                                c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                invertData = False
                                if 'invert' in figParams['export_traceSets'][cg][traceSet].keys():
                                    if figParams['export_traceSets'][cg][traceSet]['invert'][s]:
                                        invertData=True
                                if figParams['traceSet_overlay']:
                                    clustLabel = 'OVERLAY'
                                else:
                                    if 'cluster_simple_labels' in ROI_clustering[ROI_type_key][ROI_score].keys():
                                        clustLabel = ROI_clustering[ROI_type_key][ROI_score]['cluster_simple_labels'][c]
                                    else:
                                        clustLabel = 'Cluster '+str(c)
                                if 'allSess' in figParams['trialGrouping']:
                                    nMasks = ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][ttype]['cluster_summaryStats'][figParams['clustered_data']][figParams['plot_stat']]['nRepeats']
                                    label=ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+"\n"+\
                                        clustLabel+"\n"+\
                                        " All "+ttypeLabel+" Tr"+\
                                        "\n(n="+str(nMasks)+"/"+str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])+","+\
                                            str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['perROIs'])+"%)"
                                elif 'byBehavEpoch' in figParams['trialGrouping'] or 'byPrePost' in figParams['trialGrouping'] or 'deltaPrePost' in figParams['trialGrouping']:
                                    nMasks = ROI_clustering[ROI_type_key][ROI_score]['clusters'][c][figParams['traceSet_group']][align_data][figParams['trialGrouping']][be][ttype]['cluster_summaryStats'][figParams['clustered_data']][figParams['plot_stat']]['nRepeats']
                                    label=ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+"\n"+\
                                        clustLabel+"\n"+\
                                        be+" "+ttypeLabel+" Trs"+\
                                        "\n(n="+str(nMasks)+"/"+str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'])+","+\
                                            str(ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['perROIs'])+"%)"
                                plotLabels.append(label)
                                if not figParams['global_plot_contrast'] and not figParams['ROI_type_plot_contrast']:
                                    maxVal = plot_contrasts[align_data]['plot_maxVal_final']
                                    minVal = plot_contrasts[align_data]['plot_minVal_final']
                                    shiftVal = plot_contrasts[align_data]['plot_shift_final']*s
                                elif figParams['global_plot_contrast'] and not figParams['ROI_type_plot_contrast']:
                                    maxVal = plot_contrasts[align_data]['plot_maxVal_final']
                                    minVal = plot_contrasts[align_data]['plot_minVal_final']
                                    shiftVal = plot_contrasts[align_data]['plot_shift_final']*s
                                elif not figParams['global_plot_contrast'] and figParams['ROI_type_plot_contrast']:
                                    maxVal = plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_maxVal_final']
                                    minVal = plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_minVal_final']
                                    shiftVal = plot_contrasts[align_data][ROI_type_key][ROI_score]['plot_shift_final']*s
                                else:
                                    raise Exception("Error: cannot have both global_plot_contrast and ROI_type_plot_contrast set to True")
                                if not -1*shiftVal in horzLines[row,col]:
                                    ax[row,col].plot([xdata_s[0],xdata_s[-1]],[-1*shiftVal,-1*shiftVal],linewidth=0.5,alpha=0.5,color=(0,0,0))
                                    horzLines[row,col].append(-1*shiftVal)
                                if s == 0 and ((figParams['traceSet_overlay'] and ts == 0) or not figParams['traceSet_overlay']):
                                    # print("Plotting "+label+" with minVal="+str(minVal)+", maxVal="+str(maxVal)+", shiftVal="+str(shiftVal))
                                    ax[row,col].spines['top'].set_visible(False)
                                    ax[row,col].spines['right'].set_visible(False)
                                    # ax[row,col].spines['left'].set_visible(False)
                                    # if 'delta' in figParams['trialGrouping']:
                                    trial_structure_times = load_clustering_trial_structure_times(ROI_clustering,ROI_type_key,ROI_score,figParams['traceSet_group'],align_data,fix_overlaps=True,verbose = False)
                                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['traceSet_group'],\
                                        align_data,ttype,figParams['export_traceSets'][cg][traceSet]['ttypes'],\
                                        False,True,True,False,(0,0,0),1,0.5,2,1,0,minVal,maxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
                                if figParams['traceSet_overlay']:
                                    tempLabel = overlayLabels[ts]+' '+ttypeLabel+" n="+str(nMasks)
                                    tempColor = overlayColors[ts]
                                else:
                                    tempLabel = ttypeLabel+" n="+str(nMasks)
                                    if figParams['anatColor']:
                                        tempColor = figParams['anatColors'][rt]
                                    else:
                                        if 'colors' in figParams['export_traceSets'][cg][traceSet].keys():
                                            tempColor = figParams['export_traceSets'][cg][traceSet]['colors'][s]
                                        else:
                                            tempColor = ROI_clustering[ROI_type_key][ROI_score]['cluster_colors'][c]
                                if figParams['anatColor']:
                                    tempColor1 = figParams['anatColors'][rt]
                                else:
                                    if 'colors' in figParams['export_traceSets'][cg][traceSet].keys():
                                        tempColor1 = figParams['export_traceSets'][cg][traceSet]['colors'][s]
                                    else:
                                        tempColor1 = ROI_clustering[ROI_type_key][ROI_score]['cluster_colors'][c]
                                if ROI_clustering[ROI_type_key][ROI_score]['clusters'][c]['nROIs'] > 0:
                                    mainTrace,upperTrace,lowerTrace = extract_cluster_traces(ROI_clustering,ROI_type_key,ROI_score,c,figParams['traceSet_group'],align_data,figParams['trialGrouping'],\
                                        figParams['clustered_data'],ttype,be,figParams['plot_stat'],figParams['plot_data'],figParams['plot_data_errors'],figParams['plot_data_error_pos'],figParams['plot_data_error_neg'],shiftVal,invertData)
                                    # print("upperTrace: "+str(upperTrace[plotIdxs]))
                                    # print("mainTrace: "+str(mainTrace[plotIdxs]))
                                    # print("lowerTrace: "+str(lowerTrace[plotIdxs]))
                                    if np.any(np.isfinite(mainTrace)):
                                        ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs],lowerTrace[plotIdxs],color=tempColor,ls='-',alpha = figParams['alpha'],linewidth=0,edgecolor='none')
                                        ax[row,col].plot(xdata_s[plotIdxs],mainTrace[plotIdxs],linewidth=figParams['traceSet_plot_lineWidth'],alpha=figParams['traceSet_plot_line_alpha'],color=tempColor,label=tempLabel)
                                    if ((figParams['traceSet_overlay'] and ts == 0) or not figParams['traceSet_overlay']):
                                        if not figParams['traceSet_plot_full_legend'] and np.isfinite(np.nanmean(mainTrace[plotIdxs])) and not figParams['traceSet_trace_plotScalar'][a] == 0:
                                            if not figParams['clean']:
                                                ax[row,col].text(xlim_s[0],np.nanmean(mainTrace[plotIdxs]),plotLabels[s],\
                                                    color = tempColor1,fontsize=figParams['traceSet_fontSize']-2,\
                                                    va=figParams['traceSet_plot_ytick_verticalalignment'],\
                                                    ha=figParams['traceSet_plot_ytick_horizontalalignment'],\
                                                    rotation=figParams['traceSet_plot_ytick_rotation'])
                            ####################################################################
                            if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    binLabel = ""
                                else:
                                    binLabel = "bin"+str(int(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key]))
                                
                                if figParams['traceSet_overlay']:
                                    clustLabel = "OVERLAY "+figParams['export_traceSets'][cg][traceSet]['overlayLabel1']
                                else:
                                    clustLabel = figParams['export_traceSets'][cg][traceSet]['label']
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    label = ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+clustLabel+" | "+\
                                        figParams['traceSet_export_align_data_label'][a]+" | "+figParams['traceSet_group']+" "+figParams['plot_data_label_byROI']
                                else:
                                    label = ROI_clustering[ROI_type_key][ROI_score]['ROI_type']+" "+clustLabel+" | "+\
                                        figParams['traceSet_export_align_data_label'][a]+" ("+binLabel+") | "+\
                                            figParams['traceSet_group']+" "+figParams['plot_data_label_byROI']
                                # if not figParams['traceSet_split_ttypes']:
                                # print(f'minVal {minVal:.3f} maxVal {maxVal:.3f}')
                                if not figParams['clean']:
                                    ax[row,col].text(xlim_s[0],maxVal,label,ha='left',va='bottom',color = (0,0,0),fontsize=figParams['traceSet_fontSize'])
    
                                if figParams['traceSet_plot_full_legend']:
                                    for s in range(len(figParams['export_traceSets'][cg][traceSet]['ttypes'])):
                                        ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                        if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                            ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                        else:
                                            ttypeLabel = ttype
                                        c = figParams['export_traceSets'][cg][traceSet]['clusters'][s]
                                        be = figParams['export_traceSets'][cg][traceSet]['be'][s]
                                        if figParams['anatColor'] and not figParams['traceSet_overlay']:
                                            tempColor = figParams['anatColors'][rt]
                                        elif figParams['traceSet_overlay'] and not figParams['anatColor']:
                                            tempColor = overlayColors[ts]
                                        elif not figParams['anatColor'] and not figParams['traceSet_overlay']:
                                            if 'colors' in figParams['export_traceSets'][cg][traceSet].keys():
                                                tempColor = figParams['export_traceSets'][cg][traceSet]['colors'][s]
                                            else:
                                                tempColor = ROI_clustering[ROI_type_key][ROI_score]['cluster_colors'][c]
                                        else:
                                            raise Exception("Error: cannot have both anatColor and traceSet_overlay set to True")
                                        ax[row,col].plot(-1000,-1000,linewidth=1,alpha=0.8,color=tempColor,\
                                            label=plotLabels[s])
                                    ax[row,col].legend(loc=figParams['legendLoc'],frameon=False, bbox_to_anchor=figParams['legendBboxToAnchor'], ncol = 1,fontsize=figParams['traceSet_fontSize']-2,\
                                                    handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 
                                else:
                                    ax[row,col].legend(loc=figParams['legendLoc'],frameon=False, bbox_to_anchor=figParams['legendBboxToAnchor'], ncol = 1,fontsize=figParams['traceSet_fontSize']-2,\
                                                    handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 
                                ax[row,col].set_xlim(xlim_s)
                                if not figParams['clean']:
                                    ax[row,col].set_xticks(xticks_s)
                                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['traceSet_fontSize']-2,'verticalalignment': 'top','horizontalalignment': 'center'})
                                else:

                                    ax[row,col].set_yticks([])
                                    ax[row,col].set_yticklabels([])
                                    ax[row,col].set_xticks([])
                                    ax[row,col].spines['bottom'].set_visible(False)
                                    ax[row,col].spines['left'].set_visible(False)
                                # print(f'setting ylims for row {row} col {col} to minVal {minVal:.3f} and maxVal {maxVal:.3f}')
                                ax[row,col].set_ylim([minVal,maxVal])
                                if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                    ax[row,col].set_xlabel(groupLabel+"(s)",fontsize=figParams['traceSet_fontSize'],color=(0,0,0))
                                else:
                                    ax[row,col].set_xlabel(groupLabel+"(s; B"+str(int(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key]))+")",\
                                        fontsize=figParams['traceSet_fontSize'],color=(0,0,0))

                                if not figParams['clean']:
                                    if ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key] == 1:
                                        ax[row,col],figParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['plot_scaleBar'],(0,0,0),True,\
                                            figParams['traceSet_export_align_data_scaleBar_label'][a]+"\n"+figParams['plot_data_axis_label'])
                                    else:
                                        ax[row,col],figParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['plot_scaleBar'],(0,0,0),True,\
                                            figParams['traceSet_export_align_data_scaleBar_label'][a]+\
                                            "\n(bin"+str(int(ROI_clustering[ROI_type_key][ROI_score]['allParams'][figParams['traceSet_group']]['all_traceAlignParams'][align_data]['traceAlignParams']['align_reduce_factor'][ROI_type_key]))+")\n"+\
                                                figParams['plot_data_axis_label'])
                                else:
                                    ax[row,col],figParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['plot_scaleBar'],(0,0,0),True,\
                                        figParams['traceSet_export_align_data_scaleBar_label'][a])
                                
                                ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['traceSet_fontSize']-2)

                                ####################################################################
                                pos1 = ax[row,col].get_position()
                                # print(pos1.height)
                                ax[row,col].set_position([pos1.x0-(1/nrealCols)*figParams['widthAdjust']*realCol, pos1.y0+pos1.height*(1-figParams['trace_heightAdjust']), \
                                    pos1.width+(1/nrealCols)*figParams['widthAdjust'], pos1.height*figParams['trace_heightAdjust']])
                                pos1 = ax[row,col].get_position()
                                # print(pos1.height)
            

        # if figParams['clean']:
        #     pos = ax[row,col].get_position()
        #     print(f'pos for row {row} col {col} relwidth {pos.width:.3f} width {figsize[0]*pos.width:.3f} relheight {pos.height:.3f} height {figsize[1]*pos.height:.3f}')
        ################################################################################################################################################################################################
        #Spacer
        row = 0
        col = -1
        realCol = -1
        for a,align_data in enumerate(figParams['traceSet_export_align_data']):
            ##################################
            align_data_short = copy.deepcopy(align_data)
            align_data_short = align_data_short.replace("allTC_","")
            ##################################
            for rt,ROI_type_key in enumerate(figParams['export_ROI_type_keys']):
                for ROI_score in figParams['export_ROI_scores']:
                    for cg in range(len(figParams['export_traceSets'].keys())):
                        col+=2
                        realCol+=1
                        realRow = -1
                        if figParams['traceSet_overlay']:
                            realRow+=1
                        for ts in range(nTS):
                            if not figParams['traceSet_overlay']:
                                realRow+=1
                            traceSet = list(figParams['export_traceSets'][cg].keys())[ts]
                            if ((figParams['traceSet_overlay'] and ts == nTS-1) or not figParams['traceSet_overlay']):
                                if figParams['traceSet_split_ttypes']:
                                    if figParams['traceSet_overlay']:
                                        row = subRows - 1
                                    else:
                                        row = (ts*subRows)+s + 1
                                    ttype = figParams['export_traceSets'][cg][traceSet]['ttypes'][s]
                                    if 'ttypeLabels' in figParams['export_traceSets'][cg][traceSet].keys():
                                        ttypeLabel = figParams['export_traceSets'][cg][traceSet]['ttypeLabels'][s]
                                    else:
                                        ttypeLabel = ttype
                                    fig.delaxes(ax[row,col])
                                else:
                                    if figParams['traceSet_overlay']:
                                        row = 1
                                    else:
                                        row = (ts*subRows)+1
                                    fig.delaxes(ax[row,col])
        ################################################################################################################################################################################################
        
        # for row in range(nRows):
        #     for col in range(nCols):
        #         if not figParams['manColWidth'] == 0:
        #             pos = ax[row,col].get_position()
        #             print(f'OLD pos for row {row} col {col} relwidth {pos.width:.3f} width {figsize[0]*pos.width:.3f} relheight {pos.height:.3f} height {figsize[1]*pos.height:.3f}')
        #             # ax[row,col].set_position([pos.x0, pos.y0, figsize[0]/figParams['manColWidth'], pos.height])
        #             # print(f'NEW pos for row {row} col {col} width {figsize[0]*pos.width:.3f} height {figsize[1]*pos.height:.3f}')


        
        if delRow:
            for c in range(nCols):
                fig.delaxes(ax[1,c])
        if delCol:
            for r in range(nRows):
                fig.delaxes(ax[r,1])
        if figParams['clean']:
            for row in range(nRows):
                for col in range(nCols):
                    try:
                        pos = ax[row,col].get_position()
                        # print(f'pos for row {row} col {col} relwidth {pos.width:.3f} width {figsize[0]*pos.width:.3f} relheight {pos.height:.3f} height {figsize[1]*pos.height:.3f}')
                    except:
                        pass

        if not 'suptitle_margin' in figParams:
            figParams['suptitle_margin'] = 0.1
        # set_dynamic_suptitle(fig,batchID+" | "+figName+"\n"+str(figDir)+"\n\n\n\n",margin = figParams['suptitle_margin'],\
        #             fontsize=figParams['traceSet_fontSize']+2,color=(0,0,0),va='bottom')
        pdf.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a pdf page
        if close_all_figs:
            plt.close()
        else:
            fig.set_facecolor((1,1,1))
            fig.patch.set_alpha(1)
            plt.show(fig);      
    print("Finshed!")
    return plot_contrasts, im_contrasts

def prep_manual_contrasts(figParams,export_groups,export_trialGroups,export_align_data):
    if not 'manual_im_contrasts' in figParams:
        figParams['manual_im_contrasts'] = {}
    for group in export_groups:
        if not group in figParams['manual_im_contrasts']:
            figParams['manual_im_contrasts'][group] = {}
        for trialGrouping in export_trialGroups:
            if not trialGrouping in figParams['manual_im_contrasts'][group]:
                figParams['manual_im_contrasts'][group][trialGrouping] = {}
            for align_data in export_align_data:
                if not align_data in figParams['manual_im_contrasts'][group][trialGrouping]:
                    figParams['manual_im_contrasts'][group][trialGrouping][align_data] = {}
                if not 'im_minVal' in figParams['manual_im_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_im_contrasts'][group][trialGrouping][align_data]['im_minVal'] = -1e6
                if not 'im_maxVal' in figParams['manual_im_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_im_contrasts'][group][trialGrouping][align_data]['im_maxVal'] = 1e6
                if not 'im_minVal_cont' in figParams['manual_im_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_im_contrasts'][group][trialGrouping][align_data]['im_minVal_cont'] = 0
                if not 'im_maxVal_cont' in figParams['manual_im_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_im_contrasts'][group][trialGrouping][align_data]['im_maxVal_cont'] = 1
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    if not ROI_type_key in figParams['manual_im_contrasts'][group][trialGrouping][align_data].keys():
                        figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key] = {}
                    for ROI_score in figParams['export_ROI_scores']:
                        if not ROI_score in figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key].keys():
                            figParams['manual_im_contrasts'][group][trialGrouping][align_data][ROI_type_key][ROI_score] = {}
    if not 'manual_plot_contrasts' in figParams:
        figParams['manual_plot_contrasts'] = {}
    for group in export_groups:
        if not group in figParams['manual_plot_contrasts']:
            figParams['manual_plot_contrasts'][group] = {}
        for trialGrouping in export_trialGroups:
            if not trialGrouping in figParams['manual_plot_contrasts'][group]:
                figParams['manual_plot_contrasts'][group][trialGrouping] = {}
            for align_data in export_align_data:
                if not align_data in figParams['manual_plot_contrasts'][group][trialGrouping]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data] = {}
                if not 'plot_shift_init' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_shift_init'] = 0
                if not 'plot_shift_final' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_shift_final']=-1e6
                if not 'plot_maxVal' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_maxVal']=1e6
                if not 'plot_minVal' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_minVal'] = 0
                if not 'plot_maxVal_final' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_maxVal_final'] = 0
                if not 'plot_minVal_final' in figParams['manual_plot_contrasts'][group][trialGrouping][align_data]:
                    figParams['manual_plot_contrasts'][group][trialGrouping][align_data]['plot_minVal_final'] = 1
                for ROI_type_key in figParams['export_ROI_type_keys']:
                    if not ROI_type_key in figParams['manual_plot_contrasts'][group][trialGrouping][align_data].keys():
                        figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key] = {}
                    for ROI_score in figParams['export_ROI_scores']:
                        if not ROI_score in figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key].keys():
                            figParams['manual_plot_contrasts'][group][trialGrouping][align_data][ROI_type_key][ROI_score] = {}

    return figParams

def imshow_cleanup(ax,fullClean=True,faceClean=True):
    if fullClean:
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
    if faceClean:
        ax.set_facecolor('none')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    return ax
##################################################################################################

#Trial Alignment by ROI pages
def trial_alignments_byROI_pages(mergePDF, summaryInfo,\
    traceAlignParams, aligned_traces, figParams, trace_type, align_data, zoom, all_ttypes, ROI_score, ROI_type_key, \
        anmIdx, anm, ROI_idx, page_count, cluster, cluster_idx, saveFigs=True,verbose=False,nPreviewPages = 1):

    if not 'clean' in figParams:
        figParams['clean'] = False
    if not 'facecolor' in figParams:
        figParams['facecolor'] = 'none'
    if not 'constrained_layout' in figParams:
        figParams['constrained_layout'] = False  
    if not 'legendReversed' in figParams.keys():
        figParams['legendReversed'] = False
    if not 'legendLoc' in figParams.keys():
        figParams['legendLoc'] = 'upper right'
    if not 'legendBboxToAnchor' in figParams.keys():
        figParams['legendBboxToAnchor'] = (1, 1)
    if not 'max_nTrials_display' in figParams.keys():
        figParams['max_nTrials_display'] = 0
    if not 'legendFontsize' in figParams.keys():
        figParams['legendFontsize'] = figParams['scale_font']
    if not 'lineWidth' in figParams.keys():
        figParams['lineWidth'] = 1
    if not 'filt_trace_byTrial' in figParams.keys():
        figParams['filt_trace_byTrial'] = False
    if not 'splitPrePost' in figParams.keys():
        figParams['splitPrePost'] = True
    if not 'im_horzCueLocation' in figParams:
        figParams['im_horzCueLocation'] = 'top'
    if not 'plot_horzCueLocation' in figParams:
        figParams['plot_horzCueLocation'] = 'top'

    if 'lick' in traceAlignParams['align_data'] or 'Lick' in traceAlignParams['align_data']:
        lickMode = True
    else:
        lickMode = False
    if not 'match_trial_counts' in figParams:
        figParams['match_trial_counts'] = False
    if not 'formatMode' in figParams:
        figParams['formatMode'] = 'horz'
    sess_trace_type = 'bySess_byTtype'
    if 'allSess_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap,sess_colors,_=generate_cmap(figParams['sess_cmap'],len(curr_sessions)+1)
    elif 'bySess_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap,sess_colors,_=generate_cmap(figParams['sess_cmap'],len(curr_sessions)+1)
    elif 'byBehavSess_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap,sess_colors,_=generate_cmap(figParams['sess_cmap'],len(curr_sessions)+1)
    elif 'byBehavEpoch_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap = []
        sess_colors = [[0,0.5,1],[1,0.5,0],[1,0,0.5]]
    elif 'byPrePost_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap = []
        sess_colors = [[0,0.5,1],[1,0.5,0]]
        sess_colors = [[0.3,0.3,0.3],[0,0,0]]
    elif 'preShiftOnly_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap = []
        sess_colors = [[0.5,0.5,0.5]]
    elif 'postShiftOnly_byTtype' == trace_type:
        sess_trace_type = 'bySess_byTtype'
        curr_sessions = list(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type].keys())
        sess_cmap = []
        sess_colors = [[0.5,0.5,0.5]]


    if len(all_ttypes)==0:
        all_ttypes = list(traceAlignParams['align_ttype_codes'])
    alignLabel = copy.deepcopy(traceAlignParams['align_label'])
    alignLabel1 = copy.deepcopy(traceAlignParams['align_label'])

    if not figParams['clean']:

        if "NMF Decon. Spikes Events" in alignLabel:
            alignLabel = alignLabel.replace("NMF Decon. Spikes Events", "Decon. Events")
            alignLabel1 = alignLabel1.replace("NMF Decon. Spikes Events", "Decon. Events")
        if " (SD Norm." in alignLabel:
            alignLabel = alignLabel.replace(" (SD Norm", " \n(SD Norm")
            alignLabel1 = alignLabel1.replace(" (SD Norm", " \n(SD Norm")
        if " (z-score" in alignLabel:
            alignLabel = alignLabel.replace(" (z-score", " \n(z-score")
            alignLabel1 = alignLabel1.replace(" (z-score", " \n(z-score")
        if 'Dx2,Sx1' in alignLabel:
            alignLabel = alignLabel.replace('Dx2,Sx1', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
            alignLabel1 = alignLabel1.replace('Dx2,Sx1', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
        if 'Dx4,Sx2' in alignLabel:
            alignLabel = alignLabel.replace('Dx4,Sx2', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
            alignLabel1 = alignLabel1.replace('Dx4,Sx2', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
        if 'Dx8,Sx4' in alignLabel:
            alignLabel = alignLabel.replace('Dx8,Sx4', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
            alignLabel1 = alignLabel1.replace('Dx8,Sx4', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
        if 'Dx10,Sx5' in alignLabel:
            alignLabel = alignLabel.replace('Dx10,Sx5', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
            alignLabel1 = alignLabel1.replace('Dx10,Sx5', "Bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key]))
    else:
        if "NMF Decon. Spikes Events" in alignLabel:
            alignLabel = alignLabel.replace("NMF Decon. Spikes Events", "Imp")
            alignLabel1 = alignLabel1.replace("NMF Decon. Spikes Events", "Imp")
        if " (SD Norm." in alignLabel:
            alignLabel = alignLabel.replace(" (SD Norm", "")
            alignLabel1 = alignLabel1.replace(" (SD Norm", "")
        if " (z-score" in alignLabel:
            alignLabel = alignLabel.replace(" (z-score", "")
            alignLabel1 = alignLabel1.replace(" (z-score", "")
        if 'Dx2,Sx1' in alignLabel:
            alignLabel = alignLabel.replace('Dx2,Sx1', "")
            alignLabel1 = alignLabel1.replace('Dx2,Sx1', "")
        if 'Dx4,Sx2' in alignLabel:
            alignLabel = alignLabel.replace('Dx4,Sx2', "")
            alignLabel1 = alignLabel1.replace('Dx4,Sx2', "")
        if 'Dx8,Sx4' in alignLabel:
            alignLabel = alignLabel.replace('Dx8,Sx4', "")
            alignLabel1 = alignLabel1.replace('Dx8,Sx4', "")
        if 'Dx10,Sx5' in alignLabel:
            alignLabel = alignLabel.replace('Dx10,Sx5', "")
            alignLabel1 = alignLabel1.replace('Dx10,Sx5', "")
        alignLabel = alignLabel.replace('; ', "")
        alignLabel1 = alignLabel1.replace('; ', "")
        alignLabel = alignLabel.replace(')', "")
        alignLabel1 = alignLabel1.replace(')', "")

    # count = 0
    if not figParams['clean']:
        if figParams['plot_error']:
            alignLabel = alignLabel + "\n"+figParams['plot_key_label']

        figParams['vert_scale_label']=alignLabel+" "
        figParams['vert_scale_label'] = figParams['vert_scale_label'].replace(". ",".\n")
    else:
        figParams['vert_scale_label']="Mean "+alignLabel1
    ###################################################################################################################################################################################################
    if 'allSess_byTtype' in trace_type:
        if verbose:
            print("Adding anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
            " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count),end='...')
        else:
            print('',end='.')
        #################################################################
        if figParams['formatMode'] == 'horz':
            nRows = 2
            nCols = len(all_ttypes)
        else:
            nRows = len(all_ttypes) + 1
            nCols = 2
        if zoom:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im_zoom'],nRows*figParams['vertScalar_im_zoom']))
        else:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im'],nRows*figParams['vertScalar_im']))
        #################################################################
        row = 0
        col = 0
        maxVal=-1e6
        minVal=1e6
        max_nFr=0
        max_s=0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for sess in curr_sessions:
                nTrials = nTrials+aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['align_framerate']])
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        maxVal=np.nanmax([maxVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['data'],figParams['maxPer'])])
                        minVal=np.nanmin([minVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['data'],figParams['minPer'])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])


        if minVal>=maxVal:
            minVal=0
            maxVal=1
        maxCont=np.round(maxVal*figParams['highROIImCont'][ROI_type_key],decimals=4)
        minCont=np.round(minVal*figParams['lowROIImCont'][ROI_type_key],decimals=4)
        if minCont>=maxCont:
            minCont=minVal
            maxCont=maxVal
        if figParams['manual_imCont']:
            final_maxCont=figParams['manual_max_imCont'][ROI_type_key]
            final_minCont=figParams['manual_min_imCont'][ROI_type_key]
        else:
            final_maxCont=maxCont
            final_minCont=minCont
        tempRange = []
        shiftMarkers={}
        for t1,ttype in enumerate(all_ttypes):

            days = []
            days1 = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['dayRelativeToShift']))
            for d in days1:
                days.append(d)

            dayLabels = []
            dayRanges = []
            currentDay = 0
            for day in days:
                n = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['dayRelativeToShift']) == day)
                dayLabels.append("Day "+str(day)+"\n(n = "+str(n)+")")
                dayRanges.append([currentDay,currentDay+n])
                currentDay+=n
        

            days = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['dayRelativeToShift']))
            dayLabels = []
            dayRanges = []
            currentDay = 0
            for day in days:
                n = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['dayRelativeToShift']) == day)
                dayLabels.append("Day "+str(day)+"\n(n = "+str(n)+")")
                dayRanges.append([currentDay,currentDay+n])
                currentDay+=n


            if figParams['formatMode'] == 'horz':
                row = 0
                col = t1
            else:
                row = t1
                col = 0
            shiftMarkers[ttype] = np.nan
            tempRange=np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])+\
                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['align_shift_fr']
            if len(tempRange) > 0:
                tempRangeAdjust = tempRange[np.argmin(np.abs(tempRange))]
                xRange = np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])-tempRangeAdjust
                xRange_labels = xRange + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']

                emptyRow = np.ones((1,max_nFr),dtype='float32')*np.nan
                edges = []
                edgeLabels = []
                allSessData = np.zeros((0,max_nFr),dtype='float32')
                if figParams['addSpacers']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                    edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                else:
                    edges = [0]
                for sess in curr_sessions:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['data'].shape[0] > 0:
                        if 0 in aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['dayRelativeToShift']:
                            shiftMarkers[ttype] = allSessData.shape[0] - 1 + \
                                np.argmin(np.absolute(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['imaging_trial_idxs'])-figParams['shiftDay_nTrials']))
                        else:
                            if np.isnan(shiftMarkers[ttype]) and np.any(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['dayRelativeToShift'])>0):
                                shiftMarkers[ttype] = allSessData.shape[0] - 1
                        tempTrace = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['data'])
                        if figParams['revert_mask']:
                            tempTrace[np.isnan(tempTrace)] = 0
                        allSessData = np.concatenate((allSessData,tempTrace),axis=0)
                        if figParams['addSpacers']:
                            for n in range(figParams['nSpacers']):
                                allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                        edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                        edgeLabels.append(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']))
                if figParams['addSpacers'] and figParams['endSpacer']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                if figParams['nan_empty_trials']:
                    for i in range(allSessData.shape[0]):
                        if np.sum(allSessData[i,:] == 0) == allSessData.shape[1]:
                            allSessData[i,:] = np.nan
                if lickMode:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['lickColors'][ROI_idx],figParams['colorScalar'])
                else:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['cmap'],figParams['colorScalar'])
                cmap.set_bad(color=figParams['nan_color'])
                tempImage = convert2RGB(allSessData[:,plotIdxs],\
                    final_minCont,final_maxCont,cmap,figParams['colorScalar'],figParams['nan_color'])
                tempImage = tempImage[:,:,0:3]
                ax[row,col].imshow(tempImage,extent=[xlim_s[0],xlim_s[1],allSessData.shape[0],0],interpolation='none',origin = 'lower',clip_on=False,aspect = 'auto')
                ax[row,col].text(xlim_s[0],0,ttype,fontsize=figParams['scale_font']-2,color=figParams['export_ttype_colors'][t1],ha='left',va='bottom')
                for e0,el in enumerate(edgeLabels):
                    e = edges[e0]
                    if figParams['addSpacers']:
                        ax[row,col].axhline(e,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                    else:
                        ax[row,col].axhline(e-0.5,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                for e,label in enumerate(edgeLabels):
                    ax[row,col].text(xlim_s[1],edges[e],label,\
                                    fontsize=figParams['scale_font']-4,color=sess_colors[e],ha='right',va='top')
                if np.isfinite(shiftMarkers[ttype]):
                    ax[row,col].axhline(shiftMarkers[ttype],color=figParams['shiftTrial_color'],linestyle=figParams['shiftTrial_lineStyle'],linewidth=figParams['shiftTrial_lineWidth'],alpha=figParams['shiftTrial_alpha'])
                    ax[row,col].text(xlim_s[0],shiftMarkers[ttype],'Shift',\
                                    fontsize=figParams['scale_font']-4,color=figParams['shiftTrial_color'],ha='left',va='top')

                ax[row,col].set_xlim(xlim_s)
                if figParams['match_trial_counts']:
                    ax[row,col].set_ylim([-3,max_nTrials+2])
                    ax[row,col].set_yticks(trial_ticks)
                    if t1 == 0:
                        ax[row,col].set_yticklabels([str(t) for t in trial_ticks],fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'center','horizontalalignment': 'right'})
                    else:
                        ax[row,col].set_yticklabels([])
                else:
                    ax[row,col].set_ylim([-3,allSessData.shape[0]+2])
                    
                if not figParams['yTicksOn']:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                if figParams['yTicksOn'] and t1==0:
                    ax[row,col].set_ylabel("Trial Count",fontsize=figParams['scale_font']+2)
                if figParams['dayLabels']:
                    for d in range(len(days)):
                        ax[row,col].plot([xlim_s[0]+1/framerate*0.25,xlim_s[0]+1/framerate*0.25],[dayRanges[d][0],dayRanges[d][1]],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                        ax[row,col].text(xlim_s[0],np.mean(dayRanges[d]),dayLabels[d],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),ha='right',va='center',fontsize=figParams['scale_font']-5)

                if (figParams['formatMode'] == 'horz' and t1+1 == len(all_ttypes)) or (figParams['formatMode'] == 'vert' and t1 == 0):
                    if figParams['match_trial_counts']:
                        ax[row,col].text(xlim_s[1],max_nTrials/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                    else:
                        ax[row,col].text(xlim_s[1],allSessData.shape[0]/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                ax[row,col].invert_yaxis()
                if figParams['formatMode'] == 'horz':
                    ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
                    ax[row,col].set_xticks(xticks_s)
                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'top','horizontalalignment': 'center'})
                else:
                    ax[row,col].set_xticks([])
                if t1 == 0:
                    if len(cluster):
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count)+"\n"+cluster+" "+str(cluster_idx),fontsize=figParams['scale_font']+2)
                    else:
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count),fontsize=figParams['scale_font']+2)
                ax[row,col].spines['top'].set_visible(False)
                ax[row,col].spines['bottom'].set_visible(False)
                ax[row,col].spines['left'].set_visible(False)
                ax[row,col].spines['right'].set_visible(False)
                sess = 0
                if t1 == 0:
                    figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],(1,1,1),True,' Trials',str(figParams['plot_scaleBar']['length_s'])+" s")
                if figParams['formatMode'] == 'horz':
                    tempTtype = ttype
                else:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                    True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        0,-2,figParams['horzCueLineOn'],figParams['horzCueImScalar'])

        #################################################################
        row=1
        maxVal=-1e6
        minVal=1e6
        max_nFr = 0
        max_s = 0
        for t1,ttype in enumerate(all_ttypes):
            max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['align_length_fr']])
            max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['align_framerate']])
            if ttype in figParams['autoScale_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                if figParams['plot_error']:
                    if 'boot_CI' in figParams['plot_key_error_pos']:
                        maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                        minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                    else:
                        maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                        minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                else:
                    maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
                    minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        plotShift = 0
        plotmaxVal=-1e6
        plotminVal=1e6
        for t1,ttype in enumerate(all_ttypes):
            if ttype in figParams['autoScale_ttypes']:
                if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats']>0:
                    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                        load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                    plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        if figParams['plot_error']:
                            if 'boot_CI' in figParams['plot_key_error_pos']:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                            np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                            np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                        plotmaxVal=np.nanmax([plotmaxVal,np.nanmax(upperTrace[plotIdxs])])
                        plotminVal=np.nanmin([plotminVal,np.nanmin(lowerTrace[plotIdxs])])
        plotmaxVal = plotmaxVal+(plotmaxVal-plotminVal)*figParams['vertUpperPlotBuffer']
        plotminVal = plotminVal-(plotmaxVal-plotminVal)*figParams['vertLowerPlotBuffer']
        if plotminVal>=plotmaxVal:
            plotminVal=-1
            plotmaxVal=1
            plotShift = 0
        plotminVal = plotminVal-plotShift/2
        if figParams['manual_plotLims']:
            final_plotmaxVal=figParams['manual_plotmaxVal'][ROI_type_key]
            final_plotminVal=figParams['manual_plotminVal'][ROI_type_key]
            final_plotShift=figParams['manual_plotShift'][ROI_type_key]
        else:
            final_plotmaxVal=plotmaxVal
            final_plotminVal=plotminVal
            final_plotShift=plotShift
        vert_scale=smart_round(figParams['plot_scaleBar']['height_per']*(final_plotmaxVal-final_plotminVal))
        horz_scale_coord=[np.round(max_s)+traceAlignParams['align_shift_s'],final_plotmaxVal*0.95]
        for t11,ttype in enumerate(reversed(all_ttypes)):
            sess = 0
            t1 = len(all_ttypes) - t11 - 1 


            if figParams['formatMode'] == 'horz':
                row = 1
                col = t1
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                    False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            else:
                row = nRows-1
                col = 0
                if t11 == 0:
                    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                        load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                    if zoom:
                        ax[row,col].set_xlim(xlim_s)
                    else:
                        ax[row,col].set_xlim(xlim_s)
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                    trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                        False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                            final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats']>0:
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)

                plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                if figParams['plot_error']:
                    if 'boot_CI' in figParams['plot_key_error_pos']:
                        upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                        lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                    else:
                        upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                        lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                else:
                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                if 'export_ttype_colors' in figParams and final_plotShift == 0:
                    tempColor = list(figParams['export_ttype_colors'][t1])
                else:
                    if lickMode:
                        tempColor =  list(figParams['lickColors'][ROI_idx])
                    else:
                        tempColor = list(summaryInfo['ROI_types_colors'][ROI_type_key][ROI_score])
                for c in range(len(tempColor)):
                    if tempColor[c]>0.25:
                        tempColor[c] = tempColor[c]-0.25
                if figParams['plot_error']:
                    ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                                            lowerTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                                                color=tempColor,alpha = figParams['alpha'],lw=0)
                if 'export_ttype_colors' in figParams and final_plotShift == 0:
                    tempColor = list(figParams['export_ttype_colors'][t1])
                else:
                    if lickMode:
                        tempColor =  list(figParams['lickColors'][ROI_idx])
                        tempLabel = ROI_type_key+" "+summaryInfo['lickTypes'][ROI_idx]
                    else:
                        tempColor = list(summaryInfo['ROI_types_colors'][ROI_type_key][ROI_score])
                        tempLabel = summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)
                if lickMode:
                    tempLabel = ROI_type_key+" "+summaryInfo['lickTypes'][ROI_idx]
                else:
                    tempLabel = summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)
                ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                color=tempColor,linewidth=1,\
                                label="All Sess "+ttype+"\n(n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['orig_nRepeats'])+\
                                    " | "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats'])+")")
                if zoom:
                    ax[row,col].set_xlim(xlim_s)
                    text_xPos = xlim_s[0]
                    vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
                    horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]
                else:
                    ax[row,col].set_xlim(xlim_s)
                    text_xPos = xlim_s[0]
                    vert_scale_coord=[xlim_s[-1]-0.2,final_plotmaxVal*0.7]
                    horz_scale_coord=[xlim_s[-1]-0.2,final_plotmaxVal*0.95]
            if figParams['formatMode'] == 'horz':
                text_yPos = -final_plotShift*sess
                if not np.isfinite(text_yPos):
                    print(text_yPos)
                else:
                    ax[row,col].text(text_xPos,text_yPos,\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype]['day_label']+\
                                " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['orig_nRepeats'])+\
                                    " | "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][ttype][figParams['plot_stat']]['nRepeats']),\
                                    fontsize=figParams['scale_font']-4,color=(0,0,0),ha='left',va='top')
            else:
                ax[row,col].legend(frameon = False, fontsize=figParams['legendFontsize'])
            ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
            ax[row,col].set_yticks([])
            ax[row,col].set_yticklabels([])
            ax[row,col].spines['top'].set_visible(False)
            ax[row,col].spines['right'].set_visible(False)
            ax[row,col].spines['left'].set_visible(False)
            
            if not figParams['clean']:
                if traceAlignParams['align_reduce_factor'][ROI_type_key] == 0 or traceAlignParams['align_reduce_factor'][ROI_type_key] == 1:
                    ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s)",fontsize=figParams['scale_font']+2)
                else:
                    ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s; bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key])+")",fontsize=figParams['scale_font']+2)
                ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-2)
            else:

                ax[row,col].set_xticks([])
                ax[row,col].spines['bottom'].set_visible(False)

            if figParams['formatMode'] == 'horz' or (figParams['formatMode'] == 'vert' and t1 == 0):
                ax[row,col]=add_trace_scale_lines(ax[row,col],\
                    horz_scale_coord,figParams['horz_scale'],figParams['horz_scale'],figParams['horz_scale_label'],figParams['horz_scale_loc'],(0,0,0),1,figParams['scale_font']-4,\
                    vert_scale_coord,vert_scale,vert_scale,figParams['vert_scale_label'],figParams['vert_scale_loc'],(0,0,0),1,figParams['scale_font']-4)

        plt.subplots_adjust(wspace=figParams['wspace'], hspace=figParams['hspace'])
        if figParams['formatMode'] == 'horz':
            for t1,ttype in enumerate(all_ttypes):
                row = 1
                col = t1
                tempPos = ax[row,col].get_position()
                ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.45,tempPos.width,tempPos.height*0.4])
        else:
            row = nRows-1
            col = 0
            tempPos = ax[row,col].get_position()
            ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.6,tempPos.width,tempPos.height*0.4])
        if figParams['formatMode'] == 'vert':
            col = 1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
        if saveFigs:
            mergePDF.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a mergePDF page
        if verbose:
            print("Finished!")
        if figParams['preview_figs'] and page_count <= nPreviewPages:
            display_clean_subplots(fig,ax)
        elif figParams['preview_figs'] and not page_count <= nPreviewPages:
            plt.close()
        else:
            if figParams['close_all_figs']:
                plt.close()
            else:
                display_clean_subplots(fig,ax)
    
    elif 'bySess_byTtype' in trace_type:
        raise Exception("Not set up yet for bySess_byTtype")
        if verbose:
            print("Adding anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
                " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count),end='...')
        else:
            print('',end='.')
        #################################################################
        if figParams['formatMode'] == 'horz':
            nRows = 2
            nCols = len(all_ttypes)
        else:
            nRows = len(all_ttypes) + 1
            nCols = 2
        if zoom:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im_zoom'],nRows*figParams['vertScalar_im_zoom']))
        else:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im'],nRows*figParams['vertScalar_im']))
        #################################################################
        maxVal=-1e6
        minVal=1e6
        max_nFr=0
        max_s=0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for sess in curr_sessions:
                nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_framerate']])
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        maxVal=np.nanmax([maxVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'],figParams['maxPer'])])
                        minVal=np.nanmin([minVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'],figParams['minPer'])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])

        if minVal>=maxVal:
            minVal=0
            maxVal=1
        maxCont=np.round(maxVal*figParams['highROIImCont'][ROI_type_key],decimals=4)
        minCont=np.round(minVal*figParams['lowROIImCont'][ROI_type_key],decimals=4)
        if minCont>=maxCont:
            minCont=minVal
            maxCont=maxVal
        if figParams['manual_imCont']:
            final_maxCont=figParams['manual_max_imCont'][ROI_type_key]
            final_minCont=figParams['manual_min_imCont'][ROI_type_key]
        else:
            final_maxCont=maxCont
            final_minCont=minCont
        shiftMarkers = {}
        for t1,ttype in enumerate(all_ttypes):
            if figParams['formatMode'] == 'horz':
                row = 0
                col = t1
            else:
                row = t1
                col = 0
            shiftMarkers[ttype] = np.nan

            for sess in curr_sessions:
                days = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['dayRelativeToShift']))
                dayLabels = []
                dayRanges = []
                currentDay = 0
                for day in days:
                    n = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['dayRelativeToShift']) == day)
                    dayLabels.append("Day "+str(day)+"\n(n = "+str(n)+")")
                    dayRanges.append([currentDay,currentDay+n])
                    currentDay+=n

            tempRange=np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])+\
                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][0][ttype]['align_shift_fr']
            if len(tempRange) > 0:
                tempRangeAdjust = tempRange[np.argmin(np.abs(tempRange))]
                xRange = np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])-tempRangeAdjust
                xRange_labels = xRange + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][0][ttype]['align_shift_fr']

                emptyRow = np.ones((1,max_nFr),dtype='float32')*np.nan
                edges = [][int]
                edgeLabels = []
                allSessData = np.zeros((0,max_nFr),dtype='float32')
                if figParams['addSpacers']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                    edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                else:
                    edges = [0]
                for sess in curr_sessions:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'].shape[0] > 0:
                        tempTrace = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'])
                        if 0 in aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['dayRelativeToShift']:
                            shiftMarkers[ttype] = allSessData.shape[0] - 1 + \
                                np.argmin(np.absolute(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['imaging_trial_idxs'])-figParams['shiftDay_nTrials']))
                        else:
                            if np.isnan(shiftMarkers[ttype]) and np.any(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['dayRelativeToShift'])>0):
                                shiftMarkers[ttype] = allSessData.shape[0] - 1
                        if figParams['revert_mask']:
                            tempTrace[np.isnan(tempTrace)] = 0
                        allSessData = np.concatenate((allSessData,tempTrace),axis=0)
                        if figParams['addSpacers']:
                            for n in range(figParams['nSpacers']):
                                allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                        edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                        edgeLabels.append(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']))
                if figParams['addSpacers'] and figParams['endSpacer']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                if figParams['nan_empty_trials']:
                    for i in range(allSessData.shape[0]):
                        if np.sum(allSessData[i,:] == 0) == allSessData.shape[1]:
                            allSessData[i,:] = np.nan
                if lickMode:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['lickColors'][ROI_idx],figParams['colorScalar'])
                else:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['cmap'],figParams['colorScalar'])
                cmap.set_bad(color=figParams['nan_color'])
                tempImage = convert2RGB(allSessData[:,plotIdxs],\
                    final_minCont,final_maxCont,cmap,figParams['colorScalar'],figParams['nan_color'])
                tempImage = tempImage[:,:,0:3]
                ax[row,col].imshow(tempImage,extent=[xlim_s[0],xlim_s[1],allSessData.shape[0],0],interpolation='none',origin = 'lower',clip_on=False,aspect = 'auto')
                ax[row,col].text(xlim_s[0],0,ttype,fontsize=figParams['scale_font']-2,color=figParams['export_ttype_colors'][t1],ha='left',va='bottom')
                for e0,el in enumerate(edgeLabels):
                    e = edges[e0]
                    if figParams['addSpacers']:
                        ax[row,col].axhline(e,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                    else:
                        ax[row,col].axhline(e-0.5,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                for e,label in enumerate(edgeLabels):
                    ax[row,col].text(xlim_s[1],edges[e],label,\
                                    fontsize=figParams['scale_font']-4,color=sess_colors[e],ha='right',va='top')
                if np.isfinite(shiftMarkers[ttype]):
                    ax[row,col].axhline(shiftMarkers[ttype],color=figParams['shiftTrial_color'],linestyle=figParams['shiftTrial_lineStyle'],linewidth=figParams['shiftTrial_lineWidth'],alpha=figParams['shiftTrial_alpha'])
                    ax[row,col].text(xlim_s[0],shiftMarkers[ttype],'Shift',\
                                    fontsize=figParams['scale_font']-4,color=figParams['shiftTrial_color'],ha='left',va='top')
                ax[row,col].set_xlim(xlim_s)
                if figParams['match_trial_counts']:
                    ax[row,col].set_ylim([-3,max_nTrials+2])
                    ax[row,col].set_yticks(trial_ticks)
                    if t1 == 0:
                        ax[row,col].set_yticklabels([str(t) for t in trial_ticks],fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'center','horizontalalignment': 'right'})
                    else:
                        ax[row,col].set_yticklabels([])
                else:
                    ax[row,col].set_ylim([-3,allSessData.shape[0]+2])
                if not figParams['yTicksOn']:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                if figParams['yTicksOn'] and t1==0:
                    ax[row,col].set_ylabel("Trial Count",fontsize=figParams['scale_font']+2)
                if figParams['dayLabels']:
                    for d in range(len(days)):
                        ax[row,col].plot([xlim_s[0]+1/framerate*0.25,xlim_s[0]+1/framerate*0.25],[dayRanges[d][0],dayRanges[d][1]],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                        ax[row,col].text(xlim_s[0],np.mean(dayRanges[d]),dayLabels[d],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),ha='right',va='center',fontsize=figParams['scale_font']-5)

                if (figParams['formatMode'] == 'horz' and t1+1 == len(all_ttypes)) or (figParams['formatMode'] == 'vert' and t1 == 0):
                    if figParams['match_trial_counts']:
                        ax[row,col].text(xlim_s[1],max_nTrials/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                    else:
                        ax[row,col].text(xlim_s[1],allSessData.shape[0]/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                ax[row,col].invert_yaxis()
                if figParams['formatMode'] == 'horz':
                    ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
                    ax[row,col].set_xticks(xticks_s)
                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'top','horizontalalignment': 'center'})
                else:
                    ax[row,col].set_xticks([])
                if t1 == 0:
                    if len(cluster):
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count)+"\n"+cluster+" "+str(cluster_idx),fontsize=figParams['scale_font']+2)
                    else:
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count),fontsize=figParams['scale_font']+2)
                ax[row,col].spines['top'].set_visible(False)
                ax[row,col].spines['bottom'].set_visible(False)
                ax[row,col].spines['left'].set_visible(False)
                ax[row,col].spines['right'].set_visible(False)
                sess = 0
                if t1 == 0:
                    figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],(1,1,1),True,' Trials',str(figParams['plot_scaleBar']['length_s'])+" s")

                if figParams['formatMode'] == 'horz':
                    tempTtype = ttype
                else:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                    True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        0,-2,figParams['horzCueLineOn'],figParams['horzCueImScalar'])
        #################################################################
        row=1
        maxVal=-1e6
        minVal=1e6
        max_nFr = 0
        max_s = 0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for sess in curr_sessions:
                nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_length_fr']/\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['align_framerate']])

                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim'] and figParams['plot_error']:
                            if 'boot_CI' in figParams['plot_key_error_pos']:
                                maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                                minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                            else:
                                maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                                minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                        else:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        plotShift = (maxVal-minVal)*figParams['plotScalar']
        plotmaxVal=-1e6
        plotminVal=1e6
        for sess in curr_sessions:
            for t1,ttype in enumerate(all_ttypes):
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>0:
                        xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                            load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                        plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                        if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                            if figParams['plot_error']:
                                if 'boot_CI' in figParams['plot_key_error_pos']:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                else:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*sess)
                            plotmaxVal=np.nanmax([plotmaxVal,np.nanmax(upperTrace[plotIdxs])])
                            plotminVal=np.nanmin([plotminVal,np.nanmin(lowerTrace[plotIdxs])])
        plotmaxVal = plotmaxVal+(plotmaxVal-plotminVal)*figParams['vertUpperPlotBuffer']
        plotminVal = plotminVal-(plotmaxVal-plotminVal)*figParams['vertLowerPlotBuffer']
        if plotminVal>=plotmaxVal:
            plotminVal=-1
            plotmaxVal=1
            plotShift = 1/len(curr_sessions)
        plotminVal = plotminVal-plotShift/2
        if figParams['manual_plotLims']:
            final_plotmaxVal=figParams['manual_plotmaxVal'][ROI_type_key]
            final_plotminVal=figParams['manual_plotminVal'][ROI_type_key]
            final_plotShift=figParams['manual_plotShift'][ROI_type_key]
        else:
            final_plotmaxVal=plotmaxVal
            final_plotminVal=plotminVal
            final_plotShift=plotShift
        vert_scale=smart_round(figParams['plot_scaleBar']['height_per']*(final_plotmaxVal-final_plotminVal))
        vert_scale_coord=[np.round(max_s)+traceAlignParams['align_shift_s'],final_plotmaxVal*0.90]
        for t11,ttype in enumerate(reversed(all_ttypes)):
            sess = 0
            t1 = len(all_ttypes) - t11 - 1 
            if figParams['formatMode'] == 'horz':
                row = 1
                col = t1
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                    True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            else:
                row = nRows-1
                col = 0
                if t11 == 0:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                        True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                            final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            for sess in curr_sessions:
                if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']>0:
                    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                        load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                    plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                    if figParams['plot_error']:
                        if 'boot_CI' in figParams['plot_key_error_pos']:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                    else:
                        upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                        lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                    # tempColor = list(summaryInfo['ROI_types_colors'][ROI_type_key][ROI_score])
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(sess_colors[bs])
                    for c in range(len(tempColor)):
                        if tempColor[c]>0.25:
                            tempColor[c] = tempColor[c]-0.25
                    if figParams['plot_error']:
                        ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                                                lowerTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                                                color=tempColor,alpha = figParams['alpha'],lw=0)
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(sess_colors[bs])
                    ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*sess,\
                                    color=tempColor,linewidth=1,\
                                    label="Sess"+str(sess)+" "+ttype+"\n(n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats'])+")")
                    if zoom:
                        ax[row,col].set_xlim(figParams['xlim_zoom'][traceAlignParams['align_position']] )
                        text_xPos = figParams['xlim_zoom'][traceAlignParams['align_position']][0]
                        vert_scale_coord=[figParams['xlim_zoom'][traceAlignParams['align_position']][1]-0.2,final_plotmaxVal*0.7]
                        horz_scale_coord=[figParams['xlim_zoom'][traceAlignParams['align_position']][1]-0.2,final_plotmaxVal*0.95]
                    else:
                        ax[row,col].set_xlim([np.nanmin(xdata),np.nanmax(xdata)])
                        text_xPos = np.nanmin(xdata)
                        vert_scale_coord=[xdata[-1]-0.2,final_plotmaxVal*0.7]
                        horz_scale_coord=[xdata[-1]-0.2,final_plotmaxVal*0.95]
                if figParams['formatMode'] == 'horz':
                    text_yPos = -final_plotShift*sess
                    if not np.isfinite(text_yPos):
                        print(text_yPos)
                    else:
                        ax[row,col].text(text_xPos,text_yPos,\
                                    aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][sess][ttype][figParams['plot_stat']]['nRepeats']),\
                                        fontsize=figParams['scale_font']-4,color=(0,0,0),ha='left',va='top')
                else:
                    ax[row,col].legend(frameon = False, fontsize=figParams['scale_font']-4)
            ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
            ax[row,col].set_yticks([])
            ax[row,col].set_yticklabels([])
            if traceAlignParams['align_reduce_factor'][ROI_type_key] == 0 or traceAlignParams['align_reduce_factor'][ROI_type_key] == 1:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s)",fontsize=figParams['scale_font']+2)
            else:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s; bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key])+")",fontsize=figParams['scale_font']+2)
            ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
            if figParams['formatMode'] == 'horz' or (figParams['formatMode'] == 'vert' and t1 == 0):
                ax[row,col]=add_trace_scale_lines(ax[row,col],\
                    horz_scale_coord,figParams['horz_scale'],figParams['horz_scale'],figParams['horz_scale_label'],figParams['horz_scale_loc'],(0,0,0),1,figParams['scale_font']-4,\
                    vert_scale_coord,vert_scale,vert_scale,figParams['vert_scale_label'],figParams['vert_scale_loc'],(0,0,0),1,figParams['scale_font']-4)
            ax[row,col].spines['top'].set_visible(False)
            ax[row,col].spines['right'].set_visible(False)
            ax[row,col].spines['left'].set_visible(False)
        plt.subplots_adjust(wspace=figParams['wspace'], hspace=figParams['hspace'])
        if figParams['formatMode'] == 'horz':
            for t1,ttype in enumerate(all_ttypes):
                row = 1
                col = t1
                tempPos = ax[row,col].get_position()
                ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.2,tempPos.width,tempPos.height*0.6])
        else:
            row = nRows-1
            col = 0
            tempPos = ax[row,col].get_position()
            ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.2,tempPos.width,tempPos.height*0.6])
        if figParams['formatMode'] == 'vert':
            col = 1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
        if saveFigs:
            mergePDF.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a mergePDF page
        if verbose:
            print("Finished!")
        if figParams['preview_figs'] and page_count <= nPreviewPages:
            display_clean_subplots(fig,ax)
        elif figParams['preview_figs'] and not page_count <= nPreviewPages:
            plt.close()
        else:
            if figParams['close_all_figs']:
                plt.close()
            else:
                display_clean_subplots(fig,ax)

    elif 'byBehavSess_byTtype' in trace_type:
        raise Exception("Not set up yet for byBehavSess_byTtype")
        if verbose:
            print("Adding anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
                " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count),end='...')
        else:
            print('',end='.')
        #################################################################
        if figParams['formatMode'] == 'horz':
            nRows = 2
            nCols = len(all_ttypes)
        else:
            nRows = len(all_ttypes) + 1
            nCols = 2
        if zoom:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im_zoom'],nRows*figParams['vertScalar_im_zoom']))
        else:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im'],nRows*figParams['vertScalar_im']))
        #################################################################
        maxVal=-1e6
        minVal=1e6
        max_nFr=0
        max_s=0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_framerate']])
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        maxVal=np.nanmax([maxVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'],figParams['maxPer'])])
                        minVal=np.nanmin([minVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'],figParams['minPer'])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])

        if minVal>=maxVal:
            minVal=0
            maxVal=1
        maxCont=np.round(maxVal*figParams['highROIImCont'][ROI_type_key],decimals=4)
        minCont=np.round(minVal*figParams['lowROIImCont'][ROI_type_key],decimals=4)
        if minCont>=maxCont:
            minCont=minVal
            maxCont=maxVal
        if figParams['manual_imCont']:
            final_maxCont=figParams['manual_max_imCont'][ROI_type_key]
            final_minCont=figParams['manual_min_imCont'][ROI_type_key]
        else:
            final_maxCont=maxCont
            final_minCont=minCont
        for t1,ttype in enumerate(all_ttypes):
            if figParams['formatMode'] == 'horz':
                row = 0
                col = t1
            else:
                row = t1
                col = 0

            days = []
            for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                days1 = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['dayRelativeToShift']))
                for d in days1:
                    days.append(d)

            dayLabels = []
            dayRanges = []
            currentDay = 0
            for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                for day in days:
                    n = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['dayRelativeToShift']) == day)
                    dayLabels.append("Day "+str(day)+"\n(n = "+str(n)+")")
                    dayRanges.append([currentDay,currentDay+n])
                    currentDay+=n


            tempRange=np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])+\
                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']
            if len(tempRange) > 0:
                tempRangeAdjust = tempRange[np.argmin(np.abs(tempRange))]
                xRange = np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])-tempRangeAdjust
                xRange_labels = xRange + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']

                emptyRow = np.ones((1,max_nFr),dtype='float32')*np.nan
                edges = []
                edgeLabels = []
                allSessData = np.zeros((0,max_nFr),dtype='float32')
                if figParams['addSpacers']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                    edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                else:
                    edges = [0]
                for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'].shape[0] > 0:
                        tempTrace = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'])
                        if figParams['revert_mask']:
                            tempTrace[np.isnan(tempTrace)] = 0
                        allSessData = np.concatenate((allSessData,tempTrace),axis=0)
                        if figParams['addSpacers']:
                            for n in range(figParams['nSpacers']):
                                allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                        edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                        edgeLabels.append(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']))
                if figParams['addSpacers'] and figParams['endSpacer']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                if figParams['nan_empty_trials']:
                    for i in range(allSessData.shape[0]):
                        if np.sum(allSessData[i,:] == 0) == allSessData.shape[1]:
                            allSessData[i,:] = np.nan
                if lickMode:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['lickColors'][ROI_idx],figParams['colorScalar'])
                else:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['cmap'],figParams['colorScalar'])
                cmap.set_bad(color=figParams['nan_color'])
                tempImage = convert2RGB(allSessData[:,plotIdxs],\
                    final_minCont,final_maxCont,cmap,figParams['colorScalar'],figParams['nan_color'])
                tempImage = tempImage[:,:,0:3]
                ax[row,col].imshow(tempImage,extent=[xlim_s[0],xlim_s[1],allSessData.shape[0],0],interpolation='none',origin = 'lower',clip_on=False,aspect = 'auto')
                ax[row,col].text(xlim_s[0],0,ttype,fontsize=figParams['scale_font']-2,color=figParams['export_ttype_colors'][t1],ha='left',va='bottom')
                for e0,el in enumerate(edgeLabels):
                    e = edges[e0]
                    if figParams['addSpacers']:
                        ax[row,col].axhline(e,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                    else:
                        ax[row,col].axhline(e-0.5,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                all_bsess = list(traceAlignParams['animal_sessions'][anm][ROI_type_key]['behavior_sessions'].keys())
                for e,label in enumerate(edgeLabels):
                    ax[row,col].text(xlim_s[1],edges[e],label,\
                                    fontsize=figParams['scale_font']-4,color=sess_colors[e],ha='right',va='top')
                ax[row,col].set_xlim(xlim_s)
                if figParams['match_trial_counts']:
                    ax[row,col].set_ylim([-3,max_nTrials+2])
                    ax[row,col].set_yticks(trial_ticks)
                    if t1 == 0:
                        ax[row,col].set_yticklabels([str(t) for t in trial_ticks],fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'center','horizontalalignment': 'right'})
                    else:
                        ax[row,col].set_yticklabels([])
                else:
                    ax[row,col].set_ylim([-3,allSessData.shape[0]+2])
                if not figParams['yTicksOn']:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                if figParams['yTicksOn'] and t1==0:
                    ax[row,col].set_ylabel("Trial Count",fontsize=figParams['scale_font']+2)
                if figParams['dayLabels']:
                    for d in range(len(days)):
                        ax[row,col].plot([xlim_s[0]+1/framerate*0.25,xlim_s[0]+1/framerate*0.25],[dayRanges[d][0],dayRanges[d][1]],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                        ax[row,col].text(xlim_s[0],np.mean(dayRanges[d]),dayLabels[d],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),ha='right',va='center',fontsize=figParams['scale_font']-5)

                if (figParams['formatMode'] == 'horz' and t1+1 == len(all_ttypes)) or (figParams['formatMode'] == 'vert' and t1 == 0):
                    if figParams['match_trial_counts']:
                        ax[row,col].text(xlim_s[1],max_nTrials/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                    else:
                        ax[row,col].text(xlim_s[1],allSessData.shape[0]/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                ax[row,col].invert_yaxis()
                if figParams['formatMode'] == 'horz':
                    ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
                    ax[row,col].set_xticks(xticks_s)
                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'top','horizontalalignment': 'center'})
                else:
                    ax[row,col].set_xticks([])
                if t1 == 0:
                    if len(cluster):
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count)+"\n"+cluster+" "+str(cluster_idx),fontsize=figParams['scale_font']+2)
                    else:
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count),fontsize=figParams['scale_font']+2)
                ax[row,col].spines['top'].set_visible(False)
                ax[row,col].spines['bottom'].set_visible(False)
                ax[row,col].spines['left'].set_visible(False)
                ax[row,col].spines['right'].set_visible(False)
                sess = 0
                if t1 == 0:
                    # figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s'] * aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['align_framerate']
                    figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],(1,1,1),True,' Trials',str(figParams['plot_scaleBar']['length_s'])+" s")

                if figParams['formatMode'] == 'horz':
                    tempTtype = ttype
                else:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                    True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],0,-2,\
                        figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
        #################################################################
        row=1
        maxVal=-1e6
        minVal=1e6
        max_nFr = 0
        max_s = 0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['align_framerate']])
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim'] and figParams['plot_error']:
                        if 'boot_CI' in figParams['plot_key_error_pos']:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                        else:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                    else:
                        maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
                        minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        plotShift = (maxVal-minVal)*figParams['plotScalar']
        plotmaxVal=-1e6
        plotminVal=1e6
        for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
            for t1,ttype in enumerate(all_ttypes):
                if ttype in figParams['autoScale_ttypes']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>0:
                        xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                            load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                        plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                        if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                            if figParams['plot_error']:
                                if 'boot_CI' in figParams['plot_key_error_pos']:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)
                                else:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*bsess)

                            plotmaxVal=np.nanmax([plotmaxVal,np.nanmax(upperTrace[plotIdxs])])
                            plotminVal=np.nanmin([plotminVal,np.nanmin(lowerTrace[plotIdxs])])
        plotmaxVal = plotmaxVal+(plotmaxVal-plotminVal)*figParams['vertUpperPlotBuffer']
        plotminVal = plotminVal-(plotmaxVal-plotminVal)*figParams['vertLowerPlotBuffer']
        if plotminVal>=plotmaxVal:
            plotminVal=-1
            plotmaxVal=1
            plotShift = 1/len(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys())
        plotminVal = plotminVal-plotShift/2
        if figParams['manual_plotLims']:
            final_plotmaxVal=figParams['manual_plotmaxVal'][ROI_type_key]
            final_plotminVal=figParams['manual_plotminVal'][ROI_type_key]
            final_plotShift=figParams['manual_plotShift'][ROI_type_key]
        else:
            final_plotmaxVal=plotmaxVal
            final_plotminVal=plotminVal
            final_plotShift=plotShift
        vert_scale=smart_round(figParams['plot_scaleBar']['height_per']*(final_plotmaxVal-final_plotminVal))
        vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
        horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]
        for t11,ttype in enumerate(reversed(all_ttypes)):
            sess = 0
            t1 = len(all_ttypes) - t11 - 1 
            if figParams['formatMode'] == 'horz':
                row = 1
                col = t1
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                    False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            else:
                row = nRows-1
                col = 0
                if t11 == 0:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                    trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                        False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                            final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            for bs,bsess in enumerate(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type].keys()):
                if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>0:
                    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                        load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                    plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                    if figParams['plot_error']:
                        if 'boot_CI' in figParams['plot_key_error_pos']:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                    else:
                        upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                        lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(sess_colors[bs])
                    for c in range(len(tempColor)):
                        if tempColor[c]>0.25:
                            tempColor[c] = tempColor[c]-0.25
                    if figParams['plot_error']:
                        ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*bsess,\
                                                                lowerTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*bsess,\
                                                                    color=tempColor,alpha = figParams['alpha'],lw=0)
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(sess_colors[bs])
                    ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*bsess,\
                                    color=tempColor,linewidth=1,\
                                    label="Bsess"+str(bsess)+" "+ttype+"\n(n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats'])+")")
                    ax[row,col].set_xlim(xlim_s)
                    text_xPos = xlim_s[0]
                    vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
                    horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]
                if figParams['formatMode'] == 'horz':
                    text_yPos = -final_plotShift*bsess
                    if not np.isfinite(text_yPos):
                        print(text_yPos)
                    else:
                        ax[row,col].text(text_xPos,text_yPos,\
                                    aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']),\
                                        fontsize=figParams['scale_font']-4,color=(0,0,0),ha='left',va='top')
                else:
                    ax[row,col].legend(frameon = False, fontsize=figParams['scale_font']-4)
            ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
            ax[row,col].set_yticks([])
            ax[row,col].set_yticklabels([])
            if traceAlignParams['align_reduce_factor'][ROI_type_key] == 0 or traceAlignParams['align_reduce_factor'][ROI_type_key] == 1:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s)",fontsize=figParams['scale_font']+2)
            else:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s; bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key])+")",fontsize=figParams['scale_font']+2)
            ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
            if figParams['formatMode'] == 'horz' or (figParams['formatMode'] == 'vert' and t1 == 0):
                ax[row,col]=add_trace_scale_lines(ax[row,col],\
                    horz_scale_coord,figParams['horz_scale'],figParams['horz_scale'],figParams['horz_scale_label'],figParams['horz_scale_loc'],(0,0,0),1,figParams['scale_font']-4,\
                    vert_scale_coord,vert_scale,vert_scale,figParams['vert_scale_label'],figParams['vert_scale_loc'],(0,0,0),1,figParams['scale_font']-4)
            ax[row,col].spines['top'].set_visible(False)
            ax[row,col].spines['right'].set_visible(False)
            ax[row,col].spines['left'].set_visible(False)
        plt.subplots_adjust(wspace=figParams['wspace'], hspace=figParams['hspace'])
        if figParams['formatMode'] == 'horz':
            for t1,ttype in enumerate(all_ttypes):
                row = 1
                col = t1
                tempPos = ax[row,col].get_position()
                ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.5,tempPos.width,tempPos.height*0.4])
        else:
            row = nRows-1
            col = 0
            tempPos = ax[row,col].get_position()
            ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.5,tempPos.width,tempPos.height*0.4])
        if figParams['formatMode'] == 'vert':
            col = 1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
        if saveFigs:
            mergePDF.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a mergePDF page
        if verbose:
            print("Finished!")
        if figParams['preview_figs'] and page_count <= nPreviewPages:
            display_clean_subplots(fig,ax)
        elif figParams['preview_figs'] and not page_count <= nPreviewPages:
            plt.close()
        else:
            if figParams['close_all_figs']:
                plt.close()
            else:
                display_clean_subplots(fig,ax)

    elif 'byBehavEpoch_byTtype' in trace_type:
        raise Exception("byBehavEpoch_byTtype trace type not implemented yet")
        if verbose:
            print("Adding anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
                " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count),end='...')
        else:
            print('',end='.')
        #################################################################
        if figParams['formatMode'] == 'horz':
            nRows = 2
            nCols = len(all_ttypes)
        else:
            nRows = len(all_ttypes) + 1
            nCols = 2
        if zoom:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im_zoom'],nRows*figParams['vertScalar_im_zoom']))
        else:
            fig,ax=clean_subplots(nRows,nCols,figsize=(nCols*figParams['horzScalar_im'],nRows*figParams['vertScalar_im']))
        #################################################################
        maxVal=-1e6
        minVal=1e6
        max_nFr=0
        max_s=0
        max_nTrials = 0
        for t1,ttype in enumerate(all_ttypes):
            nTrials = 0
            for be,behavior_epoch in enumerate(summaryInfo['behavior_epochs']):
                nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0]
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']/\
                                    aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_framerate']])
                if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                    maxVal=np.nanmax([maxVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'],figParams['maxPer'])])
                    minVal=np.nanmin([minVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'],figParams['minPer'])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])


        if minVal>=maxVal:
            minVal=0
            maxVal=1
        maxCont=np.round(maxVal*figParams['highROIImCont'][ROI_type_key],decimals=4)
        minCont=np.round(minVal*figParams['lowROIImCont'][ROI_type_key],decimals=4)
        if minCont>=maxCont:
            minCont=minVal
            maxCont=maxVal
        if figParams['manual_imCont']:
            final_maxCont=figParams['manual_max_imCont'][ROI_type_key]
            final_minCont=figParams['manual_min_imCont'][ROI_type_key]
        else:
            final_maxCont=maxCont
            final_minCont=minCont
        for t1,ttype in enumerate(all_ttypes):
            if figParams['formatMode'] == 'horz':
                row = 0
                col = t1
            else:
                row = t1
                col = 0

            days = []
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    days1 = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['dayRelativeToShift']))
                    for d in days1:
                        days.append(d)

            dayLabels = []
            dayRanges = []
            currentDay = 0
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    for day in days:
                        n = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['dayRelativeToShift']) == day)
                        dayLabels.append("Day "+str(day)+"\n(n = "+str(n)+")")
                        dayRanges.append([currentDay,currentDay+n])
                        currentDay+=n


            tempRange=np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])+\
                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']
            if len(tempRange)>0:
                tempRangeAdjust = tempRange[np.argmin(np.abs(tempRange))]
                xRange = np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])-tempRangeAdjust
                xRange_labels = xRange + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']

                emptyRow = np.ones((1,max_nFr),dtype='float32')*np.nan
                edges = []
                edgeLabels = []
                allSessData = np.zeros((0,max_nFr),dtype='float32')
                if figParams['addSpacers']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                    edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                else:
                    edges = [0]
                for be,behavior_epoch in enumerate(summaryInfo['behavior_epochs']):
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0] > 0:
                        tempTrace = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'])
                        if figParams['revert_mask']:
                            tempTrace[np.isnan(tempTrace)] = 0
                        allSessData = np.concatenate((allSessData,tempTrace),axis=0)
                        if figParams['addSpacers']:
                            for n in range(figParams['nSpacers']):
                                allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                        edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                        edges.append(allSessData.shape[0]-1)
                        edgeLabels.append(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']))
                if figParams['addSpacers'] and figParams['endSpacer']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                if figParams['nan_empty_trials']:
                    for i in range(allSessData.shape[0]):
                        if np.sum(allSessData[i,:] == 0) == allSessData.shape[1]:
                            allSessData[i,:] = np.nan
                if lickMode:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['lickColors'][ROI_idx],figParams['colorScalar'])
                else:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['cmap'],figParams['colorScalar'])
                cmap.set_bad(color=figParams['nan_color'])
                tempImage = convert2RGB(allSessData[:,plotIdxs],\
                    final_minCont,final_maxCont,cmap,figParams['colorScalar'],figParams['nan_color'])
                tempImage = tempImage[:,:,0:3]
                ax[row,col].imshow(tempImage,extent=[xlim_s[0],xlim_s[1],allSessData.shape[0],0],interpolation='none',origin = 'lower',clip_on=False,aspect = 'auto')
                ax[row,col].text(xlim_s[0],0,ttype,fontsize=figParams['scale_font']-2,color=figParams['export_ttype_colors'][t1],ha='left',va='bottom')
                for e0,el in enumerate(edgeLabels):
                    e = edges[e0]
                    if figParams['addSpacers']:
                        ax[row,col].axhline(e,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                    else:
                        ax[row,col].axhline(e-0.5,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                for e,label in enumerate(edgeLabels):
                    ax[row,col].text(xlim_s[1],edges[e],label,\
                                    fontsize=figParams['scale_font']-4,color=sess_colors[e],ha='right',va='top')
                ax[row,col].set_xlim(xlim_s)
                if figParams['match_trial_counts']:
                    ax[row,col].set_ylim([-3,max_nTrials+2])
                    ax[row,col].set_yticks(trial_ticks)
                    if t1 == 0:
                        ax[row,col].set_yticklabels([str(t) for t in trial_ticks],fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'center','horizontalalignment': 'right'})
                    else:
                        ax[row,col].set_yticklabels([])
                else:
                    ax[row,col].set_ylim([-3,allSessData.shape[0]+2])
                if not figParams['yTicksOn']:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                if figParams['yTicksOn'] and t1==0:
                    ax[row,col].set_ylabel("Trial Count",fontsize=figParams['scale_font']+2)
                if figParams['dayLabels']:
                    for d in range(len(days)):
                        ax[row,col].plot([xlim_s[0]+1/framerate*0.25,xlim_s[0]+1/framerate*0.25],[dayRanges[d][0],dayRanges[d][1]],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                        ax[row,col].text(xlim_s[0],np.mean(dayRanges[d]),dayLabels[d],\
                            color=tuple(np.array([0.3,0.3,0.3])+np.mod(d,2)*0.3),ha='right',va='center',fontsize=figParams['scale_font']-5)

                if (figParams['formatMode'] == 'horz' and t1+1 == len(all_ttypes)) or (figParams['formatMode'] == 'vert' and t1 == 0):
                    if figParams['match_trial_counts']:
                        ax[row,col].text(xlim_s[1],max_nTrials/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                    else:
                        ax[row,col].text(xlim_s[1],allSessData.shape[0]/2,"\n\n\n"+alignLabel1+" "+\
                                        str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont)),fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                ax[row,col].invert_yaxis()

                if figParams['formatMode'] == 'horz':
                    tempTtype = ttype
                else:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                    True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        0,-2,figParams['horzCueLineOn'],figParams['horzCueImScalar'])



                if figParams['formatMode'] == 'horz':
                    ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
                    ax[row,col].set_xticks(xticks_s)
                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'top','horizontalalignment': 'center'})
                else:
                    ax[row,col].set_xticks([])
                if t1 == 0:
                    if len(cluster):
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count)+"\n"+cluster+" "+str(cluster_idx),fontsize=figParams['scale_font']+2)
                    else:
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count),fontsize=figParams['scale_font']+2)
                ax[row,col].spines['top'].set_visible(False)
                ax[row,col].spines['bottom'].set_visible(False)
                ax[row,col].spines['left'].set_visible(False)
                ax[row,col].spines['right'].set_visible(False)
                sess = 0
                if t1 == 0:
                    figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],(1,1,1),True,' Trials',str(figParams['plot_scaleBar']['length_s'])+" s")

        #################################################################
        row=1
        maxVal=-1e6
        minVal=1e6
        max_nFr = 0
        max_s = 0
        for be,behavior_epoch in enumerate(summaryInfo['behavior_epochs']):
            for t1,ttype in enumerate(all_ttypes):
                max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']])
                max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']/\
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_framerate']])
                if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                    if figParams['plot_error']:
                        if 'boot_CI' in figParams['plot_key_error_pos']:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                        else:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])])
                    else:
                        maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
                        minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        trial_ticks = np.arange(0,max_nTrials+1,figParams['match_trial_ticks'])
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        plotShift = (maxVal-minVal)*figParams['plotScalar']
        plotmaxVal=-1e6
        plotminVal=1e6
        for be,behavior_epoch in enumerate(summaryInfo['behavior_epochs']):
            for t1,ttype in enumerate(all_ttypes):
                if figParams['formatMode'] == 'horz':
                    row = 1
                    col = t1
                else:
                    row = nRows-1
                    col = 0
                if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][bsess][ttype][figParams['plot_stat']]['nRepeats']>0:
                        xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                            load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                        plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                        if figParams['plot_error']:
                            if 'boot_CI' in figParams['plot_key_error_pos']:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                            np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-\
                                            np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])-plotShift*be)

                        plotmaxVal=np.nanmax([plotmaxVal,np.nanmax(upperTrace[plotIdxs])])
                        plotminVal=np.nanmin([plotminVal,np.nanmin(lowerTrace[plotIdxs])])
        plotmaxVal = plotmaxVal+(plotmaxVal-plotminVal)*figParams['vertUpperPlotBuffer']
        plotminVal = plotminVal-(plotmaxVal-plotminVal)*figParams['vertLowerPlotBuffer']
        if plotminVal>=plotmaxVal:
            plotminVal=-1
            plotmaxVal=1
            plotShift = 1/len(summaryInfo['behavior_epochs'])
        plotminVal = plotminVal-plotShift/2
        if figParams['manual_plotLims']:
            final_plotmaxVal=figParams['manual_plotmaxVal'][ROI_type_key]
            final_plotminVal=figParams['manual_plotminVal'][ROI_type_key]
            final_plotShift=figParams['manual_plotShift'][ROI_type_key]
        else:
            final_plotmaxVal=plotmaxVal
            final_plotminVal=plotminVal
            final_plotShift=plotShift
        vert_scale=smart_round(figParams['plot_scaleBar']['height_per']*(final_plotmaxVal-final_plotminVal))
        vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
        horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]
        for t11,ttype in enumerate(reversed(all_ttypes)):
            sess = 0
            t1 = len(all_ttypes) - t11 - 1 
            if figParams['formatMode'] == 'horz':
                row = 1
                col = t1
                trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                    False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                        final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            else:
                row = nRows-1
                col = 0
                if t11 == 0:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                    trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                        False,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                            final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'])
            for be,behavior_epoch in enumerate(summaryInfo['behavior_epochs']):
                if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']>0:
                    xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                        load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                    plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']])
                    if figParams['plot_error']:
                        if 'boot_CI' in figParams['plot_key_error_pos']:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]+\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_pos']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]-\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_error_neg']]))
                    else:
                        upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                        lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']][figParams['plot_key_main']]))
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(summaryInfo['behavior_epoch_colors'][be])
                    for c in range(len(tempColor)):
                        if tempColor[c]>0.25:
                            tempColor[c] = tempColor[c]-0.25
                    if figParams['plot_error']:
                        ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*be,\
                                                                lowerTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*be,\
                                                                    color=tempColor,alpha = figParams['alpha'],lw=0)
                    if 'export_ttype_colors' in figParams and final_plotShift == 0:
                        tempColor = list(figParams['export_ttype_colors'][t1])
                    else:
                        tempColor = copy.deepcopy(summaryInfo['behavior_epoch_colors'][be])
                    ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift*be,\
                                    color=tempColor,linewidth=1,\
                                    label=behavior_epoch+" "+str(ttype)+"\n(n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats'])+")")
                    ax[row,col].set_xlim(xlim_s)
                    text_xPos = xlim_s[0]
                    vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
                    horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]

                if figParams['formatMode'] == 'horz':
                    text_yPos = -final_plotShift*be
                    if not np.isfinite(text_yPos):
                        print(text_yPos)
                    else:
                        ax[row,col].text(text_xPos,text_yPos,\
                                    aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']+\
                                    " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                        str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']),\
                                        fontsize=figParams['scale_font']-4,color=(0,0,0),ha='left',va='top')
                else:
                    ax[row,col].legend(frameon = False, fontsize=figParams['scale_font']-4)
            ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
            ax[row,col].set_yticks([])
            ax[row,col].set_yticklabels([])
            if traceAlignParams['align_reduce_factor'][ROI_type_key] == 0 or traceAlignParams['align_reduce_factor'][ROI_type_key] == 1:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s)",fontsize=figParams['scale_font']+2)
            else:
                ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s; bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key])+")",fontsize=figParams['scale_font']+2)
            ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
            if figParams['formatMode'] == 'horz' or (figParams['formatMode'] == 'vert' and t1 == 0):
                ax[row,col]=add_trace_scale_lines(ax[row,col],\
                    horz_scale_coord,figParams['horz_scale'],figParams['horz_scale'],figParams['horz_scale_label'],figParams['horz_scale_loc'],(0,0,0),1,figParams['scale_font']-4,\
                    vert_scale_coord,vert_scale,vert_scale,figParams['vert_scale_label'],figParams['vert_scale_loc'],(0,0,0),1,figParams['scale_font']-4)
            ax[row,col].spines['top'].set_visible(False)
            ax[row,col].spines['right'].set_visible(False)
            ax[row,col].spines['left'].set_visible(False)
        plt.subplots_adjust(wspace=figParams['wspace'], hspace=figParams['hspace'])
        if figParams['formatMode'] == 'horz':
            for t1,ttype in enumerate(all_ttypes):
                row = 1
                col = t1
                tempPos = ax[row,col].get_position()
                ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.5,tempPos.width,tempPos.height*0.4])
        else:
            row = nRows-1
            col = 0
            tempPos = ax[row,col].get_position()
            ax[row,col].set_position([tempPos.x0,tempPos.y0+tempPos.height*0.5,tempPos.width,tempPos.height*0.4])
        if figParams['formatMode'] == 'vert':
            col = 1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
        if saveFigs:
            mergePDF.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a mergePDF page
        if verbose:
            print("Finished!")
        if figParams['preview_figs'] and page_count <= nPreviewPages:
            display_clean_subplots(fig,ax)
        elif figParams['preview_figs'] and not page_count <= nPreviewPages:
            plt.close()
        else:
            if figParams['close_all_figs']:
                plt.close()
            else:
                display_clean_subplots(fig,ax)

    elif 'byPrePost_byTtype' in trace_type or 'preShiftOnly_byTtype' in trace_type or 'postShiftOnly_byTtype' in trace_type:
        if 'byPrePost_byTtype' in trace_type:
            behavior_epochs = summaryInfo['byPrePost']
            trace_type = 'byPrePost_byTtype'
            if figParams['splitPrePost']:
                nPlots = 2
            else:
                nPlots = 1
        elif  'preShiftOnly_byTtype' in trace_type:
            behavior_epochs = ['preShift']
            nPlots = 1
            trace_type = 'byPrePost_byTtype'
        elif 'postShiftOnly_byTtype' in trace_type:
            behavior_epochs = ['postShift']
            trace_type = 'byPrePost_byTtype'
            nPlots = 1
        # print(behavior_epochs)
        # print(trace_type)
        if verbose:
            print("Adding anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
                " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count),end='...')
        else:
            print('',end='.')
        #################################################################
        if figParams['formatMode'] == 'horz':
            nRows = 2
            nCols = len(all_ttypes) + nPlots
        else:
            nRows = len(all_ttypes) + nPlots
            nCols = 2
        if zoom:
            figSize=(nCols*figParams['horzScalar_im_zoom'],nRows*figParams['vertScalar_im_zoom'])
        else:
            figSize=(nCols*figParams['horzScalar_im'],nRows*figParams['vertScalar_im'])

        
        fig,ax=clean_subplots(nRows,nCols,figsize=figSize,facecolor = figParams['facecolor'], constrained_layout = figParams['constrained_layout'])
        if not figParams['constrained_layout']:
            plt.subplots_adjust(wspace=figParams['wspace'], hspace=figParams['hspace'])
        if verbose:
            print(f'nRows: {nRows}, nCols: {nCols}, nPlots {nPlots}, figSize: {figSize}')
        added_scaleBar=np.zeros((nRows,nCols),dtype='bool')

        #################################################################
        #imgs/rasters
        maxVal=-1e6
        minVal=1e6
        max_nFr=0
        max_s=0
        max_nTrials = 0
        nTrial_byTtype = {}
        for t1,ttype in enumerate(all_ttypes):
            nTrial_byTtype[ttype] = 0
            nTrials = 0
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    nTrials = nTrials + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0]
                    nTrial_byTtype[ttype] = nTrial_byTtype[ttype] + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0]
                    max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']])
                    max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']/\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_framerate']])
                    if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']>=figParams['minRepeats_forYLim']:
                        maxVal=np.nanmax([maxVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'],figParams['maxPer'])])
                        minVal=np.nanmin([minVal,np.nanpercentile(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'],figParams['minPer'])])
            max_nTrials = np.nanmax([max_nTrials,nTrials])
        max_nTrials = int(np.ceil(max_nTrials/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        if figParams['max_nTrials_display']> 0:
            max_nTrials_display = figParams['max_nTrials_display']
        else:
            max_nTrials_display = max_nTrials
        
        trial_ticks = np.arange(0,max_nTrials_display+1,figParams['match_trial_ticks'])
        
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        maxCont=np.round(maxVal*figParams['highROIImCont'][ROI_type_key],decimals=4)
        minCont=np.round(minVal*figParams['lowROIImCont'][ROI_type_key],decimals=4)
        if minCont>=maxCont:
            minCont=minVal
            maxCont=maxVal
        if figParams['manual_imCont']:
            final_maxCont=figParams['manual_max_imCont'][ROI_type_key]
            final_minCont=figParams['manual_min_imCont'][ROI_type_key]
        else:
            final_maxCont=maxCont
            final_minCont=minCont
        for t1,ttype in enumerate(all_ttypes):
            if figParams['formatMode'] == 'horz':
                row = 0
                col = t1
            else:
                row = t1
                col = 0


            days = {}
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    days[behavior_epoch] = []
                    days1 = list(np.unique(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['dayRelativeToShift']))

                    # print(behavior_epoch+" "+str(days1))
                    for d in days1:
                        days[behavior_epoch].append(d)

            dayLabels = {}
            dayRanges = {}
            currentDay = 0
            if figParams['addSpacers']:
                for n in range(figParams['nSpacers']):
                    currentDay+=1
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    dayLabels[behavior_epoch] = {}
                    dayRanges[behavior_epoch] = {}
                    for day in days[behavior_epoch]:
                        dayLabels[behavior_epoch][day] = []
                        dayRanges[behavior_epoch][day] = []
                        nTr = np.sum(np.array(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['dayRelativeToShift']) == day)
                        if figParams['clean']:
                            dayLabels[behavior_epoch][day] = "D"+str(day)
                        else:
                            dayLabels[behavior_epoch][day] = "Day "+str(day)+"\n(n = "+str(nTr)+")"
                        dayRanges[behavior_epoch][day] = [currentDay,currentDay+nTr]
                        currentDay+=nTr
                        # print(behavior_epoch+" Day"+str(day)+" "+str([dayRanges[behavior_epoch][day][0],dayRanges[behavior_epoch][day][1]]))
                    if figParams['addSpacers']:
                        for n in range(figParams['nSpacers']):
                            currentDay+=1

            tempRange=np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])+\
                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']
            if len(tempRange)>0:
                tempRangeAdjust = tempRange[np.argmin(np.abs(tempRange))]
                xRange = np.arange(0,np.round((max_nFr)/figParams['tickSpacing_fr'])*figParams['tickSpacing_fr'],figParams['tickSpacing_fr'])-tempRangeAdjust
                xRange_labels = xRange + aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][0][ttype]['align_shift_fr']

                emptyRow = np.ones((1,max_nFr),dtype='float32')*np.nan
                edges = [0]
                edgeLabels = []
                edges = []
                edgeLabels = []
                allSessData = np.zeros((0,max_nFr),dtype='float32')
                if figParams['addSpacers']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0) 
                    edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                else:
                    edges = [0]
                for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                    if behavior_epoch in behavior_epochs:
                        if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0] > 0:
                            trialTraces = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'])
                            if figParams['filt_trace_byTrial']:
                                for t in range(trialTraces.shape[0]):
                                    trialTraces[t,:] = running_mean_convolve(trialTraces[t,:],figParams['filt_trace_byTrial_window'][align_data][ROI_type_key])
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data_filt'] = copy.deepcopy(trialTraces)
                            if figParams['revert_mask']:
                                trialTraces[np.isnan(trialTraces)] = 0
                            allSessData = np.concatenate((allSessData,trialTraces),axis=0)
                            if figParams['addSpacers']:
                                for n in range(figParams['nSpacers']):
                                    allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                            edges.append(allSessData.shape[0]-1-np.floor(figParams['nSpacers']/2))
                            if figParams['clean']:
                                if 'pre' in aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']:
                                    bel = "Pre "
                                elif 'post' in aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']:
                                    bel = "Post "
                                else:
                                    bel = aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']+" "
                                edgeLabels.append(bel+\
                                            " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['orig_nRepeats']))
                            else:
                                edgeLabels.append(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']+\
                                            " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['orig_nRepeats'])+" | "+\
                                                str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']]['nRepeats']))
                            if figParams['filt_trace_byTrial']:
                                print("Calculating summary stats for filtered traces for anmIdx "+str(anmIdx)+" "+anm+" "+ROI_type_key+" ROI_score "+str(ROI_score)+\
                                    " ROI"+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" ROI_idx "+str(ROI_idx)+" pg"+str(page_count)+" behavior_epoch "+behavior_epoch+" ttype "+ttype,end='...')
                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][figParams['plot_stat']+"_filt"] = \
                                    simple_trial_trace_summary_stats(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data_filt'],\
                                    additional_stats = traceAlignParams['additional_stats'], bootstrap = traceAlignParams['bootstrap'], nBoot = traceAlignParams['nBoot'], \
                                    CI = traceAlignParams['CI'], save_all_bstraps = traceAlignParams['save_all_bstraps'], \
                                    exclude_naned_trials = traceAlignParams['exclude_naned_trials'], exclusion_nan_percentage = traceAlignParams['exclusion_nan_percentage'],verbose=False,\
                                    stabilityCheck = traceAlignParams['stabilityCheck'], stabilityCheck_nRepeats = traceAlignParams['stabilityCheck_nRepeats'],\
                                        stabilityCheck_ratio = traceAlignParams['stabilityCheck_ratio'],revert_NaNs = traceAlignParams['revert_NaNs'])
                                print("Finished!")

                
                if figParams['addSpacers'] and figParams['endSpacer']:
                    for n in range(figParams['nSpacers']):
                        allSessData = np.concatenate((allSessData,emptyRow),axis=0)
                xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                    load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)

                if figParams['nan_empty_trials']:
                    for i in range(allSessData.shape[0]):
                        if np.sum(allSessData[i,:] == 0) == allSessData.shape[1]:
                            allSessData[i,:] = np.nan
                if lickMode:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['lickColors'][ROI_idx],figParams['colorScalar'])
                else:
                    cmap,tempRGB,tempBGR=generate_cmap(figParams['cmap'],figParams['colorScalar'])
                cmap.set_bad(color=figParams['nan_color'])
                tempImage = convert2RGB(allSessData[:,plotIdxs],\
                    final_minCont,final_maxCont,cmap,figParams['colorScalar'],figParams['nan_color'])
                tempImage = tempImage[:,:,0:3]
                extent=[xlim_s[0],xlim_s[1],0,allSessData.shape[0]]
                # extent=[xlim_s[0],xlim_s[1],0,allSessData.shape[0]]
                
                ax[row,col].imshow(tempImage,extent=extent,interpolation='none',origin = 'lower',clip_on=False,aspect = 'auto')
                # ax[row,col].imshow(allSessData,extent=[xlim_s[0],xlim_s[1],allSessData.shape[0],0],interpolation='none',vmin=final_minCont,vmax=final_maxCont,cmap=cmap,origin = 'lower',clip_on=False,aspect = 'auto')
                ax[row,col].text(xlim_s[0],0,ttype,fontsize=figParams['scale_font']-2,color=figParams['export_ttype_colors'][t1],ha='left',va='bottom')
                for e0,el in enumerate(edgeLabels):
                    e = edges[e0]
                    if figParams['addSpacers']:
                        ax[row,col].axhline(e,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                    else:
                        ax[row,col].axhline(e-0.5,color=sess_colors[e0],linewidth=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                for e,label in enumerate(edgeLabels):
                    ax[row,col].text(xlim_s[1],edges[e],label,\
                                    fontsize=figParams['scale_font']-4,color=sess_colors[e],ha='right',va='top')
                ax[row,col].set_xlim(xlim_s)
                
                if figParams['dayLabels']:
                    for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                        if behavior_epoch in behavior_epochs:
                            for d0,d in enumerate(list(dayRanges[behavior_epoch].keys())):
                                # print(behavior_epoch+" Day"+str(d)+" "+str([dayRanges[behavior_epoch][d][0],dayRanges[behavior_epoch][d][1]]))
                                ax[row,col].plot([xlim_s[0]+1/framerate*0.25,xlim_s[0]+1/framerate*0.25],[dayRanges[behavior_epoch][d][0],dayRanges[behavior_epoch][d][1]],\
                                    color=tuple(np.array([0.3,0.3,0.3])+np.mod(d0,2)*0.3),linewidth=figParams['imshow_lineWidth']+1,alpha=1,linestyle = '-',solid_capstyle='butt')
                                ax[row,col].text(xlim_s[0],np.mean(dayRanges[behavior_epoch][d]),dayLabels[behavior_epoch][d],\
                                    color=tuple(np.array([0.3,0.3,0.3])+np.mod(d0,2)*0.3),ha='right',va='center',fontsize=figParams['scale_font']-5)

                if (figParams['formatMode'] == 'horz' and t1+1 == len(all_ttypes)) or (figParams['formatMode'] == 'vert' and t1 == 0):
                    if figParams['clean']:
                        clabel = "\n\n\n"+str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont))
                    else:
                        clabel = "\n\n\n"+alignLabel1+" "+str(smart_round(final_minCont))+"-"+str(smart_round(final_maxCont))

                    if figParams['match_trial_counts']:
                        ax[row,col].text(xlim_s[1],max_nTrials_display/2,clabel,fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')
                    else:
                        ax[row,col].text(xlim_s[1],allSessData.shape[0]/2,clabel,fontsize=figParams['scale_font']-2,rotation=90,ha='center',va='center')

                if figParams['formatMode'] == 'horz':
                    ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-4)
                    ax[row,col].set_xticks(xticks_s)
                    ax[row,col].set_xticklabels(xticks_s,fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'top','horizontalalignment': 'center'})
                else:
                    ax[row,col].set_xticks([])
                if t1 == 0:
                    if len(cluster):
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count)+"\n"+cluster+" "+str(cluster_idx),fontsize=figParams['scale_font']+2)
                    else:
                        ax[row,col].set_title(anm+" "+summaryInfo['ROI_types'][ROI_type_key][ROI_score]+str(summaryInfo['ROI_map'][anmIdx][ROI_idx]['ROI'])+" | ROI"+str(ROI_idx)+\
                            " | All"+str(page_count),fontsize=figParams['scale_font']+2)
                ax[row,col].spines['top'].set_visible(False)
                ax[row,col].spines['bottom'].set_visible(False)
                ax[row,col].spines['left'].set_visible(False)
                ax[row,col].spines['right'].set_visible(False)
                sess = 0
                if figParams['formatMode'] == 'horz':
                    tempTtype = ttype
                else:
                    tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                if figParams['match_trial_counts']:
                    ax[row,col].set_ylim([max_nTrials_display+2,-3])
                    ax[row,col].set_ylim([-3,max_nTrials_display+2])
                    ax[row,col].set_yticks(trial_ticks)
                    if t1 == 0:
                        ax[row,col].set_yticklabels([str(t) for t in trial_ticks],fontdict = {'fontsize': figParams['scale_font']-4,'verticalalignment': 'center','horizontalalignment': 'right'})
                    else:
                        if not figParams['yTicksOn']:
                            ax[row,col].set_yticklabels([])
                else:
                    ax[row,col].set_ylim([allSessData.shape[0]+2,-3])
                    ax[row,col].set_ylim([-3,allSessData.shape[0]+2])
                if not figParams['yTicksOn']:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                if figParams['yTicksOn'] and t1==0:
                    ax[row,col].set_ylabel("Trial Count",fontsize=figParams['scale_font']+2)

                figParams['im_scaleBar']['length'] = figParams['plot_scaleBar']['length_s']
                if not added_scaleBar[row,col]:
                    ax[row,col],figParams['im_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['im_scaleBar'],\
                                figParams['imshow_lineColor'],True,figParams['im_scaleBar']['vertUnit'],str(figParams['plot_scaleBar']['length_s'])+' '+figParams['im_scaleBar']['horzUnit'])

                    ax[row,col].invert_yaxis()
                    trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                    ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,[],\
                        True,True,False,False,figParams['imshow_lineColor'],figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                            nTrial_byTtype[ttype],-2,figParams['horzCueLineOn'],figParams['horzCueImScalar'],horzLoc =figParams['im_horzCueLocation'])
                    added_scaleBar[row,col] = True

        #################################################################
        #plots
        # row=1
        if figParams['filt_trace_byTrial']:
            plotStat = figParams['plot_stat']+"_filt"
        else:
            plotStat = figParams['plot_stat']
        maxVal=-1e6
        minVal=1e6
        max_nFr = 0
        max_s = 0
        for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
            if behavior_epoch in behavior_epochs:
                for t1,ttype in enumerate(all_ttypes):
                    max_nFr = np.nanmax([max_nFr,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']])
                    max_s = np.nanmax([max_s,aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_length_fr']/\
                                    aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['align_framerate']])
                    if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats']>=figParams['minRepeats_forYLim']:
                        if figParams['plot_error']:
                            if 'boot_CI' in figParams['plot_key_error_pos']:
                                maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']])])
                                minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']])])
                            else:
                                maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]+\
                                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']])])
                                minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]-\
                                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']])])
                        else:
                            maxVal=np.nanmax([maxVal,np.nanmax(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])])
                            minVal=np.nanmin([minVal,np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])])
        # max_nTrials_display = int(np.ceil(max_nTrials_display/figParams['match_trial_ticks'])*figParams['match_trial_ticks'])
        # trial_ticks = np.arange(0,max_nTrials_display+1,figParams['match_trial_ticks'])
        if minVal>=maxVal:
            minVal=0
            maxVal=1
        plotShift = np.arange((len(all_ttypes)))*(maxVal-minVal)*figParams['plotScalar']
        plotmaxVal=-1e6
        plotminVal=1e6
        for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
            if behavior_epoch in behavior_epochs:
                for t1,ttype in enumerate(all_ttypes):
                    if ttype in figParams['autoScale_be_ttypes'] and aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats']>=figParams['minRepeats_forYLim']:
                        if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats']>0:
                            xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                                load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)

                            
                            plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])
                            if figParams['plot_error']:
                                if 'boot_CI' in figParams['plot_key_error_pos']:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])
                                else:
                                    upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]+\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])
                                    lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]-\
                                                aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']])-\
                                                    np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-\
                                                np.nanmin(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])-plotShift[t1])

                            plotmaxVal=np.nanmax([plotmaxVal,np.nanmax(upperTrace[plotIdxs])])
                            plotminVal=np.nanmin([plotminVal,np.nanmin(lowerTrace[plotIdxs])])
        plotmaxVal = plotmaxVal+(plotmaxVal-plotminVal)*figParams['vertUpperPlotBuffer']
        plotminVal = plotminVal-(plotmaxVal-plotminVal)*figParams['vertLowerPlotBuffer']
        if plotminVal>=plotmaxVal:
            plotminVal=-1
            plotmaxVal=1
            plotShift = np.zeros((len(all_ttypes)))
        plotminVal = plotminVal-np.nanmean(plotShift)/2
        if figParams['manual_plotLims']:
            print(f'Current plotminVal: {plotminVal:.3f} plotmaxVal: {plotmaxVal:.3f} plotShift: {plotShift}')
            print("Applying Manual Plot Limits")
            final_plotmaxVal=figParams['manual_plotmaxVal'][ROI_type_key]
            final_plotminVal=figParams['manual_plotminVal'][ROI_type_key]
            final_plotShift=figParams['manual_plotShift'][ROI_type_key]
            print(f'Final plotminVal: {final_plotminVal:.3f} final_plotmaxVal: {final_plotmaxVal:.3f} final_plotShift: {final_plotShift}')
            if len(final_plotShift)!=len(all_ttypes):
                raise Exception("Length of manual_plotShift does not match number of ttypes for anmIdx "+str(final_plotShift))
        else:
            final_plotmaxVal=plotmaxVal
            final_plotminVal=plotminVal
            final_plotShift=plotShift
        vert_scale=smart_round(figParams['plot_scaleBar']['height_per']*(final_plotmaxVal-final_plotminVal))
        # vert_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.7]
        # horz_scale_coord=[xlim_s[1]-0.2,final_plotmaxVal*0.95]
        for t11,ttype in enumerate(reversed(all_ttypes)):
            sess = 0
            t1 = len(all_ttypes) - t11 - 1
            xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)                                            

            if figParams['formatMode'] == 'horz':
                row = 0
                for n in range(nPlots):
                    row+=1
                    col = t1
                    trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                    ax[row,col].set_xlim(xlim_s)
                    ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])

                    if not added_scaleBar[row,col]:
                        
                        ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                            True,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                                final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'],horzLoc = figParams['plot_horzCueLocation'])
                        ax[row,col],figParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['plot_scaleBar'],(0,0,0),True,\
                            figParams['vert_scale_label'])
                        added_scaleBar[row,col] = True


            else:
                row = nRows
                for n in range(nPlots):
                    row-=1
                    col = 0
                    if t11 == 0:
                        # print("adding trial features to row "+str(row)+" col "+str(col)+" for ttype "+ttype)
                        tempTtype = ''.join(copy.deepcopy(figParams['export_ttypes']))
                        # ax[row,col] = add_trial_structure_features_old('time',ax[row,col],figParams,summaryInfo,traceAlignParams,\
                        #     traceAlignParams['align_shift_s'],\
                        #     aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][sess_trace_type][sess][ttype]['align_framerate'],\
                        #     anmIdx,ROI_type_key,ROI_score,figParams['examples_group'],tempTtype,\
                        #     False,True,True,True,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],final_plotminVal,final_plotmaxVal)
                        trial_structure_times = load_trial_structure_times(summaryInfo,traceAlignParams,anmIdx,ROI_type_key,fix_overlaps=True,verbose=False)
                        ax[row,col].set_xlim(xlim_s)
                        ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
                        if not added_scaleBar[row,col]:
                            ax[row,col] = add_trial_structure_features('time',ax[row,col],figParams,trial_structure_times,ROI_type_key,ROI_score,figParams['examples_group'],align_data,ttype,all_ttypes,\
                                True,True,True,False,(0,0,0),figParams['imshow_lineWidth'],figParams['alpha'],figParams['imshow_lineWidth']+1,1,figParams['scale_font'],\
                                    final_plotminVal,final_plotmaxVal,figParams['horzCueLineOn'],figParams['horzCuePlotScalar'],horzLoc = figParams['plot_horzCueLocation'])
                            ax[row,col],figParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],figParams['plot_scaleBar'],(0,0,0),True,\
                                figParams['vert_scale_label'])
                            added_scaleBar[row,col] = True
            
            idx = -1
            zeroLines = []
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    idx+=1
                    if figParams['splitPrePost']:
                        if figParams['formatMode'] == 'horz':
                            row = 0 + (idx+1)
                        else:
                            row = nRows - (idx+1)
                    else:
                        row = nRows - 1
                    # print(behavior_epoch+" be = "+str(be)+" t1 = "+str(t1)+" row = "+str(row))
                    # if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['data'].shape[0]>0:
                    if aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats']>0:
                        xdata_s,xlim_s,xticks_s,xdata_fr,plotIdxs,framerate = \
                            load_alignment_xdata(traceAlignParams,ROI_type_key,figParams['examples_group'],zoom,figParams)
                                                        
                        plotMain = copy.deepcopy(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']])
                        if figParams['plot_error']:
                            if 'boot_CI' in figParams['plot_key_error_pos']:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']]))
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']]))
                            else:
                                upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]+\
                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_pos']]))
                                lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]-\
                                            aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_error_neg']]))
                        else:
                            upperTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]))
                            lowerTrace = copy.deepcopy((aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat][figParams['plot_key_main']]))
                        # tempColor = list(summaryInfo['ROI_types_colors'][ROI_type_key][ROI_score])
                        # if figParams['formatMode'] == 'horz':
                        if not -np.nanmin(plotMain[plotIdxs])-final_plotShift[t1] in zeroLines:
                            ax[row,col].axhline(-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1],color=(0,0,0),lw=figParams['imshow_lineWidth'],alpha=figParams['alpha'])
                            zeroLines.append(-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1])

                        if 'export_byPrePost_ttype_colors' in figParams:
                            tempColor = list(figParams['export_byPrePost_ttype_colors'][behavior_epoch][t1])
                        else:
                            tempColor = list(copy.deepcopy(summaryInfo['byPrePost_shift_colors'][be]))
                        for c in range(len(tempColor)):
                            if tempColor[c]>0.25:
                                tempColor[c] = tempColor[c]-0.25
                        # else:
                        #     tempColor = sess_colors[be]

                        if figParams['plot_error']:
                            ax[row,col].fill_between(xdata_s[plotIdxs],upperTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1],\
                                                                    lowerTrace[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1],\
                                                                        color=tempColor,alpha = figParams['alpha'],lw=0)
                        # tempColor = list(summaryInfo['ROI_types_colors'][ROI_type_key][ROI_score])
                        # if figParams['formatMode'] == 'horz':
                        if 'export_byPrePost_ttype_colors' in figParams:
                            tempColor = list(figParams['export_byPrePost_ttype_colors'][behavior_epoch][t1])
                        else:
                            tempColor = list(copy.deepcopy(summaryInfo['byPrePost_shift_colors'][be]))
                        if 'export_byPrePost_ttype_linestyle' in figParams:
                            tempLS = figParams['export_byPrePost_ttype_linestyle'][behavior_epoch][t1]
                        else:
                            tempLS = '-'
                        if 'export_byPrePost_ttype_linewidth' in figParams:
                            tempLW = figParams['export_byPrePost_ttype_linewidth'][behavior_epoch][t1]
                        else:
                            tempLW = figParams['lineWidth']
                        if 'export_byPrePost_ttype_alpha' in figParams:
                            tempAlpha = figParams['export_byPrePost_ttype_alpha'][behavior_epoch][t1]
                        else:
                            tempAlpha = 1
                        if 'export_byPrePost_ttype_dashes' in figParams:
                            tempDashes = figParams['export_byPrePost_ttype_dashes'][behavior_epoch][t1]
                        else:
                            tempDashes = None

                        # else:
                        #     tempColor = sess_colors[be]
                        
                        if figParams['clean']:
                            if not figParams['splitPrePost']:
                                if 'pre' in behavior_epoch.lower():
                                    bel = "Pre "
                                elif 'post' in behavior_epoch.lower():
                                    bel = "Post "
                                else:
                                    bel = behavior_epoch+" "
                            else:
                                bel = ""

                            label = bel+ttype+" (n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['orig_nRepeats'])+")"
                        else:
                            label = behavior_epoch+" "+ttype+"\n(n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['orig_nRepeats'])+" | "+\
                                            str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats'])+")"
                        if tempDashes is None:
                            ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1],\
                                            color=tempColor,ls = tempLS, linewidth=tempLW, alpha = tempAlpha, label=label)
                        else:
                            ax[row,col].plot(xdata_s[plotIdxs],plotMain[plotIdxs]-np.nanmin(plotMain[plotIdxs])-final_plotShift[t1],\
                                            color=tempColor,ls = tempLS, linewidth=tempLW, dashes = tempDashes, alpha = tempAlpha, label=label)
                        ax[row,col].set_xlim(xlim_s)
                        text_xPos = xlim_s[0]
                    if figParams['formatMode'] == 'horz':
                        text_yPos = -final_plotShift[t1]
                        if not np.isfinite(text_yPos):
                            print(text_yPos)
                        else:
                            ax[row,col].text(text_xPos,text_yPos,\
                                        aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype]['day_label']+\
                                        " n = "+str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['orig_nRepeats'])+" | "+\
                                            str(aligned_traces[anm][ROI_type_key][ROI_score][ROI_idx][trace_type][behavior_epoch][ttype][plotStat]['nRepeats']),\
                                            fontsize=figParams['scale_font']-4,color=(0,0,0),ha='left',va='top')
                    else:
                        # ax[row,col].legend(frameon = False, fontsize=figParams['scale_font']-4)
                        ax[row,col].legend(loc=figParams['legendLoc'],frameon=False,reverse=figParams['legendReversed'], bbox_to_anchor=figParams['legendBboxToAnchor'], ncol = 1,fontsize=figParams['legendFontsize'],\
                                                    handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3)  

            lineCount = idx
            idx = -1
            for be,behavior_epoch in enumerate(summaryInfo['byPrePost']):
                if behavior_epoch in behavior_epochs:
                    idx+=1
                    if figParams['splitPrePost']:
                        if figParams['formatMode'] == 'horz':
                            row = 0 + (idx+1)
                        else:
                            row = nRows - (idx+1)
                    else:
                        row = nRows - 1
                    
                    ax[row,col].set_ylim([final_plotminVal,final_plotmaxVal])
                    if idx == 0 and not figParams['clean']:
                        ax[row,col].text(xlim_s[0],final_plotmaxVal*0.95,behavior_epoch,fontsize=figParams['scale_font'],color=(0,0,0),ha='left',va='top')
                    ax[row,col].set_yticks([])
                    ax[row,col].set_yticklabels([])
                    
                    ax[row,col].spines['top'].set_visible(False)
                    ax[row,col].spines['right'].set_visible(False)
                    ax[row,col].spines['left'].set_visible(False)

                    if not figParams['clean']:
                        if traceAlignParams['align_reduce_factor'][ROI_type_key] == 0 or traceAlignParams['align_reduce_factor'][ROI_type_key] == 1:
                            ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s)",fontsize=figParams['scale_font']+2)
                        else:
                            ax[row,col].set_xlabel(traceAlignParams['align_position']+" Aligned Trial Time (s; bin"+str(traceAlignParams['align_reduce_factor'][ROI_type_key])+")",fontsize=figParams['scale_font']+2)
                        ax[row,col].tick_params(axis='both', which='major', labelsize=figParams['scale_font']-2)
                    else:
                        ax[row,col].set_xticks([])
                        ax[row,col].spines['bottom'].set_visible(False)


                   
        if figParams['formatMode'] == 'vert':
            col = 1
            for row in range(nRows):
                fig.delaxes(ax[row,col])
        if saveFigs:
            mergePDF.savefig(fig,bbox_inches='tight',pad_inches=0.05,dpi=600)  # saves the current figure into a mergePDF page
        if verbose:
            print("Finished!")
        if figParams['preview_figs'] and page_count <= nPreviewPages:
            display_clean_subplots(fig,ax)
        elif figParams['preview_figs'] and not page_count <= nPreviewPages:
            plt.close()
        else:
            if figParams['close_all_figs']:
                plt.close()
            else:
                facecolor = (1,1,1)
                fig.set_facecolor(facecolor)
                display_clean_subplots(fig,ax)

    
    return mergePDF
#########################################################################################################
##################################################################################################
def simple_trial_trace_summary_stats(input_data, additional_stats = 'std sem med cv var boot_med boot_mean boot_CI_pos boot_CI_neg', \
    bootstrap = True, bootstrapMethod = 'vectorized', bootStat = 'nanmean', nBoot = 1000, CI = 0.95, save_all_bstraps = False, \
    exclude_naned_trials = False, exclusion_nan_percentage = 100, save_input_data = False, verbose = False, \
    stabilityCheck = False, stabilityCheck_nRepeats = 10, stabilityCheck_ratio = 0.5, revert_NaNs = False, trial_ddof=0, boot_ddof=1):
    """
    Calculate summary statistics for a series of trial traces.
    by default only calculates the mean but you can include additional stats

    Parameters
    ----------
    input_data : np.ndarray
        Trial traces, shape (nRepeats, nFrames) (or 1D, treated as a single repeat). Rows are
        individual trial/repeat traces, columns are frames.
    additional_stats : str
        Space-separated keywords selecting which extra stats to compute, from
        'std sem sum med cv var boot_med boot_mean boot_CI_pos boot_CI_neg'.
    bootstrap : bool
        If True, also computes bootstrap-resampled stats requested via additional_stats
        (boot_mean/boot_med/boot_std/boot_sem/boot_CI_pos/boot_CI_neg).
    bootstrapMethod : str
        Only used when bootstrap=True. 'standard' | 'parallel' | 'vectorized' resampling
        implementation (see standard_bootstrap / run_parallel_bootstrap /
        vectorized_bootstrap).
    bootStat : str
        Only used when bootstrap=True. Statistic computed per bootstrap sample, e.g. 'nanmean'.
    nBoot : int
        Only used when bootstrap=True. Number of bootstrap resamples.
    CI : float
        Only used when bootstrap=True and boot_CI_pos/boot_CI_neg are requested. Confidence
        interval width, as a fraction (<1) or a percent (>=1).
    save_all_bstraps : bool
        Only used when bootstrap=True. If True, stores the full boot_data array (nBoot x nFrames)
        under summaryStats['boot_data'].
    exclude_naned_trials : bool
        If True, drops (rows sets to NaN then removes) trials whose fraction of NaN frames meets
        or exceeds exclusion_nan_percentage before computing stats.
    exclusion_nan_percentage : float
        Only used when exclude_naned_trials=True. NaN-frame percentage threshold (0-100) for
        dropping a trial.
    save_input_data : bool
        If True, stores a copy of input_data under summaryStats['data'].
    verbose : bool
        If True, prints progress/diagnostic info (e.g. number of trials excluded).
    stabilityCheck : bool
        If True, computes split-half correlation stability metrics (sequential split plus
        stabilityCheck_nRepeats random shuffles) under summaryStats['stability'].
    stabilityCheck_nRepeats : int
        Only used when stabilityCheck=True. Number of random split-half shuffles to run in
        addition to the sequential split.
    stabilityCheck_ratio : float
        Only used when stabilityCheck=True. Fraction (0-1) of trials assigned to the first half
        of each split.
    revert_NaNs : bool
        If True, replaces all NaNs in input_data with 0 before computing any stats.
    trial_ddof : int
        Delta degrees of freedom used for the across-trial std/sem/var (np.nanstd/np.nanvar).
    boot_ddof : int
        Only used when bootstrap=True and boot_std/boot_sem are requested. Delta degrees of
        freedom used for boot_std (np.nanstd over boot_data).

    Returns
    -------
    summaryStats : dict
        Dictionary of computed summary statistics (always includes 'mean', 'orig_nRepeats',
        'orig_nFrames', 'nRepeats', 'nFrames', 'delRepeats', 'stability'; other keys present
        depend on additional_stats/bootstrap/save_input_data/save_all_bstraps).

    Example
    -------
    # comments indented under a '> only if ...' marker = used only when that other arg is active
    additional_stats = 'std sem med cv var boot_med boot_mean boot_CI_pos boot_CI_neg'  # space-separated subset of: std sem sum med cv var boot_med boot_mean boot_std boot_sem boot_CI_pos boot_CI_neg
    bootstrap = True                    # True | False  (compute bootstrap-based stats)
    bootstrapMethod = 'vectorized'      #     > only if bootstrap=True: 'standard' | 'parallel' | 'vectorized'
    bootStat = 'nanmean'                #     > only if bootstrap=True: statistic per bootstrap sample, e.g. 'nanmean'
    nBoot = 1000                        #     > only if bootstrap=True: number of bootstrap resamples
    CI = 0.95                           #     > only if bootstrap=True and boot_CI_pos/neg requested: CI width as fraction (<1) or percent (>=1)
    save_all_bstraps = False            #     > only if bootstrap=True: True | False (save full boot_data array)
    exclude_naned_trials = False        # True | False  (drop trials with too many NaN frames)
    exclusion_nan_percentage = 100      #     > only if exclude_naned_trials=True: NaN-frame percentage threshold (0-100)
    save_input_data = False             # True | False  (store a copy of input_data in the output)
    verbose = False                     # True | False  (print progress/diagnostic info)
    stabilityCheck = False              # True | False  (compute split-half correlation stability)
    stabilityCheck_nRepeats = 10        #     > only if stabilityCheck=True: number of random split-half shuffles
    stabilityCheck_ratio = 0.5          #     > only if stabilityCheck=True: fraction (0-1) of trials in the first half
    revert_NaNs = False                 # True | False  (replace NaNs with 0 before computing stats)
    trial_ddof = 0                      # delta degrees of freedom for across-trial std/sem/var
    boot_ddof = 1                       #     > only if bootstrap=True and boot_std/boot_sem requested: delta degrees of freedom for boot_std
    importlib.reload(jsm)
    summaryStats = simple_trial_trace_summary_stats(
        input_data, additional_stats=additional_stats, bootstrap=bootstrap, bootstrapMethod=bootstrapMethod,
        bootStat=bootStat, nBoot=nBoot, CI=CI, save_all_bstraps=save_all_bstraps,
        exclude_naned_trials=exclude_naned_trials, exclusion_nan_percentage=exclusion_nan_percentage,
        save_input_data=save_input_data, verbose=verbose, stabilityCheck=stabilityCheck,
        stabilityCheck_nRepeats=stabilityCheck_nRepeats, stabilityCheck_ratio=stabilityCheck_ratio,
        revert_NaNs=revert_NaNs, trial_ddof=trial_ddof, boot_ddof=boot_ddof)
    """
    summaryStats={}
    if save_input_data:
        summaryStats['data'] = input_data.copy()  # ndarray .copy() instead of copy.deepcopy (direct memcpy, much faster)
    summaryStats['orig_nRepeats'] = 0
    summaryStats['orig_nFrames'] = 0
    summaryStats['delRepeats'] = list()
    summaryStats['nRepeats'] = 0
    summaryStats['nFrames'] = 0
    summaryStats['mean'] = np.array([])
    if 'std' in additional_stats:
        summaryStats['std'] = np.array([])
    if 'sem' in additional_stats:
        summaryStats['sem'] = np.array([])
    if 'sum' in additional_stats:
        summaryStats['sum'] = np.array([])  
    if 'median' in additional_stats:    
        summaryStats['median'] = np.array([])   
    if 'cv' in additional_stats:
        summaryStats['cv'] = np.array([])
    if 'var' in additional_stats:
        summaryStats['var'] = np.array([])
    if stabilityCheck:
        summaryStats['stability'] = {}
        summaryStats['stability']['R'] = list()
        summaryStats['stability']['p'] = list()
        summaryStats['stability']['R_mean'] = np.nan
        summaryStats['stability']['p_mean'] = np.nan
    if bootstrap:
        summaryStats['boot_mean'] = np.array([])
        summaryStats['boot_med'] = np.array([])
        summaryStats['boot_std'] = np.array([])
        summaryStats['boot_sem'] = np.array([])
        summaryStats['boot_CI_pos'] = np.array([])
        summaryStats['boot_CI_neg'] = np.array([])
        if save_all_bstraps:
            summaryStats['boot_data'] = np.array([])
    anyData = False
    if np.any(np.isfinite(input_data)):
        anyData = True
    else:
        if len(input_data.shape) == 1:
            if input_data.shape[0] > 0:
                anyData = True
        else:
            if input_data.shape[0] > 0 and input_data.shape[1] > 0:
                anyData = True
    if anyData:
        #input_data should be org repeats x frames to work properly
        output_data = input_data.copy()  # ndarray .copy() instead of copy.deepcopy (direct memcpy, much faster)
        if len(output_data.shape) == 1:
            output_data = np.expand_dims(output_data,axis=0)
        summaryStats['orig_nRepeats'],summaryStats['orig_nFrames'] = output_data.shape
        if summaryStats['orig_nFrames'] == 0 and summaryStats['orig_nRepeats'] == 1:
            summaryStats['orig_nFrames'] = 0
            summaryStats['orig_nRepeats'] = 0
        elif summaryStats['orig_nFrames'] == 1 and summaryStats['orig_nRepeats'] == 0:
            summaryStats['orig_nFrames'] = 0
            summaryStats['orig_nRepeats'] = 0
        summaryStats['orig_nRepeats'] = int(summaryStats['orig_nRepeats'])
        summaryStats['orig_nFrames'] = int(summaryStats['orig_nFrames'])

        #bootstrap params
        if bootstrap:
            summaryStats['nBoot'] = int(nBoot)
            if CI < 1:
                pcrt = CI * 100
            else:
                pcrt = CI
            summaryStats['boot_CI_pcrt'] = int(pcrt)

        # revert_NaNs removes all NaNs prior to calculations
        if revert_NaNs:
            output_data[np.isnan(output_data)]=0


        #delete any repeats/rows if the number NaNs is greater than or equal to the nan percentage of the total number of frames
        summaryStats['delRepeats'] = list()
        if exclude_naned_trials:
            for row in range(output_data.shape[0]):
                if np.sum(np.isnan(output_data[row,:])) >= int(output_data.shape[1]*(exclusion_nan_percentage/100)):
                    summaryStats['delRepeats'].append(row)
                    output_data[row,:] = np.nan
            if verbose:
                print("DELETING "+str(len(summaryStats['delRepeats']))+"/"+str(output_data.shape[0])+" Repeats")
            output_data=np.delete(output_data,summaryStats['delRepeats'],axis=0)
        summaryStats['delRepeats'] = np.array(summaryStats['delRepeats'])
        #new size
        summaryStats['nRepeats'],summaryStats['nFrames'] = output_data.shape
        if summaryStats['nFrames'] == 0 and summaryStats['nRepeats'] == 1:
            summaryStats['nRepeats'] = 0
        summaryStats['nRepeats'] = int(summaryStats['nRepeats'])
        summaryStats['nFrames'] = int(summaryStats['nFrames'])

        #Collect general summary stats
        summaryStats['mean'] = np.nanmean(output_data,axis=0).astype('float32')
        if 'std' in additional_stats:
            summaryStats['std'] = np.nanstd(output_data,axis=0,ddof=trial_ddof).astype('float32')
        if 'sem' in additional_stats:
            if not 'std' in additional_stats:
                summaryStats['sem'] = np.nanstd(output_data,axis=0,ddof=trial_ddof)/np.sqrt(summaryStats['nRepeats'])
            else:
                summaryStats['sem'] = summaryStats['std']/np.sqrt(summaryStats['nRepeats'])
            summaryStats['sem'] = summaryStats['sem'].astype('float32')
        if 'sum' in additional_stats:
            summaryStats['sum'] = np.nansum(output_data,axis=0).astype('float32')
        if 'med' in additional_stats:
            summaryStats['median'] = np.nanmedian(output_data,axis=0).astype('float32')
        if 'cv' in additional_stats:
            # Vectorized cv = std/mean; non-finite entries (e.g. mean==0) set to nan, matching
            # the prior per-frame loop's behavior but without the Python loop / per-element deepcopy.
            with np.errstate(divide='ignore', invalid='ignore'):
                cv_vals = summaryStats['std'] / summaryStats['mean']
            cv_vals = np.where(np.isfinite(cv_vals), cv_vals, np.nan)
            summaryStats['cv'] = cv_vals.reshape(-1,1).astype('float32')
        if 'var' in additional_stats:
            summaryStats['var'] = np.nanvar(output_data,axis=0).astype('float32')
        summaryStats['stability'] = {}
        summaryStats['stability']['R'] = list()
        summaryStats['stability']['p'] = list()
        summaryStats['stability']['R_mean'] = np.nan
        summaryStats['stability']['p_mean'] = np.nan
        if stabilityCheck:
            
            indices = list(range(summaryStats['nRepeats']))
            split_point = int(summaryStats['nRepeats'] * stabilityCheck_ratio)
            
            #First check the sequential order
            trace1 = np.nanmean(output_data[indices[:split_point],:],axis=0).astype('float32')
            trace2 = np.nanmean(output_data[indices[split_point:],:],axis=0).astype('float32')

            trace1 = trace1.flatten()
            trace2 = trace2.flatten()
            # pearsonr does not like nans so set them to zero
            if np.any(np.isfinite(trace1)) and np.any(np.isfinite(trace2)):
                trace1[np.isnan(trace1)] = 0
                trace2[np.isnan(trace2)] = 0
                trace1[np.isinf(trace1)] = 0
                trace2[np.isinf(trace2)] = 0
                correlation, p_value = pearsonr(trace1, trace2)
            else:
                correlation = np.nan
                p_value = np.nan
            summaryStats['stability']['R'].append(correlation)
            summaryStats['stability']['p'].append(p_value)

            #Then check the random order    
            for s in range(stabilityCheck_nRepeats):
                shuffled = indices.copy()
                random.shuffle(shuffled)
                trace1 = np.nanmean(output_data[shuffled[:split_point],:],axis=0).astype('float32')
                trace2 = np.nanmean(output_data[shuffled[split_point:],:],axis=0).astype('float32')

                trace1 = trace1.flatten()
                trace2 = trace2.flatten()
                # pearsonr does not like nans so set them to zero
                if np.any(np.isfinite(trace1)) and np.any(np.isfinite(trace2)):
                    trace1[np.isnan(trace1)] = 0
                    trace2[np.isnan(trace2)] = 0
                    trace1[np.isinf(trace1)] = 0
                    trace2[np.isinf(trace2)] = 0
                    correlation, p_value = pearsonr(trace1, trace2)
                else:
                    correlation = np.nan
                    p_value = np.nan
                summaryStats['stability']['R'].append(correlation)
                summaryStats['stability']['p'].append(p_value)
            summaryStats['stability']['R'] = np.array(summaryStats['stability']['R']).astype('float32')
            summaryStats['stability']['p'] = np.array(summaryStats['stability']['p']).astype('float32')
            summaryStats['stability']['R_mean'] = np.nanmean(summaryStats['stability']['R']).astype('float32')
            summaryStats['stability']['p_mean'] = np.nanmean(summaryStats['stability']['p']).astype('float32')
            # print(summaryStats['stability']['R_mean'])
        else:
            summaryStats['stability']['R_mean'] = np.nan
            summaryStats['stability']['p_mean'] = np.nan

        if bootstrap:
            if not 'boot_mean' in additional_stats and not 'boot_med' in additional_stats:
                print("additional_stats = "+str(additional_stats))
                raise Exception("ERROR: Bootstrapping is on but you arent collecting either the mean or the median of the bootstraps")

            if bootstrapMethod == 'standard':
                boot_data = standard_bootstrap(output_data, nBoot=summaryStats['nBoot'], axis = 0, replace=True, bootStat=bootStat)
            elif bootstrapMethod == 'parallel':
                boot_data = run_parallel_bootstrap(output_data, nBoot=summaryStats['nBoot'], axis = 0, replace=True, bootStat=bootStat, n_jobs=-1)
            elif bootstrapMethod == 'vectorized':
                boot_data = vectorized_bootstrap(output_data, nBoot=summaryStats['nBoot'], axis = 0, replace=True, bootStat=bootStat)
            #calculate stats for the boot_data
            if 'boot_mean' in additional_stats:
                summaryStats['boot_mean'] = np.nanmean(boot_data,axis=0).astype('float32')
            if 'boot_med' in additional_stats:
                summaryStats['boot_med'] = np.nanmedian(boot_data,axis=0).astype('float32')
            if 'boot_std' in additional_stats:
                summaryStats['boot_std'] = np.nanstd(boot_data,axis=0,ddof=boot_ddof).astype('float32')
            if 'boot_sem' in additional_stats:
                summaryStats['boot_sem'] = summaryStats['boot_std'].astype('float32')/np.sqrt(summaryStats['nBoot'])
            if 'boot_CI_pos' in additional_stats and 'boot_CI_neg' in additional_stats:
                # Compute both CI bounds in one call so boot_data is sorted once instead of twice.
                _ci_neg, _ci_pos = np.nanpercentile(boot_data,[(100-pcrt)/2,(100-pcrt)/2+pcrt],axis=0)
                summaryStats['boot_CI_pos'] = _ci_pos.astype('float32')
                summaryStats['boot_CI_neg'] = _ci_neg.astype('float32')
            elif 'boot_CI_pos' in additional_stats:
                summaryStats['boot_CI_pos'] = np.nanpercentile(boot_data,((100-pcrt)/2+pcrt),axis=0).astype('float32')
            elif 'boot_CI_neg' in additional_stats:
                summaryStats['boot_CI_neg'] = np.nanpercentile(boot_data,(100-pcrt)/2,axis=0).astype('float32')
            if save_all_bstraps:
                summaryStats['boot_data'] = boot_data
    return summaryStats

###################################
# bootstrapping functions
def bootstrap_summary_stats(output_data, bootstrap_stats, nBoot = 10000, 
    bootStat = 'nanmean', bootstrapMethod = 'standard', boot_ddof = 0, pcrt = 95, save_boot_data = False):
    """
    Calculate bootstrap summary statistics for the provided data.
    Parameters:
    output_data: input data to analyze
    bootstrap_stats: a dictionary to store the bootstrap statistics 
    nBoot: number of bootstrap iterations
    bootStat: statistic to calculate for each bootstrap sample (e.g., 'nanmean', 'nanmedian')
    bootstrapMethod: method to use for bootstrapping ('standard', 'parallel', 'vectorized')
    boot_ddof: degrees of freedom to use when calculating standard deviation (default is 0 for
    population standard deviation)
    pcrt: percentile for confidence interval (default is 95 for 95% confidence interval)
    save_boot_data: if True, save the bootstrap samples in the bootstrap_stats dictionary
    Returns:
    bootstrap_stats: a dictionary containing the bootstrap statistics
    """
    bootstrap_stats['nBoot'] = nBoot
    bootstrap_stats['bootStat'] = bootStat  
    if bootstrapMethod == 'standard':
        boot_data = standard_bootstrap(output_data, nBoot=bootstrap_stats['nBoot'], axis = 0, replace=True, bootStat=bootStat)
    elif bootstrapMethod == 'parallel':
        boot_data = run_parallel_bootstrap(output_data, nBoot=bootstrap_stats['nBoot'], axis = 0, replace=True, bootStat=bootStat, n_jobs=-1)
    elif bootstrapMethod == 'vectorized':
        boot_data = vectorized_bootstrap(output_data, nBoot=bootstrap_stats['nBoot'], axis = 0, replace=True, bootStat=bootStat)
    bootstrap_stats['boot_mean'] = np.nanmean(boot_data,axis=0).astype('float32')
    bootstrap_stats['boot_med'] = np.nanmedian(boot_data,axis=0).astype('float32')
    bootstrap_stats['boot_std'] = np.nanstd(boot_data,axis=0,ddof=boot_ddof).astype('float32')
    bootstrap_stats['boot_sem'] = bootstrap_stats['boot_std'].astype('float32')/np.sqrt(bootstrap_stats['nBoot'])
    bootstrap_stats['boot_CI_pos'] = np.nanpercentile(boot_data,((100-pcrt)/2+pcrt),axis=0).astype('float32')
    bootstrap_stats['boot_CI_neg'] = np.nanpercentile(boot_data,(100-pcrt)/2,axis=0).astype('float32')
    if save_boot_data:
        bootstrap_stats['boot_data'] = boot_data.astype('float32')
    return bootstrap_stats

def standard_bootstrap(data, nBoot = 1000, axis = 0, replace = True, bootStat = 'nanmean'):
    """Performs standard bootstrapping by resampling the data and calculating the specified statistic for each bootstrap sample
    Parameters:
    data: input data to analyze
    nBoot: number of bootstrap iterations
    axis: axis along which to perform the resampling (default is 0 for resampling rows)
    replace: if True, resampling is done with replacement (default is True)     
    bootStat: statistic to calculate for each bootstrap sample (e.g., 'nanmean', 'nanmedian')
    Returns:
    boot_data: a 2D array containing the bootstrap statistics for each iteration
    """
    n_samples = data.shape[axis]
    if axis == 0:
        nFrames = data.shape[axis+1]
    else:
        nFrames = data.shape[axis-1]
    boot_data = np.zeros((nBoot,nFrames)).astype('float32')
    for boot in range(nBoot):
        if bootStat == 'nanmean':
            boot_data[boot,:] = np.nanmean(data[resample(range(n_samples), replace=replace),:],axis=axis)
    return boot_data

def bootstrap_task(data, axis, replace, bootStat):
    """
    Performs a single bootstrap iteration: 
    Resamples rows (axis 0) and calculates the mean.
    """
    n_samples = data.shape[axis]
    # Generate indices for resampling
    indices = np.random.choice(n_samples, size=n_samples, axis=axis, replace=replace)
    # Resample and calculate mean along axis
    if bootStat == 'mean':
        return np.mean(data[indices, :], axis=axis)
    elif bootStat == 'median':
        return np.mean(data[indices, :], axis=axis)
    elif bootStat == 'nanmean':
        return np.mean(data[indices, :], axis=axis)
    elif bootStat == 'nanmedian':
        return np.mean(data[indices, :], axis=axis)
    else:
        raise Exception("PROVIDE VALID BOOTSTRAPPING STAT")

def run_parallel_bootstrap(data, nBoot=1000, axis = 0, replace=True, bootStat='nanmean', n_jobs=-1):
    """
    Distributes bootstrap iterations across multiple CPU cores.
    n_jobs = -1 uses all avail
    """
    # Parallel execution
    results = Parallel(n_jobs=n_jobs)(
        delayed(bootstrap_task)(data, axis, replace, bootStat) for _ in range(nBoot)
    )
    
    return np.array(results)

def vectorized_bootstrap(data, nBoot=1000, axis = 0, replace=True, bootStat='nanmean'):
    """
    Highly optimized bootstrapping using NumPy broadcasting.
    """
    n_samples = data.shape[axis]

    # 1. Generate a 2D matrix of random indices all at once
    # Shape: (nBoot, n_samples)
    # randint is faster than choice for the uniform-with-replacement case (no probability/permutation setup).
    if replace:
        indices = np.random.randint(0, n_samples, size=(nBoot, n_samples))
    else:
        indices = np.random.choice(n_samples, size=(nBoot, n_samples), replace=replace)

    # 2. Calculate the statistic along the 'n_samples' axis.
    # NaN-aware reductions allocate a full bool mask over the (nBoot, n_samples, n_cols) array and
    # do extra passes; when the data contains no NaNs the plain reductions give identical results
    # much faster. Check NaNs once on the small 2D input, not the large 3D resample.
    nan_free = not np.isnan(data).any()

    # Fast path for the bootstrap-of-mean on NaN-free 2D data: a resampled mean is a weighted sum of
    # the original rows, where the weight on row i is the number of times it was drawn. So
    #   boot_mean[b] = (1/n) * sum_i counts[b,i] * data[i]   ==   (counts @ data) / n.
    # This avoids ever materializing the (nBoot, n_samples, n_cols) resample -- the heavy work is a
    # single BLAS matmul, and n_cols never appears in an intermediate. Matches the gather+mean result
    # up to floating-point summation order (negligible at these sizes). Mean only; median can't be
    # expressed as a weighted row-sum, and NaNs would break the per-row weighting.
    use_matmul_mean = bootStat in ('mean', 'nanmean') and nan_free and replace and axis == 0 and data.ndim == 2
    if use_matmul_mean:
        # counts[b,i] = number of times row i was drawn in bootstrap b, built fully vectorized via a
        # single offset bincount (one bin per (b,i) cell).
        flat = indices + n_samples * np.arange(nBoot)[:, None]
        counts = np.bincount(flat.ravel(), minlength=nBoot * n_samples).reshape(nBoot, n_samples).astype(data.dtype)
        boot_data = (counts @ data) / n_samples
    else:
        # Materialize the 3D resample only when actually needed (median, or NaN-aware reductions).
        # data[indices] creates a 3D array of shape (nBoot, n_samples, n_cols).
        resampled_data = data[indices]
        if bootStat == 'mean' or (bootStat == 'nanmean' and nan_free):
            boot_data = resampled_data.mean(axis=axis+1)
        elif bootStat == 'median' or (bootStat == 'nanmedian' and nan_free):
            boot_data = np.median(resampled_data, axis=axis+1)
        elif bootStat == 'nanmean':
            boot_data = np.nanmean(resampled_data, axis=axis+1)
        elif bootStat == 'nanmedian':
            boot_data = np.nanmedian(resampled_data, axis=axis+1)
        else:
            raise Exception("PROVIDE VALID BOOTSTRAPPING STAT")
    return boot_data

###################################
#Image and trace example plots
def image_and_trace_movie_preview(pdf,data,anat,sess,traceKeys,imageData,imageParams,cropParams,movieParams, \
                                  importFrs,ROIs,ROI_borders,ROI_colors,plotParams,trial_details,f,ROI_title,eventFrame,saveFig,clean=False):

    """
    Render a single frame of the image + trace movie into a preview PDF page.

    A non-export render path used to build preview PDFs. It assembles a one-frame
    ``imageArray`` (with optional per-channel filtering), recollects traces and
    trial-epoch labels (``collect_plot_traces``, ``trial_epoch_frameLabels``,
    ``collect_eventFrames``), builds the figure with ``image_and_trace_fig`` /
    ``image_and_trace_panels``, and optionally saves the figure as a page onto the
    open ``PdfPages`` object.

    Example
    -------
    saveFig = True
    clean = False
    importlib.reload(jsm)
    pdf = image_and_trace_movie_preview(
        pdf, data, anat, sess, traceKeys, imageData, imageParams, cropParams, movieParams,
        importFrs, ROIs, ROI_borders, ROI_colors, plotParams, trial_details, f, ROI_title, eventFrame, saveFig, clean=clean
    )

    Parameters
    ----------
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Open multi-page PDF to append the rendered frame to.
    data : dict
        Nested behavioral/trace data.
    anat : str
        Anatomy/ROI-type key into ``data``.
    sess : int or str
        Session key.
    traceKeys : dict
        Map from trace names to ``data`` keys.
    imageData : dict
        Per-channel image data dict.
    imageParams : dict
        Per-channel display parameters.
    cropParams : dict
        Image crop parameters.
    movieParams : dict
        Movie export parameters.
    importFrs : list of int
        Frame indices defining the trace x-axis.
    ROIs : list of int
        ROI indices to display.
    ROI_borders : list
        ROI contours.
    ROI_colors : dict
        Map from ROI id to RGB color.
    plotParams : dict
        Trace-plotting parameters.
    trial_details : dict
        Per-trial metadata.
    f : int
        Local index into ``importFrs`` of the frame to render.
    ROI_title : str
        Title text.
    eventFrame : bool
        Event flag (recomputed internally).
    saveFig : bool
        If True, write the figure as a page to ``pdf``.
    clean : bool
        If True, suppress titles.

    Returns
    -------
    pdf : matplotlib.backends.backend_pdf.PdfPages
        The (possibly appended) PDF object.

    Local Dependencies
    -------------------
    Other functions defined in this file (py) that this
    function calls, directly or transitively — duplicate all of these into a
    pared-down repo alongside image_and_trace_movie_preview() itself:
      - collect_plot_traces
      - trial_epoch_frameLabels
      - collect_eventFrames
      - image_and_trace_fig
      - image_and_trace_panels
          - update_image_panels
              - frames_to_seconds
              - get_frameLabels
              - frameLabel_text_kwargs
    """
    preview = False
    frame = importFrs[f]
    fidx = f
    nFrs = len(importFrs)
    if not 'exportChs' in movieParams:
        movieParams['exportChs'] = list(imageParams.keys())
    nCh = len(movieParams['exportChs'])
    imageArrayShape = (1,imageData[movieParams['exportChs'][0]]['data'].shape[1],imageData[movieParams['exportChs'][0]]['data'].shape[2],nCh)
    imageArray = np.ones(imageArrayShape,dtype='float32')*np.nan
    for c,ch in enumerate(movieParams['exportChs']):
        if 'data' in imageData[ch]:
            if np.any(imageData[ch]['data']):
                tempData = copy.deepcopy(imageData[ch]['data'][f,:,:])
                if imageParams[ch]['filter']:
                    tempData=gaussian_filter(tempData, imageParams[ch]['filterSigma_px'])
                imageArray[0,:,:,c] = copy.deepcopy(tempData)
    plotData,plotParams['plotLims'],eventTraces,plotParams['bouts'],plotParams['events'],plotParams['licks'] =  \
        collect_plot_traces(data,traceKeys,anat,sess,ROIs,plotParams,importFrs)
    
    movieParams = trial_epoch_frameLabels(movieParams,trial_details)

    eventFrame = collect_eventFrames(movieParams,plotParams,fidx)

    fig,ax,movieParams,xdata = image_and_trace_fig(movieParams,imageParams,cropParams,importFrs,plotParams,ROI_borders,ROI_title)
    fig,ax,frame_trackers = image_and_trace_panels(fig,ax,xdata,imageArray,plotData,imageParams,cropParams,ROI_borders,movieParams,ROI_colors,\
                                                   plotParams,trial_details,frame,fidx,0,ROI_title,eventFrame,preview,clean)
    if saveFig:
        pdf.savefig(fig, bbox_inches='tight', dpi=600)
    if movieParams['exportClosePreviewPDF']:
        plt.close()
    else:
        plt.show()
    return pdf

def collect_plot_traces(data,traceKeys,anat,sess,ROIs,plotParams,importFrs):
    """
    Gather all per-ROI trace data, events, licks, and bouts for the exported frames.

    Slices the session's trace arrays at ``importFrs`` for every ROI and every trace
    listed in ``plotParams['traceGroups']``, producing the ``plotData`` dict keyed by
    (ROI column, row, trace index). Also derives per-ROI event frame indices (from the
    'eventKey' trace), NaN-masked bout markers (from the 'boutKey' trace), and left/
    right lick frames, and computes per-panel y-limits padded by the configured
    top/bottom buffers. Consumed by ``image_and_trace_panels``.

    Example
    -------
    importlib.reload(jsm)
    plotData, plotLims, eventTraces, bouts, events, licks = collect_plot_traces(
        data, traceKeys, anat, sess, ROIs, plotParams, importFrs
    )

    Parameters
    ----------
    data : dict
        Nested behavioral/trace data indexed by anatomy -> sessions -> session -> key.
    traceKeys : dict
        Map from trace names to keys inside the ``data`` structure.
    anat : str
        Anatomy/ROI-type key into ``data``.
    sess : int or str
        Session key into ``data[anat]['sessions']``.
    ROIs : list of int
        ROI indices to extract (define the trace columns).
    plotParams : dict
        Trace-plotting parameters ('traceGroups', 'eventKey', 'boutKey',
        'topBuffer', 'bottomBuffer').
    importFrs : list of int
        Frame indices to slice from each trace.

    Returns
    -------
    plotData : dict
        Trace values keyed by (col, row, trace_idx).
    plotLims : dict
        Per-(col, row) [ymin, ymax] limits, buffered.
    eventTraces : dict
        Per-ROI raw event trace slices.
    bouts : dict
        Per-ROI NaN-masked bout markers.
    events : dict
        Per-ROI list of frame indices where an event occurs.
    licks : dict
        'right'/'left' lists of frame indices with licks.
    """
    licks = {}
    licks['right'] = []
    licks['left'] = []
    if np.any(data[anat]['sessions'][sess]['licks_per_fr'][0][importFrs]>0):
        licks['right'] = np.argwhere(data[anat]['sessions'][sess]['licks_per_fr'][0][importFrs]>0).flatten().tolist()
    if np.any(data[anat]['sessions'][sess]['licks_per_fr'][1][importFrs]>0):
        licks['left'] = np.argwhere(data[anat]['sessions'][sess]['licks_per_fr'][1][importFrs]>0).flatten().tolist()
    
    eventTraces = {}
    bouts = {}
    events = {}
    for col, ROI in enumerate(ROIs):
        eventTraces[ROI] = copy.deepcopy(data[anat]['sessions'][sess][traceKeys[plotParams['eventKey']]][ROI,importFrs])
        bouts[ROI] = copy.deepcopy(data[anat]['sessions'][sess][traceKeys[plotParams['boutKey']]][ROI,importFrs]).astype('float32')
        bouts[ROI][bouts[ROI]==0] = np.nan
        events[ROI] = []
        for f in range(len(eventTraces[ROI])):
            if eventTraces[ROI][f]>0:
                events[ROI].append(f)
    plotData = {}
    plotLims = {}
    for col, ROI in enumerate(ROIs):
        for row, traces in enumerate(plotParams['traceGroups']):
            plotLims[col,row] = [1e6,-1e6]
    for col, ROI in enumerate(ROIs):
        for row, traces in enumerate(plotParams['traceGroups']):
            for t,trace in enumerate(traces):
                plotData[col,row,t] = copy.deepcopy(data[anat]['sessions'][sess][traceKeys[trace]][ROI,importFrs])
                if 'events' in traceKeys[trace]:
                    plotData[col,row,t][plotData[col,row,t]==0] = np.nan
                plotLims[col,row][0] = np.nanmin([plotLims[col,row][0],np.nanmin(plotData[col,row,t])])
                plotLims[col,row][1] = np.nanmax([plotLims[col,row][1],np.nanmax(plotData[col,row,t])])
    for col, ROI in enumerate(ROIs):
        for row, traces in enumerate(plotParams['traceGroups']):
            plotLims[col,row][0] = plotLims[col,row][0] - np.absolute(plotLims[col,row][1]-plotLims[col,row][0])*plotParams['bottomBuffer'][row]
            plotLims[col,row][1] = plotLims[col,row][1] + np.absolute(plotLims[col,row][1]-plotLims[col,row][0])*plotParams['topBuffer'][row]

    return plotData,plotLims,eventTraces,bouts,events,licks

def image_and_trace_fig(movieParams,imageParams,cropParams,importFrs,plotParams,ROI_borders,ROI_title):
    """
    Create the matplotlib grid for an image + trace movie frame.

    Determines the grid size and creates the axes via
    ``clean_subplots``. Rows span the larger of the channel
    count and the number of trace groups; columns hold the image panel(s) plus one
    column per ROI when traces are on. ``movieParams['splitImgCol']`` splits the
    image channels across two columns instead of stacking them. Unused axes are
    deleted, and the image crop offsets ('im_xadjust'/'im_yadjust') are computed.

    Example
    -------
    ROI_title = movieParams['saveName']
    importlib.reload(jsm)
    fig, ax, movieParams, xdata = image_and_trace_fig(
        movieParams, imageParams, cropParams, importFrs, plotParams, ROI_borders, ROI_title
    )

    Parameters
    ----------
    movieParams : dict
        Movie export parameters. Reads 'splitImgCol', 'exportChs', 'tracesOn',
        'cropData'; writes 'nRows', 'nCols', 'im_xadjust', 'im_yadjust'.
    imageParams : dict
        Per-channel display parameters (used for channel iteration/parity).
    cropParams : dict
        Image crop bounds with 'xlim'/'ylim'.
    importFrs : list of int
        Frame indices being exported (defines the trace x-axis length).
    plotParams : dict
        Trace-plotting parameters; reads 'traceGroups', 'colScalar', 'rowScalar'.
    ROI_borders : list
        ROI contours; its length sets the number of trace columns.
    ROI_title : str
        Title text (kept for signature parity).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    ax : np.ndarray
        2-D array of axes handles.
    movieParams : dict
        Updated movieParams (see above).
    xdata : np.ndarray
        Integer x-axis (frame indices) for the trace panels.
    """
    if movieParams['splitImgCol']:
        movieParams['nRows'] = int(np.max([len(movieParams['exportChs'])/2,len(plotParams['traceGroups'])]))
        if movieParams['tracesOn']:
            movieParams['nCols'] = int(len(ROI_borders)+1+1)
        else:
            movieParams['nCols'] = int(len(ROI_borders)+1)
    else:
        movieParams['nRows'] = int(np.max([len(movieParams['exportChs']),len(plotParams['traceGroups'])]))
        if movieParams['tracesOn']:
            movieParams['nCols'] = int(len(ROI_borders)+1)
        else:
            movieParams['nCols'] = int(len(ROI_borders))
    xdata = np.arange(len(importFrs))
    if movieParams['cropData']:
        movieParams['im_xadjust'] = cropParams['xlim'][0]
        movieParams['im_yadjust'] = cropParams['ylim'][0]

    else:
        movieParams['im_xadjust'] = 0
        movieParams['im_yadjust'] = 0


    fig,ax = clean_subplots(movieParams['nRows'],movieParams['nCols'],\
        figsize=(movieParams['nCols']*plotParams['colScalar'],movieParams['nRows']*plotParams['rowScalar']))
    ax = np.atleast_1d(ax)
    if ax.ndim == 1:
        if movieParams['nRows']>movieParams['nCols']:
            ax = ax[:, np.newaxis]
        else:
            ax = ax[np.newaxis, :]

    plt.subplots_adjust(wspace=0, hspace=0)   
    if movieParams['splitImgCol']:
        if np.mod(len(movieParams['exportChs']),2) == 1:
            row = int(len(movieParams['exportChs'])/2)-1
            col = 1
            # print(f'row {row} col {col}')
            fig.delaxes(ax[row,col])
        if np.ceil(len(movieParams['exportChs'])/2)<movieParams['nRows']:
            nrow2del = movieParams['nRows']-int(np.ceil(len(movieParams['exportChs'])/2))
            for rr in range(nrow2del):
                row = int(len(movieParams['exportChs'])/2)+rr
                col = 0
                fig.delaxes(ax[row,col])
                col = 1
                fig.delaxes(ax[row,col])
    else:
        if len(movieParams['exportChs'])<movieParams['nRows']:
            col = 0
            row = movieParams['nRows']
            for i in range(movieParams['nRows']-len(movieParams['exportChs'])):
                row-=1
                fig.delaxes(ax[row,col])
    # print(frame_trackers)
    # print(fig)
    # print(ax)
    # plt.suptitle(ROI_title,y=0.9,fontsize=20)       
    return fig,ax,movieParams,xdata

def trial_epoch_frameLabels(movieParams,trial_details):
    """
    Precompute trial-epoch frame labels (cue and GO windows) for overlay.

    Builds ``movieParams['frameLabels']`` as a list of (start_fr, end_fr, label,
    color) tuples: a left/right "Cue" window starting at each trial's sample start
    and a green "GO Cue" window at the go start, with durations derived from the
    trial timing and ``movieParams['imaging_FPS']``. Any user-supplied
    ``movieParams['manualFrameLabels']`` are appended. Consumed later by
    ``get_frameLabels`` during rendering.

    Example
    -------
    importlib.reload(jsm)
    movieParams = trial_epoch_frameLabels(movieParams, trial_details)

    Parameters
    ----------
    movieParams : dict
        Movie export parameters; reads 'imaging_FPS' and 'manualFrameLabels',
        writes 'frameLabels'.
    trial_details : dict
        Per-trial metadata (epoch frame indices, durations, left/right cue flags).

    Returns
    -------
    movieParams : dict
        Updated movieParams with the populated 'frameLabels' list.
    """
    movieParams['frameLabels'] = []
    for trial in reversed(trial_details.keys()):
        if trial_details[trial]['right']:
            temp = (trial_details[trial]['sample_start_fr'],trial_details[trial]['sample_start_fr']+int(np.ceil(trial_details[trial]['sample_s']*movieParams['imaging_FPS'])),"Right Cue (0.1 s)",(1,1,1))
        else:
           temp = (trial_details[trial]['sample_start_fr'],trial_details[trial]['sample_start_fr']+int(np.ceil(trial_details[trial]['sample_s']*movieParams['imaging_FPS'])),"Left Cue (0.1 s)",(1,1,1))
        if not temp in movieParams['manualFrameLabels']:
            movieParams['frameLabels'].append(temp)
        temp = (trial_details[trial]['go_start_fr'],trial_details[trial]['go_start_fr']+int(np.ceil(trial_details[trial]['go_actual_s']*movieParams['imaging_FPS'])),"GO Cue (2.6 s)",(0,1,0))
        if not temp in movieParams['frameLabels']:
            movieParams['frameLabels'].append(temp)
    for temp in movieParams['manualFrameLabels']:
        if not temp in movieParams['frameLabels']:
            movieParams['frameLabels'].append(temp)
    # print(movieParams['frameLabels'])
    return movieParams

def collect_eventFrames(movieParams,plotParams,fidx):
    """
    Test whether a frame is at or near an event, for the event marker.

    Returns True if any ROI has an event at ``fidx`` or within
    ``movieParams['eventPersistenceFr']`` frames before/after it (using the per-ROI
    event indices in ``plotParams['events']``). Used by the generator/exporter to set
    the ``eventFrame`` flag that adds a '*' marker to the channel label.

    Example
    -------
    importlib.reload(jsm)
    eventFrame = collect_eventFrames(movieParams, plotParams, fidx)

    Parameters
    ----------
    movieParams : dict
        Movie export parameters; reads 'eventPersistenceFr'.
    plotParams : dict
        Trace-plotting parameters; reads 'events' (per-ROI event frame indices).
    fidx : int
        Trial-relative frame index to test.

    Returns
    -------
    eventFrame : bool
        True if an event occurs at/near ``fidx``.
    """
    eventFrame = False
    for r, ROI in enumerate(plotParams['events'].keys()):
        if fidx in plotParams['events'][ROI]:
            eventFrame = True
        for ff in range(movieParams['eventPersistenceFr']+1):
            if fidx+ff in plotParams['events'][ROI]:
                eventFrame = True
        for ff in range(movieParams['eventPersistenceFr']+1):
            if fidx-ff in plotParams['events'][ROI]:
                eventFrame = True
    return eventFrame

def image_and_trace_panels(fig,ax,xdata,imageArray,plotData,imageParams,cropParams,\
    ROI_borders,movieParams,ROI_colors,plotParams,trial_details,frame,fidx,f,ROI_title,eventFrame,preview=False,clean=False):
    """
    Populate the full figure for one frame: image panels plus static trace panels.

    Called once to build the initial frame. It realizes each channel's colormap
    (``generate_cmap``), draws the image panels via ``update_image_panels``, then
    repositions the image axes and, when ``movieParams['tracesOn']``, builds each
    ROI's trace column: trial-epoch background patches (sample/delay/response
    Rectangles), cue/response vertical lines, the trace curves themselves, event /
    lick / bout overlays, and the current-frame marker line whose handles are
    returned in ``frame_trackers`` for cheap per-frame updates by
    ``update_panels_only``.

    Example
    -------
    eventFrame = False
    preview = False
    clean = False
    importlib.reload(jsm)
    fig, ax, frame_trackers = image_and_trace_panels(
        fig, ax, xdata, imageArray, plotData, imageParams, cropParams, ROI_borders, movieParams,
        ROI_colors, plotParams, trial_details, frame, fidx, f, ROI_title, eventFrame, preview=preview, clean=clean
    )

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure created by ``image_and_trace_fig``.
    ax : np.ndarray
        2-D array of axes handles.
    xdata : np.ndarray
        Integer x-axis (frame indices) for the trace panels.
    imageArray : np.ndarray
        Image stack (nFrames, H, W, nChannels).
    plotData : dict
        Per-(ROI, trace group, trace) data from ``collect_plot_traces``.
    imageParams : dict
        Per-channel display parameters; 'cmap1' is set here via ``generate_cmap``.
    cropParams : dict
        Image crop bounds (offsets computed here if absent).
    ROI_borders : list
        ROI contours; its length sets the number of trace columns.
    movieParams : dict
        Movie export parameters (reads 'tracesOn', 'splitImgCol', layout/spacing).
    ROI_colors : dict
        Map from ROI id to RGB color.
    plotParams : dict
        Trace-plotting parameters (groups, colors, patches, scale bars, events).
    trial_details : dict
        Per-trial metadata (epoch frame indices, cue info).
    frame : int
        Original (absolute) frame index.
    fidx : int
        Trial-relative frame index.
    f : int
        Frame index into imageArray's first axis.
    ROI_title : str
        Title text.
    eventFrame : bool
        Whether this frame is at/near an event.
    preview : bool
        If True, skip the current-frame marker (static preview render).
    clean : bool
        If True, suppress titles.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The populated figure.
    ax : np.ndarray
        The updated axes array.
    frame_trackers : list
        Handles of the current-frame marker lines (empty in preview mode).
    """
    if not cropParams:
        cropParams['xlim']=[0,imageArray.shape[2]]
        cropParams['ylim']=[0,imageArray.shape[1]]
    cropParams['xoffset'] = np.absolute(cropParams['xlim'][1]-cropParams['xlim'][0])*0.02
    cropParams['yoffset'] = np.absolute(cropParams['ylim'][1]-cropParams['ylim'][0])*0.02
    for d0,d in enumerate(imageParams):
        imageParams[d]['cmap1'], _, _ = generate_cmap(imageParams[d]['cmap'],movieParams['colorScalar'],cmap_name=d)
        imageParams[d]['cmap1'].set_bad(color=movieParams['nan_color'])

    frame_trackers = list()
    col = 0
    ax = update_image_panels(ax,imageArray,imageParams,cropParams,plotParams,movieParams,ROI_borders,ROI_colors,frame,fidx,f,ROI_title,eventFrame,preview,clean)
    imEdge = 0
    count = -1
    for d0,d in enumerate(imageParams):
        if movieParams['splitImgCol']:
            count+=1
            if d0 < np.ceil(len(movieParams['exportChs'])/2):
                col = 0
            else:
                col = 1
            row = copy.deepcopy(count)
            if count+1 == np.ceil(len(movieParams['exportChs'])/2):
                count = -1
            # pos1 = ax[row,col].get_position() # get the original position
            if movieParams['tracesOn']:
                if col == 0:
                    pos0 = [0,1-(row+1)*(1/np.ceil(len(movieParams['exportChs'])/2)),\
                        1/(len(ROI_borders)+2)-(movieParams['imageSpacer']),(1/np.ceil(len(movieParams['exportChs'])/2)-(movieParams['imageSpacer']))]
                else:
                    pos0 = [0.25,1-(row+1)*(1/np.ceil(len(movieParams['exportChs'])/2)),\
                        1/(len(ROI_borders)+2)-(movieParams['imageSpacer']),(1/np.ceil(len(movieParams['exportChs'])/2)-(movieParams['imageSpacer']))]
            else:
                if col == 0:
                    pos0 = [0,1-(row+1)*(1/np.ceil(len(movieParams['exportChs'])/2)),\
                        0.5-(movieParams['imageSpacer']),(1/np.ceil(len(movieParams['exportChs'])/2)-(movieParams['imageSpacer']))]
                else:
                    pos0 = [0.5,1-(row+1)*(1/np.ceil(len(movieParams['exportChs'])/2)),\
                        0.5-(movieParams['imageSpacer']),(1/np.ceil(len(movieParams['exportChs'])/2)-(movieParams['imageSpacer']))]
            # print(f'{d} {d0} {row} {col} {pos0}')
        else:
            row = d0
            col = 0
            # pos1 = ax[row,col].get_position() # get the original position 
            pos0 = [0,1-(row+1)*(1/len(movieParams['exportChs'])),1/(len(ROI_borders)+1),(1/len(movieParams['exportChs']))]
        # print(f'{d} {d0} {row} {col} {pos0}')
        imEdge = np.max([imEdge,pos0[0]+pos0[2]])
        if col<movieParams['nCols'] and row<movieParams['nRows']:
            ax[row,col].set_position(pos0) # set a new position
    imEdge = imEdge+0.02
    # print(f'imEdge {imEdge}')
    # print(imEdge)
    if movieParams['tracesOn']:
        traceWidth = ((1-imEdge)/len(ROI_borders))*plotParams['widthAdjust'] 
        for r, ROI in enumerate(plotParams['events'].keys()):
            if movieParams['splitImgCol']:
                col = r+2
            else:
                col = r+1
            for row1, traces in enumerate(reversed(plotParams['traceGroups'])):
                row = len(plotParams['traceGroups'])-1-row1
                # for row, traces in enumerate(plotParams['traceGroups']):
                # print(f'row {row} col {col} traceGroup {traces}')
                # ax[row,col].cla()
                ax[row,col].set_ylim([plotParams['plotLims'][r,row][0],plotParams['plotLims'][r,row][1]])
                for t in trial_details.keys():
                    # ts = ax[row,col].axvline(trial_details[t]['sample_start_fr'],color=(0,0,0),lw=1,alpha=0.5)
                    SampleRectangle=Rectangle(((trial_details[t]['sample_start_fr']),plotParams['plotLims'][r,row][0]+np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0),\
                        trial_details[t]['sample_len_fr'],np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])-np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0)
                    DelayRectangle=Rectangle(((trial_details[t]['delay_start_fr']),plotParams['plotLims'][r,row][0]+np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0),\
                        trial_details[t]['delay_len_fr'],np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])-np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0)
                    ResponseRectangle=Rectangle(((trial_details[t]['go_start_fr']),plotParams['plotLims'][r,row][0]+np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0),\
                        trial_details[t]['go_len_fr'],np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])-np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0)
                    patches = []
                    colors = []
                    if plotParams['samplePatchOn']:
                        patches.append(SampleRectangle)
                        if trial_details[t]['left']:
                            colors.append(plotParams['facecolors_left'][0])
                        elif trial_details[t]['right']:
                            colors.append(plotParams['facecolors_right'][0])
                        else:
                            colors.append(plotParams['facecolors_generic'][0])
                    if plotParams['delayPatchOn']:
                        patches.append(DelayRectangle)
                        if trial_details[t]['left']:
                            colors.append(plotParams['facecolors_left'][1])
                        elif trial_details[t]['right']:
                            colors.append(plotParams['facecolors_right'][1])
                        else:
                            colors.append(plotParams['facecolors_generic'][1])
                    if plotParams['responsePatchOn']:
                        patches.append(ResponseRectangle)
                        if trial_details[t]['left']:
                            colors.append(plotParams['facecolors_left'][2])
                        elif trial_details[t]['right']:
                            colors.append(plotParams['facecolors_right'][2])
                        else:
                            colors.append(plotParams['facecolors_generic'][2])

                    GeneralPatches = PatchCollection(patches, facecolor=colors, alpha=plotParams['patch_alpha'], edgecolor=plotParams['edgecolor'])
                    ax[row,col].add_collection(GeneralPatches);
                    
                    if plotParams['vertSampleLineOn']:
                        ax[row,col].axvline(trial_details[t]['sample_start_fr'],color=(0,0,0),lw=plotParams['vert_lw'],alpha=0.5)
                    if plotParams['vertReponseLineOn']:
                        ax[row,col].axvline(trial_details[t]['go_start_fr'],color=(0,0,0),lw=plotParams['vert_lw'],alpha=0.5)
                    

                    # if plotParams['horzCueLineOn']:
                    SampleTimes = trial_details[t]['sample_len_fr']/movieParams['imaging_FPS']
                    DelayTimes = trial_details[t]['delay_len_fr']/movieParams['imaging_FPS']
                    ResponseTimes = trial_details[t]['go_len_fr']/movieParams['imaging_FPS']
                    StartTime = trial_details[t]['trial_start_fr']/movieParams['imaging_FPS']
                    SampleTime = trial_details[t]['sample_start_fr']/movieParams['imaging_FPS']
                    DelayTime = trial_details[t]['delay_start_fr']/movieParams['imaging_FPS']
                    ResponseTime = trial_details[t]['go_start_fr']/movieParams['imaging_FPS']
                    SampleFrames = trial_details[t]['sample_len_fr']
                    DelayFrames = trial_details[t]['delay_len_fr']
                    ResponseFrames = trial_details[t]['go_len_fr']
                    StartFrame = trial_details[t]['trial_start_fr']
                    SampleFrame = trial_details[t]['sample_start_fr']
                    DelayFrame = trial_details[t]['delay_start_fr']
                    ResponseFrame = trial_details[t]['go_start_fr']


                    mixedCue = False
                    if trial_details[t]['left']:
                        cue = 'left'
                    elif trial_details[t]['right']:
                        cue = 'right'
                    else:
                        cue = 'left'
                        mixedCue = True
                    sampleTrace_ydata,sampleTrace_xdata,sampleTrace_env = generate_pulsed_tone(
                        carrier_hz=plotParams['cueInfo'][cue]['sim']['carrier_hz'] ,       # set your tone frequency here
                        pulse_rate_hz=plotParams['cueInfo'][cue]['sim']['pulse_rate_hz'] ,       # 4 pulses per second
                        on_duration_s=plotParams['cueInfo'][cue]['sim']['on_duration_s'] ,   # ON duration per pulse
                        duration_s=float(trial_details[t]['sample_len_fr'])/movieParams['imaging_FPS'] ,          # total length
                        sr=plotParams['cueInfo'][cue]['sim']['sr'] ,
                        amplitude=plotParams['cueInfo'][cue]['sim']['amplitude'] ,
                        fade_ms=plotParams['cueInfo'][cue]['sim']['fade_ms'] ,             # tweak for softer/harder edges
                        return_envelope=True)
                    cue = 'go'
                    responseTrace_ydata,responseTrace_xdata,responseTrace_env = generate_pulsed_tone(
                        carrier_hz=plotParams['cueInfo'][cue]['sim']['carrier_hz'] ,       # set your tone frequency here
                        pulse_rate_hz=plotParams['cueInfo'][cue]['sim']['pulse_rate_hz'] ,       # 4 pulses per second
                        on_duration_s=plotParams['cueInfo'][cue]['sim']['on_duration_s'] ,      # ON duration per pulse
                        duration_s=float(trial_details[t]['go_len_fr'])/movieParams['imaging_FPS']  ,          # total length
                        sr=plotParams['cueInfo'][cue]['sim']['sr'] ,
                        amplitude=plotParams['cueInfo'][cue]['sim']['amplitude'] ,
                        fade_ms=plotParams['cueInfo'][cue]['sim']['fade_ms'] ,             # tweak for softer/harder edges
                        return_envelope=True)
                    responseTrace_ydata[responseTrace_env<=0] = np.nan

                    sampleTrace_xdata_frame = sampleTrace_xdata*movieParams['imaging_FPS']
                    responseTrace_xdata_frame = responseTrace_xdata*movieParams['imaging_FPS']
                    cueHeight = plotParams['horzCueScalar'] *np.abs(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])
                    if plotParams['horzLoc'] == 'top':
                        horzCueLineMax = np.nanmax([np.nanmax(sampleTrace_ydata*cueHeight+plotParams['plotLims'][r,row][1]),np.nanmax(responseTrace_ydata*cueHeight+plotParams['plotLims'][r,row][1])])
                        # plotParams['plotLims'][r,row][1] = plotParams['plotLims'][r,row][1] - (horzCueLineMax - plotParams['plotLims'][r,row][1])
                        horzPlotLoc = plotParams['plotLims'][r,row][1]-cueHeight
                        va = 'top'
                    else:
                        horzCueLineMin = np.nanmin([np.nanmin(sampleTrace_ydata*cueHeight+plotParams['plotLims'][r,row][0]),np.nanmin(responseTrace_ydata*cueHeight+plotParams['plotLims'][r,row][0])])
                        # plotParams['plotLims'][r,row][0] = plotParams['plotLims'][r,row][0] + (plotParams['plotLims'][r,row][0] - horzCueLineMin)
                        horzPlotLoc = plotParams['plotLims'][r,row][0]+cueHeight
                        va = 'bottom'
                    # print(f'row {row} col {col} trial {t} plotParams["horzLoc"] {plotParams["horzLoc"]} horzPlotLoc {horzPlotLoc} StartFrame {StartFrame} SampleFrame {SampleFrame} DelayFrame {DelayFrame} ResponseFrame {ResponseFrame}')
                    # print(f'sampleTrace_xdata {sampleTrace_xdata[0]} - {sampleTrace_xdata[-1]}')
                    # print(f'sampleTrace_xdata_frame {sampleTrace_xdata_frame[0]} - {sampleTrace_xdata_frame[-1]}')
                    # print(f'{sampleTrace_xdata_frame[0]+(StartFrame)+(SampleFrame-StartFrame)}')
                    # print(f'{responseTrace_xdata_frame[0]+(StartFrame)+(ResponseFrame-StartFrame)}')
                    # print(plotParams['horzLineAlpha'])
                    if trial_details[t]['right']:
                        if plotParams['startMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(sampleTrace_xdata_frame+(StartFrame)+(SampleFrame-StartFrame),sampleTrace_ydata*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_right'][0],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([SampleFrame,SampleFrame+SampleFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_right'][0],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['delayMark']:
                            ax[row,col].plot([DelayFrame,DelayFrame+DelayFrames],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=plotParams['facecolors_right'][1],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['goMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(responseTrace_xdata_frame+(StartFrame)+(ResponseFrame-StartFrame),responseTrace_ydata*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_right'][2],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([ResponseFrame,ResponseFrame+ResponseFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_right'][2],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                    elif trial_details[t]['left']:
                        if plotParams['startMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(sampleTrace_xdata_frame+(StartFrame)+(SampleFrame-StartFrame),sampleTrace_ydata*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_left'][0],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([SampleFrame,SampleFrame+SampleFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_left'][0],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['delayMark']:
                            ax[row,col].plot([DelayFrame-StartFrame,DelayFrame-StartFrame+DelayFrames],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=plotParams['facecolors_left'][1],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['goMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(responseTrace_xdata_frame+(StartFrame)+(ResponseFrame-StartFrame),responseTrace_ydata*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_left'][2],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([ResponseFrame,ResponseFrame+ResponseFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_left'][2],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                    else:
                        if plotParams['startMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(np.insert(sampleTrace_xdata_frame, 0, -1*sampleTrace_xdata[1])+(StartFrame)+(SampleFrame-StartFrame),np.insert(sampleTrace_env, 0, 0)*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_generic'][0],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([SampleFrame,SampleFrame+SampleFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_generic'][0],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['delayMark']:
                            ax[row,col].plot([DelayFrame-StartFrame,DelayFrame-StartFrame+DelayFrames],[horzPlotLoc,horzPlotLoc],\
                                linestyle='-',color=plotParams['facecolors_generic'][1],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                        if plotParams['goMark']:
                            if plotParams['horzCueLineOn']:
                                ax[row,col].plot(responseTrace_xdata_frame+(StartFrame)+(ResponseFrame-StartFrame),responseTrace_ydata*cueHeight+horzPlotLoc+plotParams['horzCueShift'],\
                                    linestyle='-',color=plotParams['facecolors_generic'][2],linewidth=plotParams['horzCueLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')
                            else:
                                ax[row,col].plot([ResponseFrame,ResponseFrame+ResponseFrames],[horzPlotLoc,horzPlotLoc],\
                                    linestyle='-',color=plotParams['facecolors_generic'][2],linewidth=plotParams['horzLineWidth'],alpha=plotParams['horzLineAlpha'],solid_capstyle='butt')

                    if row == 0:
                        # -np.absolute(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0.02
                        ax[row,col].text(trial_details[t]['sample_start_fr'],plotParams['plotLims'][r,row][1],\
                                        trial_details[t]['short_label'],ha='left',va='bottom',color=trial_details[t]['textColor'],fontsize=plotParams['fontSize'])
                if not preview:
                    if plotParams['currFrame']:
                        ft = ax[row,col].axvline(fidx,ymin = plotParams['currFrame_ymin'],ymax = plotParams['currFrame_ymax'],\
                            color=plotParams['currFrame_color'],lw=plotParams['currFrame_lw'],alpha=plotParams['currFrame_alpha'])
                        frame_trackers.append(ft)
                if plotParams['showBouts'][row]:
                    traceBout = copy.deepcopy(plotParams['bouts'][ROI])
                    if plotParams['boutPosition'] == 'top':
                        for i in range(len(traceBout)):
                            traceBout[i] = plotParams['plotLims'][r,row][1] -(traceBout[i]+1) * (plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['boutHeightSplit']
                    elif plotParams['boutPosition'] == 'bottom':
                        for i in range(len(traceBout)):
                            traceBout[i] = plotParams['plotLims'][r,row][0] + (traceBout[i]+1) * (plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['boutHeightSplit']
                    
                    ax[row,col].plot(xdata,traceBout,\
                        plotParams['boutLineStyle'],color=plotParams['boutColor'],lw=plotParams['boutLineWidth'],alpha = plotParams['boutAlpha'], label = plotParams['boutLabel'])

                for t,trace in enumerate(traces):
                    ax[row,col].plot(xdata,plotData[r,row,t],plotParams['traceGroupStyles'][row][t],\
                        color=plotParams['traceGroupColors'][row][t],lw=plotParams['traceGroupLineWidths'][row][t],\
                            ms=plotParams['traceGroupMarkerSizes'][row][t],label = plotParams['traceLabels'][row][t])
                if plotParams['showEventHashes'][row]:
                    # print("row = "+str(row))
                    # print("plotParams['plotLims'][r,row][0] = "+str(plotParams['plotLims'][r,row][0]))
                    # print("plotParams['plotLims'][r,row][1] = "+str(plotParams['plotLims'][r,row][1]))
                    # print("plotParams['eventHeightScalar'] = "+str(plotParams['eventHeightScalar']))
                    count = 0
                    for e in plotParams['events'][ROI]:
                        # ax[row,col].axvline(e,ymin=plotParams['plotLims'][r,row][0],ymax=plotParams['plotLims'][r,row][0]+(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['eventHeightScalar'],color=plotParams['event_color'],lw=1,alpha=0.5)
                        count += 1
                        if count == 1:
                            ax[row,col].plot([e,e],[plotParams['plotLims'][r,row][0],\
                                plotParams['plotLims'][r,row][0]+(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['eventHeightScalar']],\
                                color=plotParams['eventColor'],lw=plotParams['eventLineWidth'],alpha=plotParams['eventAlpha'], label = plotParams['eventLabel'])
                        else:
                            ax[row,col].plot([e,e],[plotParams['plotLims'][r,row][0],\
                                plotParams['plotLims'][r,row][0]+(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['eventHeightScalar']],\
                                color=plotParams['eventColor'],lw=plotParams['eventLineWidth'],alpha=plotParams['eventAlpha'])
                if plotParams['showLicks'][row]:
                    count = 0
                    if plotParams['lickPos'] == 'top':
                        lickYdata = [plotParams['plotLims'][r,row][1],\
                                    plotParams['plotLims'][r,row][1]-\
                                        (plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['lickHeightScalar']]
                    else:
                        lickYdata = [plotParams['plotLims'][r,row][0],\
                                    plotParams['plotLims'][r,row][0]+\
                                        (plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*plotParams['lickHeightScalar']]
                    for e in plotParams['licks']['right']:
                        count += 1
                        if count == 1:
                            ax[row,col].plot([e,e],lickYdata,\
                                color=plotParams['rightLickColor'],lw=plotParams['lickLineWidth'],alpha=plotParams['lickAlpha'], label = 'Right Licks')
                        else:
                            ax[row,col].plot([e,e],lickYdata,\
                                color=plotParams['rightLickColor'],lw=plotParams['lickLineWidth'],alpha=plotParams['lickAlpha'])
                    count = 0
                    for e in plotParams['licks']['left']:
                        count += 1
                        if count == 1:
                            ax[row,col].plot([e,e],lickYdata,\
                                color=plotParams['leftLickColor'],lw=plotParams['lickLineWidth'],alpha=plotParams['lickAlpha'], label = 'Left Licks')
                        else:
                            ax[row,col].plot([e,e],lickYdata,\
                                color=plotParams['leftLickColor'],lw=plotParams['lickLineWidth'],alpha=plotParams['lickAlpha'])
                if clean:
                    ax[row,col].set_yticks([])
                    ax[row,col].set_xticks([])
                    ax[row,col].spines['top'].set_visible(False)
                    ax[row,col].spines['right'].set_visible(False)
                    ax[row,col].spines['left'].set_visible(False)
                    ax[row,col].spines['bottom'].set_visible(False)
                else:
                    if row+1<len(plotParams['traceGroups']):
                        ax[row,col].set_xticks([])
                    else:
                        ax[row,col].set_xlabel('Frames',fontsize=plotParams['fontSize'])
                    ax[row,col].set_ylabel(plotParams['traceGroupYLabels'][row],fontsize=plotParams['fontSize'])
                    ax[row,col].set_xlim([xdata[0],xdata[-1]])
                    ax[row,col].set_ylim([plotParams['plotLims'][r,row][0],plotParams['plotLims'][r,row][1]])
                    ax[row,col].spines['top'].set_visible(False)
                    ax[row,col].spines['right'].set_visible(False)
                    if row == len(plotParams['traceGroups'])-1:
                        ax[row,col].spines['bottom'].set_visible(True)
                    else:
                        ax[row,col].spines['bottom'].set_visible(False)
                        ax[row,col].axhline(plotParams['plotLims'][r,row][0]+(plotParams['plotLims'][r,row][1]-plotParams['plotLims'][r,row][0])*0.005,color=(0,0,0),lw=1)
                if plotParams['plot_scaleBar']['length_s'] <= 1 and (xdata[-1]-xdata[0])/ movieParams['imaging_FPS'] > 60:
                    plotParams['plot_scaleBar']['length_s'] = 10
                plotParams['plot_scaleBar']['length_fr'] = np.round(float(plotParams['plot_scaleBar']['length_s'] * movieParams['imaging_FPS']))
                ax[row,col],plotParams['plot_scaleBar'] = add_plot_scaleBar(ax[row,col],plotParams['plot_scaleBar'],(0,0,0),True,plotParams['traceGroupUnits'][row])
                ax[row,col].legend(frameon=False,fontsize=plotParams['legendFontSize'],loc=plotParams['legendLoc'],bbox_to_anchor=plotParams['legendBoxToAnchor'],\
                    handlelength=0.75, markerscale=0.75, handletextpad=0.2, labelspacing=0.1,  borderpad=0.3) 
                # pos1 = ax[row,col].get_position() # get the original position 
                pos0 = [imEdge,1-(row+1)*(1/len(plotParams['traceGroups'])),traceWidth,(1/len(plotParams['traceGroups']))*1]
                ax[row,col].set_position(pos0) # set a new position
    return fig,ax,frame_trackers

def update_image_panels(ax,imageArray,imageParams,cropParams,plotParams,movieParams,ROI_borders,ROI_colors,frame,fidx,f,ROI_title,eventFrame,preview,clean=False):
    """
    Redraw the image panels for a single frame (called every frame).

    For each exported channel: clears the axis, colorizes the frame with
    ``convert2RGB`` (using the channel's contrast and colormap),
    optionally crops, then draws the image plus its overlays — scale bar
    (``add_image_scaleBar``), title, channel label, timestamp/frame labels
    (``get_frameLabels``), contrast readout, ROI borders (``draw_ROI_borders``),
    manual markers (``annotate_image``), and an optional border rectangle —
    finishing with ``imshow_cleanup``. Panel placement follows
    ``movieParams['splitImgCol']``.

    Example
    -------
    eventFrame = False
    preview = False
    clean = False
    importlib.reload(jsm)
    ax = update_image_panels(
        ax, imageArray, imageParams, cropParams, plotParams, movieParams, ROI_borders, ROI_colors,
        frame, fidx, f, ROI_title, eventFrame, preview, clean=clean
    )

    Parameters
    ----------
    ax : np.ndarray
        2-D array of axes handles.
    imageArray : np.ndarray
        Image stack (nFrames, H, W, nChannels).
    imageParams : dict
        Per-channel display parameters ('cont', 'cmap1', 'filter', 'filterSigma_px',
        'label', 'ROI_color').
    cropParams : dict
        Image crop parameters with 'xlim', 'ylim', 'xoffset', 'yoffset'.
    plotParams : dict
        Trace-plotting parameters (unused here, kept for call-signature parity).
    movieParams : dict
        Movie export parameters (display toggles, scale bar, panel layout, markers).
    ROI_borders : list
        ROI contours drawn when 'ROI_borders' display is on.
    ROI_colors : dict
        Map from ROI id to RGB color.
    frame : int
        Original (absolute) frame index for labels/timestamps.
    fidx : int
        Trial-relative frame index for labels/timestamps.
    f : int
        Frame index into imageArray's first axis.
    ROI_title : str
        Title drawn on panel (0,0) when not clean.
    eventFrame : bool
        Whether this frame is at/near an event (adds a marker to the label).
    preview : bool
        Preview flag (kept for parity).
    clean : bool
        If True, suppress the title.

    Returns
    -------
    ax : np.ndarray
        The updated axes array.
    """
    count = -1
    # for d0,d in enumerate(imageParams):
    for d0,d in enumerate(movieParams['exportChs']):
        temp_ROI_colors = {}
        if imageParams[d]['ROI_color']:
            for ROI in ROI_colors:
                temp_ROI_colors[ROI] = imageParams[d]['ROI_color']
        if len(movieParams['exportChs'])>1:
            textColor = (1,1,1)
        else:
            textColor = imageParams[d]['ROI_color']

        if movieParams['splitImgCol']:
            count+=1
            if d0 < np.ceil(len(movieParams['exportChs'])/2):
                col = 0
            else:
                col = 1
            row = copy.deepcopy(count)
            if count+1 == np.ceil(len(movieParams['exportChs'])/2):
                count = -1
        else:
            row = d0
            col = 0
        ax[row,col].cla()
        imgRGB = convert2RGB(imageArray[f,:,:,d0],\
            imageParams[d]['cont'][0],imageParams[d]['cont'][1],\
            imageParams[d]['cmap1'],movieParams['colorScalar'],movieParams['nan_color'])
        imgRGB = imgRGB[:,:,0:3]
        if movieParams['cropData']:
            imgRGB = imgRGB[cropParams['ylim'][0]:cropParams['ylim'][1]+1,cropParams['xlim'][0]:cropParams['xlim'][1]+1,:]

        # im = ax[row,col].imshow(imageArray[f,:,:,d0],vmin=imageParams[d]['cont'][0],vmax=imageParams[d]['cont'][1],interpolation = 'none',cmap = imageParams[d]['cmap'])
        im = ax[row,col].imshow(imgRGB,interpolation = 'none',clip_on=False)
        if movieParams['cropData']:
            ax[row,col].set_xlim([0,imgRGB.shape[1]])
            ax[row,col].set_ylim([imgRGB.shape[0],0])
        else:
            ax[row,col].set_xlim([cropParams['xlim'][0]-movieParams['im_xadjust'],cropParams['xlim'][1]-movieParams['im_xadjust']])
            ax[row,col].set_ylim([cropParams['ylim'][1]-movieParams['im_yadjust'],cropParams['ylim'][0]-movieParams['im_yadjust']])
        if textColor == (0,0,0):
            scaleColor = (0,0,0)
        else:
            scaleColor = (1,1,1)
        if row==0 and col==0:
            includeLabel = True
        else:
            includeLabel = False
        movieParams['image_scaleBar']['fontsize'] = movieParams['fontSize']
        movieParams['image_scaleBar']['length_px'] = float(movieParams['image_scaleBar']['length_um'] / movieParams['pixel_sizes_umpx'])
        ax[row,col],movieParams['image_scaleBar'] = add_image_scaleBar(ax[row,col], movieParams['image_scaleBar'], color = scaleColor,includeLabel = includeLabel)
        if row == 0 and col == 0 and not clean:
            # ax[row,col].set_title(ROI_title+" | Frames: Trial "+str(fidx)+" Orig "+str(frame),color=(0,0,0),\
            #                       fontsize=plotParams['fontSize'],ha='left',va='top')
        
            ax[row,col].set_title(ROI_title+"\nFrames: Trial "+str(fidx)+" Orig "+str(frame),color=(0,0,0),\
                                  fontsize=movieParams['fontSize'],ha='left',va='bottom')
        
        if imageParams[d]['filter']:
            label = imageParams[d]['label']+"(σ"+str(imageParams[d]['filterSigma_px'])+"px Filt.)"        
        else:
            label = imageParams[d]['label']
        if movieParams['displayEventFrames'] and eventFrame:
            label = label+"\n*"
        contLabel = str(imageParams[d]['cont'][1])+"\n"+str(imageParams[d]['cont'][0])
        if movieParams['displayLabel']:
            ax[row,col].text(cropParams['xlim'][0]+cropParams['xoffset']-movieParams['im_xadjust'],cropParams['ylim'][0]+cropParams['yoffset']-movieParams['im_yadjust'],\
                            label,color=textColor,fontsize=movieParams['fontSize'],ha='left',va='top')
        
        all_frameLabels = []
        if movieParams['displayTimestamp'] and col == movieParams['timestamp_col'] and row == movieParams['timestamp_row']:
            ts_unit = movieParams.get('timestamp_unit', 's')
            ts_zero = movieParams.get('timestamp_zeroFrame', 0)
            ts_dec = movieParams.get('timestamp_decimals', 2)
            if movieParams['displayOverallTime']:
                timestamp = frames_to_seconds(frame, movieParams['imaging_FPS'],ts_dec,ts_unit,ts_zero)
            else:
                timestamp = frames_to_seconds(fidx, movieParams['imaging_FPS'],ts_dec,ts_unit,ts_zero)
            all_frameLabels.append((timestamp,(1,1,1),False))
        else:
            timestamp = ""
        if movieParams['displayFrameLabels'] and col == movieParams['frameLabel_col'] and row == movieParams['frameLabel_row']:

            frameLabels = get_frameLabels(fidx, movieParams)
            for fl in frameLabels:
                all_frameLabels.append((fl[0],fl[1],True))
        # print(all_frameLabels)
        if movieParams['displayFrameLabels'] or movieParams['displayTimestamp']:
            x = cropParams['xlim'][1]-cropParams['xoffset']-movieParams['im_xadjust']
            y = cropParams['ylim'][0]+cropParams['yoffset']-movieParams['im_yadjust']
            fl_kwargs = frameLabel_text_kwargs(movieParams)
            fl_lineHeight = movieParams.get('frameLabelFormat',{}).get('lineHeight')
            if fl_lineHeight is None:
                fl_lineHeight = movieParams['labelLineHeight']
            for text, color, is_fl in all_frameLabels:
                if is_fl:
                    ax[row,col].text(x,y,str(text),color=color,**fl_kwargs)
                    y += fl_lineHeight
                else:
                    ax[row,col].text(x,y,str(text),\
                                color=color,fontsize=movieParams['fontSize'],ha='right',va='top')
                    y += movieParams['labelLineHeight']


        if movieParams['displayContrast']:
            ax[row,col].text(cropParams['xlim'][1]-cropParams['xoffset']-movieParams['im_xadjust'],cropParams['ylim'][1]-cropParams['yoffset']-movieParams['im_yadjust'],\
                            contLabel,color=textColor,fontsize=movieParams['fontSize'],ha='right',va='bottom')
        if movieParams['ROI_borders']:
            ax[row,col] = draw_ROI_borders(ax[row,col],ROI_borders,temp_ROI_colors,ls='-',lw=movieParams['ROI_lw'],\
                alpha=movieParams['ROI_alpha'],xadjust=movieParams['im_xadjust'],yadjust=movieParams['im_yadjust'])
        # fig.colorbar(im, ax=ax[row,col])
            
        for a in movieParams['manualImageMarkers']:
            annotate_image(ax[row,col], a["x"], a["y"],
                        mode=a["mode"],
                        u= a["u"], v= a["v"],
                        fmt=a["fmt"], label=a["label"])
            
        
        ax[row,col].set_axis_off();
        if movieParams['imshowBorder']:
            xlim = ax[row, col].get_xlim()
            ylim = ax[row, col].get_ylim()
            x0 = min(xlim)
            y0 = min(ylim)
            width = max(xlim) - min(xlim)
            height = max(ylim) - min(ylim)
            rect = patches.Rectangle((x0, y0), width, height,
                                    linewidth=movieParams['imshowBorder_lw'],
                                    edgecolor=movieParams['imshowBorder_color'], facecolor='none', alpha=movieParams['imshowBorder_alpha'],
                                    transform=ax[row, col].transData,
                                    clip_on=False)
            ax[row, col].add_patch(rect)
        # ax[row,col].set_xlim([-0.5,sess_data['im_w']+0.5]);
        # ax[row,col].set_ylim([sess_data['im_h']+0.5,-0.5]);
        ax[row,col]=imshow_cleanup(ax[row,col])
    return ax

def frameLabel_text_kwargs(movieParams):
    """
    Build matplotlib ``Text`` kwargs for rendering trial-epoch frame labels.

    Reads the formatting overrides in ``movieParams['frameLabelFormat']`` and
    resolves the fall-throughs: a None 'fontsize' uses ``movieParams['fontSize']``
    and any None-valued entry (e.g. 'fontfamily', 'rotation') is dropped so
    matplotlib applies its own default. The returned dict is splatted into
    ``ax.text(...)`` by ``update_image_panels`` / ``render_image_panels``. The
    separate 'lineHeight' override (consumed by the callers, not a Text kwarg) is
    not included here.

    Example
    -------
    importlib.reload(jsm)
    kwargs = frameLabel_text_kwargs(movieParams)

    Parameters
    ----------
    movieParams : dict
        Movie export parameters; reads 'frameLabelFormat' (a dict of Text-kwarg
        overrides) and 'fontSize' (fallback font size).

    Returns
    -------
    kwargs : dict
        Keyword arguments ready to splat into ``ax.text(...)`` (e.g. 'fontsize',
        'fontweight', 'fontstyle', 'alpha', 'ha', 'va', and optionally
        'fontfamily'/'rotation').
    """
    fmt = movieParams.get('frameLabelFormat', {}) or {}
    kwargs = {
        'fontsize': fmt.get('fontsize') if fmt.get('fontsize') is not None else movieParams['fontSize'],
        'fontweight': fmt.get('fontweight', 'normal'),
        'fontstyle': fmt.get('fontstyle', 'normal'),
        'alpha': fmt.get('alpha', 1.0),
        'ha': fmt.get('ha', 'right'),
        'va': fmt.get('va', 'top'),
    }
    if fmt.get('fontfamily') is not None:
        kwargs['fontfamily'] = fmt['fontfamily']
    if fmt.get('rotation') is not None:
        kwargs['rotation'] = fmt['rotation']
    return kwargs

def get_frameLabels(frame, movieParams):
    """
    Return the (text, color) frame labels active at a given frame.

    Looks up ``movieParams['frameLabels']`` (a list of (start, end, label, color)
    tuples built by ``trial_epoch_frameLabels``) and returns those active at
    ``frame``. When ``movieParams['persistentFrameLabels']`` is True a label stays
    active for all frames at/after its start; otherwise it is active only within
    [start, end]. Used by ``update_image_panels`` to overlay trial-epoch labels.

    Example
    -------
    importlib.reload(jsm)
    matches = get_frameLabels(frame, movieParams)

    Parameters
    ----------
    frame : int
        Trial-relative frame index to test against label windows.
    movieParams : dict
        Movie export parameters; reads 'frameLabels' and 'persistentFrameLabels'.

    Returns
    -------
    matches : list of tuple
        (label, color) pairs active at ``frame``.
    """
    if movieParams['persistentFrameLabels']:
        matches = [(label, color) for start, end, label, color in movieParams['frameLabels'] if frame >= start]
    else:
        matches = [(label, color) for start, end, label, color in movieParams['frameLabels'] if start <= frame <= end]
    return matches

def frames_to_seconds(frame: int, fps: float, decimals: int = 1, unit: str = "s", zero_frame: int = 0) -> str:
    """
    Convert a frame number to a formatted elapsed-time string.

    The elapsed time is measured relative to ``zero_frame`` so that the frame
    equal to ``zero_frame`` reads as 0; earlier frames read negative. ``unit``
    selects the displayed units: 's' (seconds) or 'ms' (milliseconds).

    Example
    -------
    frame = 120
    fps = 22.8
    decimals = 2
    unit = 'ms'
    zero_frame = 60
    importlib.reload(jsm)
    label = frames_to_seconds(
        frame, fps=fps, decimals=decimals, unit=unit, zero_frame=zero_frame
    )

    Parameters
    ----------
    frame : int
        Frame index to convert.
    fps : float
        Acquisition frame rate (frames per second).
    decimals : int
        Number of decimal places in the formatted output.
    unit : str
        Display unit: 's' for seconds or 'ms' for milliseconds.
    zero_frame : int
        Frame index treated as t = 0; subtracted from ``frame`` before
        conversion so the timeline can be zeroed to a user-defined frame.

    Returns
    -------
    label : str
        Formatted time string (value plus unit), e.g. '0.50 s' or '-44.84 ms'.
    """
    seconds = (frame - zero_frame) / fps
    value = seconds * 1000.0 if unit == "ms" else seconds
    return f"{value:.{decimals}f} {unit}"
def add_image_scaleBar(ax, scaleBarParams, color = (1,1,1), includeLabel = False):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    sizeX = xlim[1] - xlim[0]
    sizeY = ylim[1] - ylim[0]
    if not 'orientation' in scaleBarParams:
        scaleBarParams['orientation'] = 'horz'
    vertRatio=sizeY*scaleBarParams['vertAdjust']
    # Handle inverted x-axes (xlim high->low) while preserving "visual corner" behavior.
    xdir = 1 if sizeX >= 0 else -1
    horzRatio=abs(sizeX)*scaleBarParams['horzAdjust']
    if scaleBarParams['orientation'] == 'horz':
        if scaleBarParams['corner'] == 'BL':
            x0 = xlim[0] + xdir*horzRatio
            scaleBarParams['xcoords']=[x0, x0 + xdir*scaleBarParams['length_px']]
            scaleBarParams['ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio]
            scaleBarParams['labelcoord']=[(scaleBarParams['xcoords'][0]+scaleBarParams['xcoords'][1])/2,ylim[0]+vertRatio*1.5]
            ha='center'
            va='bottom'
        elif scaleBarParams['corner'] == 'BR':
            x0 = xlim[1] - xdir*horzRatio
            scaleBarParams['xcoords']=[x0, x0 - xdir*scaleBarParams['length_px']]
            scaleBarParams['ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio]
            scaleBarParams['labelcoord']=[(scaleBarParams['xcoords'][0]+scaleBarParams['xcoords'][1])/2,ylim[0]+vertRatio*1.5]
            ha='center'
            va='bottom'
        elif scaleBarParams['corner'] == 'TL':
            x0 = xlim[0] + xdir*horzRatio
            scaleBarParams['xcoords']=[x0, x0 + xdir*scaleBarParams['length_px']]
            scaleBarParams['ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio]
            scaleBarParams['labelcoord']=[(scaleBarParams['xcoords'][0]+scaleBarParams['xcoords'][1])/2,ylim[1]-vertRatio*1.5]
            ha='center'
            va='top'
        elif scaleBarParams['corner'] == 'TR':
            x0 = xlim[1] - xdir*horzRatio
            scaleBarParams['xcoords']=[x0, x0 - xdir*scaleBarParams['length_px']]
            scaleBarParams['ycoords']=[ylim[1]-vertRatio,ylim[1]-vertRatio]
            scaleBarParams['labelcoord']=[(scaleBarParams['xcoords'][0]+scaleBarParams['xcoords'][1])/2,ylim[1]-vertRatio*1.5]
            ha='center'
            va='top'
    elif scaleBarParams['orientation'] == 'vert':
        if scaleBarParams['corner'] == 'BL':
            x0 = xlim[0] + xdir*horzRatio
            scaleBarParams['xcoords']=[x0,x0]
            scaleBarParams['ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio-scaleBarParams['length_px']]
            scaleBarParams['labelcoord']=[xlim[0]+xdir*horzRatio*1.5,ylim[0]+vertRatio-scaleBarParams['length_px']/2]
            ha='left'
            va='center'
        elif scaleBarParams['corner'] == 'BR':
            x0 = xlim[1] - xdir*horzRatio
            scaleBarParams['xcoords']=[x0,x0]
            scaleBarParams['ycoords']=[ylim[0]+vertRatio,ylim[0]+vertRatio-scaleBarParams['length_px']]
            scaleBarParams['labelcoord']=[xlim[1]-xdir*horzRatio*1.5,ylim[0]+vertRatio-scaleBarParams['length_px']/2]
            ha='right'
            va='center'
        elif scaleBarParams['corner'] == 'TL':
            x0 = xlim[0] + xdir*horzRatio
            scaleBarParams['xcoords']=[x0,x0]
            scaleBarParams['ycoords']=[ylim[1]+vertRatio,ylim[1]+vertRatio+scaleBarParams['length_px']]
            scaleBarParams['labelcoord']=[xlim[0]+xdir*horzRatio*1.5,ylim[1]+vertRatio+scaleBarParams['length_px']/2]
            ha='left'
            va='center'
        elif scaleBarParams['corner'] == 'TR':
            x0 = xlim[1] - xdir*horzRatio
            scaleBarParams['xcoords']=[x0,x0]
            scaleBarParams['ycoords']=[ylim[1]+vertRatio,ylim[1]+vertRatio+scaleBarParams['length_px']]
            scaleBarParams['labelcoord']=[xlim[1]-xdir*horzRatio*1.5,ylim[1]+vertRatio+scaleBarParams['length_px']/2]
            ha='right'
            va='center'
    # print(scaleBarParams['xcoords'])
    ax.plot(scaleBarParams['xcoords'],scaleBarParams['ycoords'],linestyle = '-',color=color,linewidth=scaleBarParams['lw'],solid_capstyle='butt')
    # 'displayLabel' (default True) is a master on/off for the scale-bar text. When
    # False, the text is suppressed regardless of the per-call/panel-based includeLabel
    # or the scaleBarParams['includeLabel'] flag (which still select WHICH panel(s)
    # would otherwise show it).
    if scaleBarParams.get('displayLabel', True) and (scaleBarParams['includeLabel'] or includeLabel):
        ax.text(scaleBarParams['labelcoord'][0],scaleBarParams['labelcoord'][1],\
                str(scaleBarParams['length_um'])+' μm',color=color,fontsize=scaleBarParams['fontsize'],ha=ha,va=va)
    return ax,scaleBarParams
##################################################################################################
##################################################################################################
def draw_ROI_borders(ax,ROI_borders,ROI_colors, ls = '-', lw = 1, alpha = 1,xadjust=0,yadjust=0):
    for r,ROI in enumerate(ROI_borders):
        for b in ROI_borders[r].keys():
            ax.plot(ROI_borders[r][b]['border'][0,:]-xadjust,ROI_borders[r][b]['border'][1,:]-yadjust,ls = ls,lw = lw,alpha = alpha, color=ROI_colors[r])
    return ax
def draw_crop_box(ax, crop_coords, lw=1, ls=':', color=(1,1,1), alpha=1, mode=2, dashes=None):
    x1 = crop_coords['x1']
    y1 = crop_coords['y1']
    x2 = crop_coords['x2']
    y2 = crop_coords['y2']
    w = np.absolute(x2 - x1)
    h = np.absolute(y2 - y1)
    if mode == 1:
        lines = [
            ax.plot([x1,x1], [y1,y2], ls, lw=lw, color=color, alpha=alpha)[0],
            ax.plot([x2,x2], [y1,y2], ls, lw=lw, color=color, alpha=alpha)[0],
            ax.plot([x1,x2], [y1,y1], ls, lw=lw, color=color, alpha=alpha)[0],
            ax.plot([x1,x2], [y2,y2], ls, lw=lw, color=color, alpha=alpha)[0],
        ]
        if dashes is not None:
            for line in lines:
                line.set_dashes(dashes)
    elif mode == 2:
        rect = Rectangle((x1, y1), w, h, ls=ls, lw=lw, edgecolor=color, alpha=alpha, facecolor='none')
        if dashes is not None:
            rect.set_dashes(dashes)
        ax.add_patch(rect)
    return ax

