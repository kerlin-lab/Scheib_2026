# This contains all the supporting functions utilized in the Scheib et. al., 2026 manuscript


import jsimg as jsi
import ttracking as tt
###################################################
import numpy as np
from tqdm.auto import tqdm
import copy
import statsmodels.api as sm
from scipy import stats
import time
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import skew,mode,gmean,ranksums,ttest_ind, zscore, kstest
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["Arial"]
plt.rcParams['pdf.fonttype'] = 42
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
                #         conglom2 = np.hstack([conglom, np.ones([conglom.shape[0], frMissing])*np.NaN])
                #     else:
                #         conglom2 = np.hstack([conglom, np.ones([conglom.shape[0], frMissing])*np.NaN])
                #     aligned.append(conglom2)
                # else:
                if not downSample:
                    aligned.append(np.ones([tsep.shape[1], nFr])*np.NaN)  
                else:
                    aligned.append(np.ones([tsep.shape[1], int(nFr/downSampleFactor)])*np.NaN)  
        else:
            if not downSample:
                aligned.append(np.ones([tsep.shape[1], nFr])*np.NaN)  
            else:
                aligned.append(np.ones([tsep.shape[1], int(nFr/downSampleFactor)])*np.NaN)  
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
                    aligned.append(np.ones([tsep.shape[1], totalFrames])*np.NaN)
                else:
                    aligned.append(np.ones([tsep.shape[1], expectedBins])*np.NaN)
            else:
                if not downSample:
                    pre = tsep[t,:,time-preFrames:time]
                    post = []
                    for tp in range(time, time+postFrames):
                        if tp < tsep.shape[2]:
                            post.append(tsep[t,:,tp])
                        else:
                            post.append(np.ones(tsep.shape[1]) * np.NaN)
                    aligned.append(np.vstack([list(pre.T),post]).T)
                else:
                    pre = []
                    for tp in range(time-preFrames, time, downSampleFactor):
                        pre.append(np.nanmean(tsep[t,:,tp:tp + downSampleFactor], axis=1))
                    post = []
                    for tp in range(time, time+postFrames, downSampleFactor):
                        if tp+downSampleFactor<tsep.shape[2]:
                            post.append(np.nanmean(tsep[t,:,tp:tp + downSampleFactor], axis=1))
                        else:
                            post.append(np.ones(tsep.shape[1])*np.NaN)
                    aligned.append(np.vstack([pre, post]).T)
        else:
            if not downSample:
                aligned.append(np.ones((tsep.shape[1], totalFrames))*np.NaN)
            else:
                aligned.append(np.ones((tsep.shape[1], expectedBins))*np.NaN)
    return np.array(aligned)

def get_secondLickSide(beh, returnTime = False):
    bpodLicks = tt.get_lickTimes(beh)
    cTimes = tt.get_cTimes2(beh)
    licks = tt.get_contacts(beh['trajs'])
    azi,ele = tt.get_aziEle(beh, licks)
    tracker = tt.get_tracker(beh)
    goCues = tt.get_goCues(beh)
    sLickTrajs = tt.get_consumption_licks(beh['trajs'],0)
    sLickTimes = sLickTrajs[:,:,3]
    sLick = np.ones(tracker.shape[0])*np.NaN
    slTimes = np.ones(tracker.shape[0])*np.NaN
    
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
    licks = tt.get_contacts(beh['trajs'])
    azi,ele = tt.get_aziEle(beh, licks)
    tracker = tt.get_tracker(beh)
    me, de = tt.get_errors_LDA(azi, ele, tracker)
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
    bpodLicks = tt.get_lickTimes(beh)
    cTimes = tt.get_cTimes2(beh)
    licks = tt.get_contacts(beh['trajs'])
    azi,ele = tt.get_aziEle(beh, licks)
    tracker = tt.get_tracker(beh)
    goCues = tt.get_goCues(beh)
    sLickTrajs = tt.get_consumption_licks(beh['trajs'],0)
    sLickTimes = sLickTrajs[:,:,3]
    sLick = np.ones(tracker.shape[0])*np.NaN
    slTimes = np.ones(tracker.shape[0])*np.NaN
    slNHits = np.ones(tracker.shape[0])*np.NaN
    nextIdent = np.ones(tracker.shape[0])*np.NaN
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
                projs.append([np.NaN, np.NaN])
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
    #     irfs = np.ones([tlength,2])*np.NaN
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
                FL = np.ones_like(CL)*np.NaN
    
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
                DE = np.ones_like(CR)*np.NaN
            if hasME:
                if binary:
                    ME = np.nanmean(psTsep[MEt,:,:]>0, axis = 0)
                else:
                    ME = np.nanmean(psTsep[MEt,:,:], axis = 0)
            else:
                ME = np.ones_like(CR)*np.NaN
            
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
            FL = np.ones_like(CL)*np.NaN

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
            ME = np.ones_like(CL)*np.NaN

        if np.sum(rDEm)>=5:
            if binary:
                DE = np.nanmean(psTsep[rDEm,:,:]>0, axis = 0)
            else:
                DE = np.nanmean(psTsep[rDEm,:,:], axis = 0)
        else:
            DE = np.ones_like(CL)*np.NaN

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
            NSH = np.ones_like(CR)*np.NaN

        if hasHit:
            if binary:
                HIT = np.nanmean(psTsep[hitT,:,:]>0, axis = 0)
            else:
                HIT = np.nanmean(psTsep[hitT,:,:], axis = 0)
        else:
            HIT = np.ones_like(CR)*np.NaN
        
        if hasDE:
            if binary:
                DE = np.nanmean(psTsep[DEt,:,:]>0, axis = 0)
            else:
                DE = np.nanmean(psTsep[DEt,:,:], axis = 0)
        else:
            DE = np.ones_like(CR)*np.NaN
        
        if hasME:
            if binary:
                ME = np.nanmean(psTsep[MEt,:,:]>0, axis = 0)
            else:
                ME = np.nanmean(psTsep[MEt,:,:], axis = 0)
        else:
            ME = np.ones_like(CR)*np.NaN
            
        if hasAR:
            if binary:
                AR = np.nanmean(psTsep[ARt,:,:]>0, axis = 0)
            else:
                AR = np.nanmean(psTsep[ARt,:,:], axis = 0)
        else:
            AR = np.ones_like(CR)*np.NaN

        if hasAL_ME:
            if binary:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:]>0, axis = 0)
            else:
                AL_ME = np.nanmean(psTsep[al_MEt,:,:], axis = 0)
        else:
            AL_ME = np.ones_like(CR)*np.NaN

        if hasAL_noME:
            if binary:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:]>0, axis = 0)
            else:
                AL_noME = np.nanmean(psTsep[al_noMEt,:,:], axis = 0)
        else:
            AL_noME = np.ones_like(CR)*np.NaN

        if hasAN:
            if binary:
                AN = np.nanmean(psTsep[ANt,:,:]>0, axis = 0)
            else:
                AN = np.nanmean(psTsep[ANt,:,:], axis = 0)
        else:
            AN = np.ones_like(CR)*np.NaN

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
        DE = np.ones_like(CR)*np.NaN

    if np.sum(rMEm)>=5:
        if binary:
            ME = np.nanmean(psTsep[rMEm,:,:]>0, axis = 0)
        else:
            ME = np.nanmean(psTsep[rMEm,:,:], axis = 0)
    else:
        ME = np.ones_like(CR)*np.NaN

    if np.sum(ARm)>=5:
        if binary:
            AR = np.nanmean(psTsep[ARm,:,:]>0, axis = 0)
        else:
            AR = np.nanmean(psTsep[ARm,:,:], axis = 0)
    else:
        AR = np.ones_like(CR)*np.NaN

    if np.sum(ALnoMEm)>=5:
        if binary:
            AL_noME = np.nanmean(psTsep[ALnoMEm,:,:]>0, axis = 0)
        else:
            AL_noME = np.nanmean(psTsep[ALnoMEm,:,:], axis = 0)
    else:
        AL_noME = np.ones_like(CR)*np.NaN

    if np.sum(AL_MEm)>=5:
        if binary:
            AL_ME = np.nanmean(psTsep[AL_MEm,:,:]>0, axis = 0)
        else:
            AL_ME = np.nanmean(psTsep[AL_MEm,:,:], axis = 0)
    else:
        AL_ME = np.ones_like(CR)*np.NaN

    if np.sum(ANm)>=5:
        if binary:
            AN = np.nanmean(psTsep[ANm,:,:]>0, axis = 0)
        else:
            AN = np.nanmean(psTsep[ANm,:,:], axis = 0)
    else:
        AN = np.ones_like(CR)*np.NaN

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
    ortho=False, orthoOrder=False, fixedCDs=False, dffKey = 'consensus_NMFtc_dff_bc_decon_sp_events',useDS = True, binary = False):


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
        CR = np.ones_like(CR)*np.NaN
        CL = np.ones_like(CR)*np.NaN
        IR = np.ones_like(CR)*np.NaN
        IL = np.ones_like(CR)*np.NaN
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
            CR.append(cd*psTsep[trial,:,:])
        CR=np.nanmean(np.dstack(CR),axis=2)

        CL=[]
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
            CL.append(cd*psTsep[trial,:,:])
        CL=np.nanmean(np.dstack(CL),axis=2)

        IR=[]
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
            IR.append(cd*psTsep[trial,:,:])
        IR=np.nanmean(np.dstack(IR),axis=2)

        IL=[]
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
            IL.append(cd*psTsep[trial,:,:])
        IL=np.nanmean(np.dstack(IL),axis=2)
  
    tsep3 = np.array([CR, CL, IR, IL])
        
    mT_preTType.append(tsep3)

    mT_preTType = np.dstack(mT_preTType)
    return mT_preTType

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

def get_SMO_fast(megaTsepsS, megaTsepsM, megaTsepsO, nBoot, anat, mod, calcPhase='pre1', center = False):
    allData = {'projs': {'pre': []}}
  
    for i in range(nBoot):
        s  = megaTsepsS["pre"][mod][anat][:,:,:,i]
        m  = megaTsepsM["pre"][mod][anat][:,:,:,i]
        o  = megaTsepsO["pre"][mod][anat][:,:,:,i]
        
        nTime = s.shape[1]
        
        allD = [s, m, o]
        if center:
            for d in allD:
                d = d-np.nanmean(d, axis = 1, keepdims=True)
                
        for d in allD:
            d[np.isnan(d)] = 0

        data = [s,m,o]
        projs = []
        for a in range(3):
            axProjs = {
                'pre': np.nanmean(data[a],axis=2)
            }
            projs.append(axProjs)
        
        for key in list(axProjs.keys()):
            allData['projs'][key].append(np.dstack([projs[a][key] for a in range(3)]))

    for key in list(allData['projs'].keys()):
        allData['projs'][key] = np.stack(allData['projs'][key])
    
    return allData


def get_CRMER_CD2_fast(megaTsepsHH, anat, mod, calcPhase='post1', square = False, cdTrials = [0,3], nBoot = 1000):
    def proj(mat, axes, sqr = True):
        prjn = np.einsum('gtr,tr->gt', mat, axes)
        if sqr:
            prjn = np.sqrt(abs(prjn))*np.sign(prjn)
        return prjn
    allData = {'projs': {'pre1': [], 'pre2': [], 'post1': [], 'post2': []}, 
               'axes': [], 'axesOrig': []}
    
    for i in range(nBoot):
        post1 = megaTsepsHH["post1"][f'err{mod}'][anat][:,:,:,i]
        post2 = megaTsepsHH["post2"][f'err{mod}'][anat][:,:,:,i]
        pre1  = megaTsepsHH["pre1"][f'err{mod}'][anat][:,:,:,i]
        pre2  = megaTsepsHH["pre2"][f'err{mod}'][anat][:,:,:,i]

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
        post1 = megaTsepsHH["post1"][mod][anat][:,:,:,i]
        post2 = megaTsepsHH["post2"][mod][anat][:,:,:,i]
        pre1  = megaTsepsHH["pre1"][mod][anat][:,:,:,i]
        pre2  = megaTsepsHH["pre2"][mod][anat][:,:,:,i]
        # post1 = megaTsepsHH["early1"][mod][anat][:,:,:,i]
        # post2 = megaTsepsHH["early2"][mod][anat][:,:,:,i]

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
        early1 = megaTsepsHH["early1"][mod][anat][:,:,:,i]
        early2 = megaTsepsHH["early2"][mod][anat][:,:,:,i]
        late1 = megaTsepsHH["late1"][mod][anat][:,:,:,i]
        late2 = megaTsepsHH["late2"][mod][anat][:,:,:,i]
        # pre1 = megaTsepsHH["pre1"][mod][anat][:,:,:,i]
        # pre2 = megaTsepsHH["pre2"][mod][anat][:,:,:,i]
        pre1 = megaTsepsHH["preShift_late1"][mod][anat][:,:,:,i]
        pre2 = megaTsepsHH["preShift_late2"][mod][anat][:,:,:,i]

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







# Wrapper functions

def get_allAnmParams(masterData, anms, frs, factors, Z, phases2, dffKeys, windowSize = [-3.6,3.6]):
    allAnmParams = {}
    for anmIdx,anm in enumerate(anms):
        print("=========================================================================================================")
        print(f'{anmIdx} {anm}')

        params = {}
        anmData = masterData[anm]
        beh = anmData['behavior']
        
        params['zavail'] = np.array(anmData['behavior']['twoCams'])
        dendDays = jsi.get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = tt.get_goCues(anmData['behavior'])
        cTimes = tt.get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = tt.get_lastLick(anmData['behavior'])
        params['lastLicks'] = lastLicks
        slSides, slTimes = get_secondLickSide(beh, returnTime = True)
        params['sTimes'] = slTimes
        slTimes2 = slTimes+goCues

        tracker = tt.get_tracker(anmData['behavior'])
        params['tracker'] = tracker
        
        try:
            params['aw'] = tt.get_autoWaterTrials(anmData['behavior'])
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
                
                #tsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, Z = Z)#, gatherPreTrial = True)
                preTsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
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

        cLicks = tt.get_contacts(anmData['behavior']['trajs'])
        cazi,cele = tt.get_aziEle(anmData['behavior'], cLicks)
        cdists, lda = tt.get_relDistances(cazi, cele, tracker, returnLDA = True)
        params['cLicks'] = cLicks
        params['cDists'] = cdists
        params['cAzi'] = cazi
        params['cEle'] = cele
        params['lda'] = lda
        

        fLicks = tt.get_firsts(anmData['behavior']['trajs'])
        fazi,fele = tt.get_aziEle(anmData['behavior'], fLicks)
        fdists2 = tt.get_relDistances(fazi, fele, tracker) ########################## need to check if lda is redefined here
        fdists = get_relDistFromLDA(fazi, fele, lda)
        params['fLicks'] = fLicks
        params['fDists_orig'] = fdists2
        params['fDists'] = fdists
        params['fAzi'] = fazi
        params['fEle'] = fele



        params['shiftMask'] = {}
        # for phase in phases2:
        #     params['shiftMask'][phase] = tt.get_behaviorShiftMask(tracker, phase)
        relDays = np.unique(tracker[:,0])
        print(f'relDays: {relDays}')
        for phase in phases2:
            if phase == 'late':
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
                params['shiftMask'][phase] = tt.get_behaviorShiftMask(tracker, phase)

        
        me, de = tt.get_errors_LDA(cazi, cele, tracker, returnLDA=False)
        minDists = tt.get_minDistToPort(beh, 'contacts', 'post_shift', 'left')
        params['me'] = me
        params['de'] = de
        params['minDist'] = minDists

        ar, al, an = get_MER_bySecondLick(beh)
        params['afterRight'] = ar
        params['afterLeft'] = al
        params['afterNone'] = an
        
        sLicks = tt.get_consumption_licks(beh['trajs'],0)
        aziS,eleS = tt.get_aziEle(beh, sLicks)
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

        actual = np.ones(nanMask.shape[0])*np.NaN
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


def get_allAnmParams2(masterData, anms, frs, factors, Z, phases2, dffKeys, windowSize = [-3.6,3.6]):
    allAnmParams = {}
    for anmIdx,anm in enumerate(anms):
        print("=========================================================================================================")
        print(f'{anmIdx} {anm}')

        params = {}
        anmData = masterData[anm]
        beh = anmData['behavior']
        
        params['zavail'] = np.array(anmData['behavior']['twoCams'])
        dendDays = jsi.get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = tt.get_goCues(anmData['behavior'])
        cTimes = tt.get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = tt.get_lastLick(anmData['behavior'])
        params['lastLicks'] = lastLicks
        slSides, slTimes = get_secondLickSide(beh, returnTime = True)
        params['sTimes'] = slTimes
        slTimes2 = slTimes+goCues

        tracker = tt.get_tracker(anmData['behavior'])
        params['tracker'] = tracker

        try:
            params['aw'] = tt.get_autoWaterTrials(anmData['behavior'])
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
                
                #tsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, Z = Z)#, gatherPreTrial = True)
                preTsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
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

        cLicks = tt.get_contacts(anmData['behavior']['trajs'])
        cazi,cele = tt.get_aziEle(anmData['behavior'], cLicks)
        cdists, lda = tt.get_relDistances(cazi, cele, tracker, returnLDA = True)
        params['cLicks'] = cLicks
        params['cDists'] = cdists
        params['cAzi'] = cazi
        params['cEle'] = cele
        params['lda'] = lda
        

        fLicks = tt.get_firsts(anmData['behavior']['trajs'])
        fazi,fele = tt.get_aziEle(anmData['behavior'], fLicks)
        fdists2 = tt.get_relDistances(fazi, fele, tracker) ########################## need to check if lda is redefined here
        fdists = get_relDistFromLDA(fazi, fele, lda)
        params['fLicks'] = fLicks
        params['fDists_orig'] = fdists2
        params['fDists'] = fdists
        params['fAzi'] = fazi
        params['fEle'] = fele



        params['shiftMask'] = {}
        # for phase in phases2:
        #     params['shiftMask'][phase] = tt.get_behaviorShiftMask(tracker, phase)
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
                params['shiftMask'][phase] = tt.get_behaviorShiftMask(tracker, phase)

        
        me, de = tt.get_errors_LDA(cazi, cele, tracker, returnLDA=False)
        minDists = tt.get_minDistToPort(beh, 'contacts', 'post_shift', 'left')
        params['me'] = me
        params['de'] = de
        params['minDist'] = minDists

        ar, al, an = get_MER_bySecondLick(beh)
        params['afterRight'] = ar
        params['afterLeft'] = al
        params['afterNone'] = an
        
        sLicks = tt.get_consumption_licks(beh['trajs'],0)
        aziS,eleS = tt.get_aziEle(beh, sLicks)
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

        actual = np.ones(nanMask.shape[0])*np.NaN
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


def update_allAnmParams(allAnmParams, masterData, anms, frs, factors, Z, dffKeys,  windowSize = [-3.6,3.6]):

    for anm in anms:
        
        params = {}
        anmData = masterData[anm]

        dendDays = jsi.get_dendDays(anmData)
        params['dendDays'] = dendDays
        if np.sum(dendDays)>0:
            if np.sum(~dendDays)>0:
                anmA = ['dendrites', 'somas']
            else:
                anmA = ['dendrites']
        else:
            anmA = ['somas']

        goCues = tt.get_goCues(anmData['behavior'])
        cTimes = tt.get_cTimes2(anmData['behavior'])
        params['goCues'] = goCues
        params['cTimes'] = cTimes
        lastLicks = tt.get_lastLick(anmData['behavior'])
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
                
                #tsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, Z = Z)#, gatherPreTrial = True)
                preTsep = jsi.get_tsep(anmData[anat], dffKey, nanMaskFrs = True, gatherPreTrial = True, Z = Z)
                preTsep = preTsep[:,:,int(fr*3):]
                # params2['fullTsep'][anat] = {}
                # params2['fullTsep'][anat]['full'] = preTsep
                
                # ctps = psTsep_func(preTsep, cTimes[anatMask]+1, fr, windowSize = [-3.6,3.6])#, downSample = True, downSampleFactor = factors[ak])
                # gops = psTsep_func(preTsep, goCues[anatMask]+1, fr, windowSize = [-3.6,3.6])#, downSample = True, downSampleFactor = factors[ak])
                # sps = psTsep_func(preTsep, slTimes2[anatMask]+1, fr, windowSize = [-3.6,3.6])
                # params2['fullTsep'][anat]['ct-psTsep'] = ctps
                # params2['fullTsep'][anat]['go-psTsep'] = gops
                # params2['fullTsep'][anat]['sl-psTsep'] = sps

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
            # goodRois = jsi.getGoodRois(anmData, dffKey = 'NMF')
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


    # iAnms = {}
    # aAnms = {}
    # minTrials = 100
    # allMinTrials = []
    # anmsMinTrial = {}
    # for anat in anatKeys:
    #     iAnms[anat] = []
    #     aAnms[anat] = []
    #     anmsMinTrial[anat] = {}
    #     for anm in anms:
    #         if anat in list(allAnmParams[anm][dffKey]['dsTsep'].keys()):
    #             aAnms[anat].append(anm)
    #             tracker = allAnmParams[anm]['tracker']
    #             dendDays = allAnmParams[anm]['dendDays']
    #             preMask = allAnmParams[anm]['shiftMask']['pre']
                
    #             if anat == 'dendrites':
    #                 anatMask = dendDays
    #             else:
    #                 anatMask = ~dendDays

    #             preAnatMask = np.logical_and(anatMask, preMask)
    #             preCR = np.sum(tracker[preAnatMask,2]==0)
    #             preCL = np.sum(tracker[preAnatMask,2]==1)
    #             preFR = np.sum(tracker[preAnatMask,2]==5)
    #             preFL = np.sum(tracker[preAnatMask,2]==6)

    #             if np.all(np.hstack([item>minTrials for item in [preCR, preCL]])):
    #                 if anm != 'B00002121777':
    #                     if tsep.shape[1]>=10:
    #                         iAnms[anat].append(anm)
    #                         print(anat, anm, preCR, preCL)
    #                         allMinTrials.append(np.hstack([preCR, preCL]))
    #                         anmsMinTrial[anat][anm] = np.min(np.hstack([preCR, preCL]))
    # print(np.min(np.hstack(allMinTrials)))

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
    ##################################################################################################################################

    return aAnms, iAnms, anmsMinTrial, allMinTrials

def full_get_SMO(anatKeys, mods, phases, allAnmParams, anmsToUse, dffKey, alignTimes, frs, keepRoisDS, label, preComp, nBoot = 1000, shuffAnms = True,
             binary = False, enforceMinTrial = True, anmsMinTrial = None):
    
    megaTsepsS = {}
    megaTsepsM = {}
    megaTsepsO = {}

    megaTsepsHH = {}
    emptyMegaTsep = {'dendrites': [], 'somas': []}

    for p, phase in enumerate(phases):
        megaTsepsS[phase] = {}
        megaTsepsM[phase] = {}
        megaTsepsO[phase] = {}
        megaTsepsHH[f'{phase}1'] = {}
        megaTsepsHH[f'{phase}2'] = {}
        for mod in mods:
            megaTsepsS[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsM[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsO[phase][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsHH[f'{phase}1'][mod] = copy.deepcopy(emptyMegaTsep)
            megaTsepsHH[f'{phase}2'][mod] = copy.deepcopy(emptyMegaTsep)

    #making holder matracies
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
                for cdKey in label:
                    mT[cdKey] = gather_anmShuff_megaTsepFullSMO(allAnmParams, anat, ak, mod, phases, frs, keepRoisDS, cdKey, dffKey, alignTimes[m], anmsToUse, nRois,
                                                    anmIdxs, preComp, enforceMinTrial = enforceMinTrial, 
                                                    binary = binary, returnErrors = returnErrors, anmsMinTrial = anmsMinTrial)
                # mT1, mT2 = gather_anmShuff_megaTsepHH(allAnmParams, anat, ak, mod, phases, frs,  keepRoisDS, dffKey, alignTimes[m], anmsToUse, nRois,
                                                    # anmIdxs, enforceMinTrial = enforceMinTrial,
                                                    # binary = binary, returnErrors = returnErrors, anmsMinTrial = anmsMinTrial)
                for phase in phases:
                    megaTsepsS[phase][mod][anat][:,:,:,i] = mT['s'][phase]
                    megaTsepsM[phase][mod][anat][:,:,:,i] = mT['m'][phase]
                    megaTsepsO[phase][mod][anat][:,:,:,i] = mT['o'][phase]
                    # megaTsepsHH[f'{phase}1'][mod][anat][:,:,:,i] = mT1[phase]
                    # megaTsepsHH[f'{phase}2'][mod][anat][:,:,:,i] = mT2[phase]

    return megaTsepsS, megaTsepsM, megaTsepsO, megaTsepsHH



def plot_CD_projections(nRows,nCols,mod,colors,labels,niceAnatLabels,plot_scaleBar, plotTtypes, plotTtypeLabels, dataSource, normMethod, colScaling, normFactors, ratio, CR_idx, CA_idx, AP_idx, CL_idx, CRCA_HHscroll, APCA_HHscroll, CRCL_HHscroll, anatKeys, frs, factors, window, xspace, xlim, error, CI, nBoots, sharex, sharey, gridspec_kw, vscalar, hscalar):
    nGroups = len(plotTtypes)
    delCol = 0
    if not nCols == nGroups:
        delCol=nCols-nGroups
    figSize = [nCols*hscalar,nRows*vscalar]
    fig,ax=sled_general_tools.clean_subplots(nRows,nCols,figsize=figSize,sharex=sharex,sharey=sharey,gridspec_kw=gridspec_kw,verbose=True)
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
                pText = sled_stat_tools.pstar(pVal,inclP = False)
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
            ax[row,col],plot_scaleBar = sled_ROI_tools.add_plot_scaleBar(ax[row,col],plot_scaleBar,(0,0,0),True,vertLabel)
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
#Clean Matplotlib Figures
def enforce_equal_axes(fig, axes):
    """
    Force all axes to have identical positions after layout/rendering.
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
    sharey = False, sharex = False, gridspec_kw = None, verbose = False):
    

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


    if verbose:
        display_nested_dict(info)


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
    """Turn off all clipping in axes ax; call immediately before drawing/showing"""
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
##################################################################################################
