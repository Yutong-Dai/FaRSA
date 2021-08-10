'''
File: params.py
Author: Yutong Dai (rothdyt@gmail.com)
File Created: 2019-10-31 15:51
Last Modified: 2020-12-22 17:18
--------------------------------------------
Description:
'''
params = {}
params["max_iter"] = 100000
params["max_time"] = 1000
params["printlevel"] = 2
params["printHeadEvery"] = 1000
params["printevery"] = 1000
# try this no lager [1e-10, 1e-3]
params["eta"] = 1e-3
params["zeta"] = 0.8
params["beta"] = 1 / 0.9
params["maxback"] = 100
params["stepsizeRule"] = "bkt"
params["stepsizeStrategy"] = "frac"
params["tfocs_t"] = True

fileTypeDict = {}
fileTypeDict['a9a'] = 'txt'
fileTypeDict['australian'] = 'txt'
fileTypeDict['breast_cancer'] = 'txt'
fileTypeDict['cod_rna'] = 'txt'
fileTypeDict['colon_cancer'] = 'bz2'
fileTypeDict['covtype'] = 'txt'
fileTypeDict['diabetes'] = 'txt'
fileTypeDict['duke'] = 'bz2'
fileTypeDict['fourclass'] = 'txt'
fileTypeDict['german_numer'] = 'txt'
fileTypeDict['gisette'] = 'bz2'
fileTypeDict['heart'] = 'txt'
fileTypeDict['ijcnn1'] = 'txt'
fileTypeDict['ionosphere'] = 'txt'
fileTypeDict['leu'] = 'bz2'
fileTypeDict['liver_disorders'] = 'txt'
fileTypeDict['madelon'] = 'txt'
fileTypeDict['mushrooms'] = 'txt'
fileTypeDict['phishing'] = 'txt'
fileTypeDict['rcv1'] = 'txt'
fileTypeDict['skin_nonskin'] = 'txt'
fileTypeDict['sonar'] = 'txt'
fileTypeDict['splice'] = 'txt'
fileTypeDict['svmguide1'] = 'txt'
fileTypeDict['svmguide3'] = 'txt'
fileTypeDict['w8a'] = 'txt'
fileTypeDict['avazu_app_tr'] = 'txt'
fileTypeDict['HIGGS'] = 'txt'
fileTypeDict['news20'] = 'txt'
fileTypeDict['real-sim'] = 'txt'
fileTypeDict['SUSY'] = 'txt'
fileTypeDict['url_combined'] = 'txt'

# regression
fileTypeDict['abalone_scale'] = 'txt'
fileTypeDict['bodyfat_scale'] = 'txt'
fileTypeDict['cadata'] = 'txt'
fileTypeDict['cpusmall_scale'] = 'txt'
fileTypeDict['eunite2001'] = 'txt'
fileTypeDict['housing_scale'] = 'txt'
fileTypeDict['mg_scale'] = 'txt'
fileTypeDict['mpg_scale'] = 'txt'
fileTypeDict['pyrim_scale'] = 'txt'
fileTypeDict['space_ga_scale'] = 'txt'
fileTypeDict['E2006.train'] = 'txt'
fileTypeDict['log1p.E2006.train'] = 'txt'
fileTypeDict['YearPredictionMSD'] = 'txt'
fileTypeDict['blogData_train'] = 'txt'
fileTypeDict['UJIIndoorLoc'] = 'txt'
fileTypeDict['driftData'] = 'txt'
fileTypeDict['virusShare'] = 'txt'
fileTypeDict['triazines_scale'] = 'txt'


fileTypeDict['bodyfat_scale_expanded7'] = 'mat'
fileTypeDict['pyrim_scale_expanded5'] = 'mat'
fileTypeDict['triazines_scale_expanded4'] = 'mat'
fileTypeDict['housing_scale_expanded7'] = 'mat'

gammas = {}
gammas['bodyfat_scale_expanded7'] = [1e-4, 1e-5, 1e-6]
gammas['pyrim_scale_expanded5'] = [1e-2, 1e-3, 1e-4]
gammas['triazines_scale_expanded4'] = [1e-2, 1e-3, 1e-4]
gammas['housing_scale_expanded7'] = [1e-2, 1e-3, 1e-4]
groups = {}
groups['bodyfat_scale_expanded7'] = 388
groups['pyrim_scale_expanded5'] = 671
groups['triazines_scale_expanded4'] = 2118
groups['housing_scale_expanded7'] = 258
