from fafbseg import flywire
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#go to this link to set flywire secret: https://globalv1.flywire-daf.com/auth/api/v1/refresh_token
flywire.set_chunkedgraph_secret("92d0b9c2bf9d8989cd25096ca28b6aa9")

def get_inputs(flywire_ID_list, threshold=5):
    all_inputs = []
    all_synapses = []
    for i, flywire_ID in enumerate(flywire_ID_list):
        print('Processing ID ' + str(i) + ' out of ' + str(len(flywire_ID_list)) + ' flywire IDs')
        syn = flywire.fetch_connectivity(flywire_ID)
        inputs = syn[syn['post'] == flywire_ID]
        if len(inputs) > 0:
            inputs = inputs[inputs['weight'] >= 5]
            synapses = list(inputs['weight'])
            inputs = list(inputs['pre'])
            all_inputs = all_inputs + inputs
            all_synapses = all_synapses + synapses

    return all_inputs, all_synapses

def get_outputs(flywire_ID_list, threshold=5):
    all_outputs = []
    all_synapses = []
    for i, flywire_ID in enumerate(flywire_ID_list):
        print('Processing ID ' + str(i) + ' out of ' + str(len(flywire_ID_list)) + ' flywire IDs')
        syn = flywire.fetch_connectivity(flywire_ID)
        outputs = syn[syn['pre'] == flywire_ID]
        if len(outputs) > 0:
            outputs = outputs[outputs['weight'] >= 5]
            synapses = list(outputs['weight'])
            outputs = list(outputs['post'])
            all_outputs = all_outputs + outputs
            all_synapses = all_synapses + synapses

    return all_outputs, all_synapses

# ORDER PRESERVING, as opposed to flywire's out of the box update function
def update_ids(ids):
    return [flywire.update_ids(id)['new_id'][0] for id in ids]

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def NormalizeData(data):
    return data/np.sum(data)

def get_cmap(n, name='rainbow'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

# Internal helper function for obtaining Neuprint data

datapath="../Data/"

#returns 8 things:
#Win, the matrix of incoming connections (post=rows, pre=columns)
#Wout, the matrix of outgoing connections (post=rows, pre=columns)
#c, dataframe containing information about the neurons (same order as rows/columns of Win/Wout)
#cin, dataframe containing information about the upstream neurons (same order as columns of Win)
#cout, dataframe containing information about the downstream neurons (same order as rows of Wout)
#conns, a list with conns[i] containing the information about all connections of neuron i
#syns_out, a list with syns_out[i] containing the information about outgoing connections of neuron i
#syns_in, a list with syns_in[i] containing the information about incoming connections of neuron i
def genW2(filestr):
    print("loading from file")
    c,conns,syns_out,syns_in,out_syntots = pickle.load(open(filestr+".pkl","rb"))
    N = len(c) #number of neurons

    print("processing connections")
    #add information about the matching neuron to each connection dataframe
    for ii in range(N):
        #rename the column containing synapse counts to 'syncount'
        conns[ii].rename(columns={conns[ii].columns[3]: 'syncount'},inplace=True)
        #add data about the neuron that this connection involves
        conns[ii]["neuronbodyId"] = c.bodyId[ii] 
        conns[ii]["neuronname"] = c.name[ii] 

    conns_flat = pd.concat(conns,sort=False)
    conns_flat.loc[conns_flat.name == 0,"name"] = ""

    conns_in = conns_flat[conns_flat.relation == "upstream"].reset_index(drop=True)
    conns_out = conns_flat[conns_flat.relation == "downstream"].reset_index(drop=True)

    M_in = len(conns_in.bodyId.unique()) #number of unique neurons upstream of KCs
    M_out= len(conns_out.bodyId.unique()) #number of unique neurons downstream of KCs

    print("generating weight matrices")
    conns_in_unique = conns_in.drop_duplicates("bodyId")
    conns_out_unique = conns_out.drop_duplicates("bodyId")
    conns_in_unique = conns_in_unique.sort_values(by="name")
    conns_out_unique = conns_out_unique.sort_values(by="name")
    conns_in_unique.reset_index(drop=True,inplace=True)
    conns_out_unique.reset_index(drop=True,inplace=True)

    #generate hash functions that convert bodyIds into the index that will be used for the weight matrix
    neuron_hash = dict(zip(c.bodyId,np.arange(N))) #goes from neuron bodyId to index
    in_hash = dict(zip(conns_in_unique.bodyId,np.arange(M_in))) #goes from upstream bodyId to index
    out_hash = dict(zip(conns_out_unique.bodyId,np.arange(M_out))) #goes from downstream bodyId to index

    #make Win
    postinds = np.vectorize(neuron_hash.get)(conns_in.neuronbodyId)
    preinds = np.vectorize(in_hash.get)(conns_in.bodyId)

    Win = np.zeros([N,M_in])
    Win[postinds,preinds] = conns_in.syncount

    #make Wout
    postinds = np.vectorize(out_hash.get)(conns_out.bodyId)
    preinds = np.vectorize(neuron_hash.get)(conns_out.neuronbodyId)

    Wout = np.zeros([M_out,N])
    Wout[postinds,preinds] = conns_out.syncount

    cin = conns_in_unique[["bodyId","name"]]
    cout = conns_out_unique[["bodyId","name"]]
    
    #add postsynaptic total to downstream neurons
    #out_syntots_flat = pd.concat(out_syntots,sort=False)
    #out_syntots_flat = out_syntots_flat.drop_duplicates()
    #postinds = np.vectorize(out_hash.get)(out_syntots_flat.bodyId)

    #cout["tot"] = None
    #cout.tot[postinds] = out_syntots_flat.post.values.astype(int)

    print("done")

    return Win,Wout,c,cin,cout,conns,syns_out,syns_in

# Get Neuprint weights and glomerulus descripters as a dataframe
def get_Neuprint():
	Win,Wout,c,cin,cout,conns,syns_out,syns_in = genW2("../Data/kcs_nosyn")
	N = len(c)

	### annotate known inputs to KCs based on Greg's data and data about visual inputs ###

	#Greg's data
	#read info about uPN bodyIds and the glomeruli they belong to
	upns = pd.read_csv(datapath+"jefferis_glomerular_data_final/FIB_uPNs.csv", sep=',')
	validinds = upns.skid.isin(cin.bodyId)
	print("found",np.sum(validinds),"of",len(upns),"uPN IDs upstream of KCs")
	# print("\nPNs in Greg's data not found in EM:")
	# print(*upns[~validinds].skid.values,sep='\n')

	upns = upns[validinds] #only keep uPNs that exist in the list of neurons upstream of KCs

	#do same for mPNs and VP-PNs
	mpns = pd.read_csv(datapath+"jefferis_glomerular_data_final/FIB_mPNs.csv", sep=',')
	vppns = pd.read_csv(datapath+"jefferis_glomerular_data_final/FIB_VP_PNs.csv", sep=',')
	#vppns.val = np.array(["thermo_hygro" for i in range(len(vppns))])
	#read info about the valence of each glomerulus, add this to the above data
	gloms_val = pd.read_csv(datapath+"jefferis_glomerular_data_final/uPNsdf_FIB_FAFB.csv", sep=',')
	Nglom = len(gloms_val)

	upns["val"] = "" #add a column to the list of uPNs that contains their valence
	upns["priority"] = 1 #this will be used to sort which come first when saving the data

	# print(*upns[~validinds].val.values,sep='\n')



	allpns = pd.concat([upns,mpns,vppns],sort=False,ignore_index=True)

	#visual neuron data
	visual = pd.read_csv(datapath+"visual_inputs_Jack.csv")#pd.read_excel(datapath+"li_visual/v3.xlsx")
	#visual.rename(columns={"Body ID": "bodyId"},inplace=True)
	visual["priority"] = 4

	#generate subset of weight matrix and cin corresponding to inputs with valence label
	all_val = pd.concat([visual.bodyId,allpns.skid]) #all bodyIds that have a valence label
	validinds = cin.bodyId.isin(all_val)

	#find PNs that are in dataset but not in Greg's data
	allpninds = cin.name.str.contains("PN")
	unmatchedpns = cin[allpninds & (~validinds)]
	# print("\nPNs not matched to Greg's data:")
	# print(*unmatchedpns.name.values,sep='\n')

	Win_val = Win[:,validinds]
	cin_val = cin[validinds].copy()
	cin_val.reset_index(drop=True,inplace=True)

	cin_val["glom"] = ""
	cin_val["val"] = ""
	cin_val["pos_neg"] = ""
	cin_val["Type"] = ""
	cin_val["priority"] = ""

	#make everything lower case
	allpns.val = allpns.val.str.lower()

	#annotate cin_val with the valence, glomerulus, and type 
	for ii in allpns.index:
		cin_val.loc[cin_val.bodyId == allpns.skid[ii],"glom"] = allpns.glomerulus[ii]
		cin_val.loc[cin_val.bodyId == allpns.skid[ii],"Type"] = allpns.type[ii]
		cin_val.loc[cin_val.bodyId == allpns.skid[ii],"priority"] = allpns.priority[ii]
	for ii in visual.index:
		cin_val.loc[cin_val.bodyId == visual.bodyId[ii],"val"] = "visual"
		cin_val.loc[cin_val.bodyId == visual.bodyId[ii],"Type"] = "visual"
		cin_val.loc[cin_val.bodyId == visual.bodyId[ii],"priority"] = visual.priority[ii]


	#reorder by valence
	cin_val.sort_values(by=["priority","val"],inplace=True)
	Win_val = Win_val[:,cin_val.index]
	cin_val.reset_index(drop=True,inplace=True)


	valences =  pd.read_csv(datapath+"jefferis_glomerular_data_final/odour_scenes.csv", sep=',')


	for i in valences.index:
		glom = valences.glomeruli[i]
		for pn in range(len(cin_val.name)):
			if glom == cin_val.glom[pn]:
				cin_val.val[pn] = valences.odour_scene[i]
				cin_val.pos_neg[pn] = valences.valence[i]



	###

	# for val in cin_val.val.unique():
	#     if pd.isna(val):
	#         print("mean W for",val,":",np.mean(Win_val[:,cin_val.val.isna()]))
	#     else:
	#         print("mean W for",val,":",np.mean(Win_val[:,cin_val.val==val]))

	Win_vis = Win_val[:,cin_val.val == "visual"]
	frac_vis = np.sum(Win_vis,1)/np.sum(Win_val,1) #fraction of synaptic input from visual neurons
	frac_vis[np.isnan(frac_vis)] = 0

	in_deg = np.sum(Win_val > 0,1)

	in_deg_vis = in_deg[frac_vis > 0.8]
	in_deg_nonvis = in_deg[frac_vis < 0.2]

	total_in_deg = np.array([np.sum(x.relation == "upstream") for x in conns])
	total_in_deg_vis = total_in_deg[frac_vis > 0.8]
	total_in_deg_nonvis = total_in_deg[frac_vis < 0.2]

	total_in_syns = np.array([np.sum((x.relation == "upstream")*x.syncount) for x in conns])
	total_in_syns_vis = total_in_syns[frac_vis > 0.8]
	total_in_syns_nonvis = total_in_syns[frac_vis < 0.2]


	vals = cin_val.val.astype(str).values
	valswitch = np.append(0,1+np.where(vals[0:-1] != vals[1:])[0]) #indices of valence boundaries
	valnames = vals[valswitch]
	for ii in range(len(valnames)):
		if len(valnames[ii]) > 32:
			valnames[ii] = valnames[ii][:32]

	W = pd.DataFrame(np.copy(Win_val))

	thermo_hygro_bodyIds = [1639243580,
	1639234609,
	1727979406,
	5901222731,
	2065197353,
	2038617453,
	2069644133,
	1225100939,
	1788676257,
	5813069447,
	5813056072,
	5813063239,
	1974846717,
	1881401277,
	1943811736,
	1881059780,
	1943812176,
	1912777226,
	5812996748,
	1858901026,
	850717220,
	509626405,
	850708783,
	603785283,
	1975187675,
	1975187554,
	1975878958,
	5813040515,
	2069648663,
	1944502935,
	1755556097,
	1039237469,
	543010474,
	634759240,
	663432544,
	1789013008,
	663787020,
	664814903,
	5813077788,
	2100010400,
	1883443112]

	for i in range(len(cin_val)):
		if cin_val.bodyId[i] in thermo_hygro_bodyIds:
			cin_val.val[i]="thermo_hygro"

	return W, cin_val, c

# Return binary FAFB weight matrix as a dataframe
def get_FAFB():
	bin_mat = np.load('../Data/200514_pn_kc_bi_conn.npy') # Shape (1356, 113)
	pns = [2863104, 57349, 57353, 11544074, 16, 23569, 57361, 43539, 57365, 11524119, 192547, 57381, 36390, 57385, 23597, 24622, 37935, 400943, 775731, 67637, 11544121, 57402, 22594, 57410, 57414, 30791, 57418, 57422, 68697, 1775706, 23134, 56424, 41578, 27246, 33903, 22132, 35447, 37513, 32399, 24726, 57499, 27295, 771242, 30891, 57003, 57516, 51886, 45242, 24251, 40637, 49865, 1785034, 28876, 186573, 46800, 73937, 22744, 30434, 39139, 65762, 27884, 39668, 39682, 22277, 36108, 23829, 61221, 40749, 55085, 56623, 54072, 45882, 58686, 61773, 755022, 67408, 55125, 39254, 41308, 40306, 22906, 53631, 60799, 37250, 23432, 51080, 52106, 22422, 57241, 46493, 57246, 581536, 53671, 27048, 35246, 42927, 42421, 165303, 65465, 22976, 32214, 23512, 27611, 57307, 57311, 62434, 38885, 57319, 57323, 21999, 57333, 57337, 57341]
	types = pd.read_csv('../Data/180613-pn_subtypes.csv')

	# Map each skid to the appropriate glom
	pns_gloms = []
	for id in pns:
		if id in types['skids'].values:
			pns_gloms.append(types[types['skids']==id]['gloms'].item())

	FAFB = pd.DataFrame(bin_mat)
	FAFB.columns = pns_gloms
	FAFB = FAFB.groupby(lambda x:x, axis=1).sum() # Group like glomeruli

	return FAFB

# Shuffle matrix returning shuffled matrix mat_W with indegree of mat_X, keeping connection probs consistent
# shufmat
def fixed(mat_W, mat_X):
    M = mat_W.shape[0]
    N = mat_W.shape[1]

    indeg = np.sum(mat_X>0,1)
    
    cprobs = np.mean(mat_W>0,0)
    cprobs = cprobs / np.sum(cprobs)

    Wshuf = np.zeros([M,N])

    for mi in range(M):
        num_inputs = np.random.choice(indeg)
        inds = np.random.choice(N,num_inputs,p=cprobs, replace=False)
        Wshuf[mi,inds] = 1

    return Wshuf