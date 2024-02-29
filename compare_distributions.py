import os
import json
import numpy as np
import matplotlib.pyplot as plt


relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'






def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	data_directory = input_json["data_directory"]
	properties = 15 #number of columns to save for every Num

	sav_lognorm = data_directory+'data/lognorm_averageData.csv'
	sav_const = data_directory+'data/const_averageData.csv'
	figure_folder = data_directory+'data/figures/'


	temps = [3,10,30,100,300,1000]
	Nums = [30,100,300]
	savs = [sav_const,sav_lognorm]

	porositiesabcavg = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	porositiesabcstd = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_abc = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	ABC_numruns = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	porositiesKBMavg = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	porositiesKBMstd = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_KBM = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	KBM_numruns = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	contactsavg = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	contactsstd = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_ca = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	contacts_numruns = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	FD_dataavg = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	FD_datastd = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_FD = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	FD_numruns = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	angmomavg = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	angmomstd = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_angmom = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	angmom_numruns = np.full(shape=(len(savs),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)



	for s_i,sav in enumerate(savs):	
		data = np.loadtxt(sav,delimiter=' ',skiprows=0,dtype=np.float64)
		for i,N in enumerate(Nums):
			porositiesabcavg[s_i,i] = data[i*properties,:]
			yerr_abc[s_i,i] = data[i*properties+1,:]
			ABC_numruns[s_i,i] = data[i*properties+2,:]

			porositiesKBMavg[s_i,i] = data[i*properties+3,:]
			yerr_KBM[s_i,i] = data[i*properties+4,:]
			KBM_numruns[s_i,i] = data[i*properties+5,:]

			FD_dataavg[s_i,i] = data[i*properties+6,:]
			yerr_FD[s_i,i] = data[i*properties+7,:]
			FD_numruns[s_i,i] = data[i*properties+8,:]

			contactsavg[s_i,i] = data[i*properties+9,:]
			yerr_ca[s_i,i] = data[i*properties+10,:]
			contacts_numruns[s_i,i] = data[i*properties+11,:]

			angmomavg[s_i,i] = data[i*properties+12,:]
			yerr_angmom[s_i,i] = data[i*properties+13,:]
			angmom_numruns[s_i,i] = data[i*properties+14,:]


	
	ratio_porabc = np.zeros(shape=(len(Nums),len(temps)))
	porabc_uncert = np.zeros(shape=(len(Nums),len(temps)))

	ratio_porkbm = np.zeros(shape=(len(Nums),len(temps)))
	porkbm_uncert = np.zeros(shape=(len(Nums),len(temps)))

	ratio_numcon = np.zeros(shape=(len(Nums),len(temps)))
	numcon_uncert = np.zeros(shape=(len(Nums),len(temps)))

	ratio_fradim = np.zeros(shape=(len(Nums),len(temps)))
	fradim_uncert = np.zeros(shape=(len(Nums),len(temps)))
	
	for n in range(len(Nums)):
		for t in range(len(temps)):
			ratio_porabc[n,t] = porositiesabcavg[0,n,t]/porositiesabcavg[1,n,t]
			porabc_uncert[n,t] += ratio_porabc[n,t]*np.sqrt((yerr_abc[0,n,t]/porositiesabcavg[0,n,t])**2 + (yerr_abc[1,n,t]/porositiesabcavg[1,n,t])**2) #propagate uncertainty 
			
			ratio_porkbm[n,t] += porositiesKBMavg[0,n,t]/porositiesKBMavg[1,n,t]
			porkbm_uncert[n,t] += ratio_porkbm[n,t]*np.sqrt((yerr_KBM[0,n,t]/porositiesKBMavg[0,n,t])**2 + (yerr_KBM[1,n,t]/porositiesKBMavg[1,n,t])**2) #propagate uncertainty 
			
			ratio_numcon[n,t] += contactsavg[0,n,t]/contactsavg[1,n,t]
			numcon_uncert[n,t] += ratio_numcon[n,t]*np.sqrt((yerr_ca[0,n,t]/contactsavg[0,n,t])**2 + (yerr_ca[1,n,t]/contactsavg[1,n,t])**2) #propagate uncertainty 

			ratio_fradim[n,t] += FD_dataavg[0,n,t]/FD_dataavg[1,n,t]
			fradim_uncert[n,t] += ratio_fradim[n,t]*np.sqrt((yerr_FD[0,n,t]/FD_dataavg[0,n,t])**2 + (yerr_FD[1,n,t]/FD_dataavg[1,n,t])**2) #propagate uncertainty 

		print(f"Ratio of const / lognorm dist for N={Nums[n]}:")
		print(f"\tpor abc: {ratio_porabc[n,:]}")
		print(f"\tpor kbm: {ratio_porkbm[n,:]}")
		print(f"\tnum con: {ratio_numcon[n,:]}")
		print(f"\tfra dim: {ratio_fradim[n,:]}")


	styles = ['-','--','-.','--.']
	colors = ['g','b','r','orange','black']

	datas = [ratio_porabc,ratio_porkbm,ratio_fradim,ratio_numcon]
	methods = ["por abc","por kbm","fra dim","num con"]
	errs = [porabc_uncert,porkbm_uncert,fradim_uncert,numcon_uncert]
	for d_i,data in enumerate(datas):
		print(data[n]) 
		print(np.isnan(data[n]))
		if np.sum(np.isnan(data)) < data.size:
			fig,ax = plt.subplots()
			for n in range(len(Nums)):
				ax.errorbar(temps,data[n],yerr=errs[d_i][n],\
							label="N={}".format(Nums[n]),color=colors[d_i],linestyle=styles[n],zorder=5)
			ax.set_xlabel('Temp (K)')
			ax.set_title('Ratio of const/lognorm')
			ax.set_ylabel('Ratio')

			ax.axhline(1)
			ax.set_xscale('log')

			fig.legend()
			plt.savefig(f"{figure_folder}distComp_{methods[d_i].replace(' ','_')}.png")
			plt.show()
			plt.close("all")





if __name__ == '__main__':
	main()

