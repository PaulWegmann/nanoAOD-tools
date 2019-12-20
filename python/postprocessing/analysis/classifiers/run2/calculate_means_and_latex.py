import numpy as np
import pandas as pd
import sys

triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
data = pd.read_csv(sys.argv[1])
col = data.columns.values
print(col)
data.drop(col[0], axis = 1, inplace = True)
data.drop("MSE", axis = 1, inplace = True)
data.drop("MAD", axis = 1, inplace = True)
data.drop("AMS", axis=1, inplace = True)

#data = data[np.array(["NTrees=500:NTrees" not in i for i in data.options]) & \
#	np.array(["MinNodeSize" not in i for i in data.options]) & np.array(["MaxDepth" not in i for i in data.options])]
#data.reset_index(drop=True, inplace=True)

final = pd.DataFrame(columns = np.append(data.columns.values, "AMSstd"))

for trigger in triggers:
    for option in np.unique(data.leftvar.values[trigger == data.trigger.values]):
        print(option)
        mean =  data[np.array([option == i for i in data.leftvar]) & np.array([trigger == i for i in data.trigger])].mean(axis = 0)
        std  =  data[np.array([option == i for i in data.leftvar]) & np.array([trigger == i for i in data.trigger])].std(axis = 0)
        mean["leftvar"] = option
        mean["trigger"] = trigger
        #print(pd.Series([round(std["AMS"],3)], index = ["AMSstd"]))
        #mean = mean.append(pd.Series([round(std["AMS"],3)], index = ["AMSstd"]))
        mean = mean.append(pd.Series([round(std["AMSrecon"],3)], index = ["AMSstd"]))

        #print(mean)
        final = final.append(mean, ignore_index = True)
        #print(final)
        
final.drop("trigger", axis = 1, inplace = True)

data = final


for i in range(len(data["leftvar"].values)):
    #data.set_value(i, "AMS", round(data["AMS"][i], 3))
    data.set_value(i, "AMSrecon", round(data["AMSrecon"][i], 3))
    data.set_value(i, "leftvar", data["leftvar"][i])#[34:])
    data.set_value(i, "ROCIntegral", round(data["ROCIntegral"][i]*0.01, 3))
    data.set_value(i, "sigevents", round(data["sigevents"][i], 0))
	#data.set_value(i, "MAD", round(data["MAD"][i], 3))
	#data.set_value(i, "MSE", round(data["MSE"][i], 3))
    data.set_value(i, "overtrain", round(data["overtrain"][i], 3))

data = data.reindex(["leftvar", "AMSrecon", "AMSstd", "ROCIntegral", "overtrain"], axis = 1)
data.columns = ["leftvar", "AMS", "AMSstd", "ROCIntegral", "overtrain"]


print(data)

pd.set_option('display.max_colwidth', -1)
data.to_latex("latex_"+sys.argv[1][:-4]+".txt", index=False)
data.to_csv("mean_"+sys.argv[1], index=False)

#print(data)


