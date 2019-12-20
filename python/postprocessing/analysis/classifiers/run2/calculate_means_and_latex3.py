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
data.drop("AMS", axis = 1, inplace = True)


#data = data[np.array(["NTrees=500:NTrees" not in i for i in data.options]) & \
#	np.array(["MinNodeSize" not in i for i in data.options]) & np.array(["MaxDepth" not in i for i in data.options])]
#data.reset_index(drop=True, inplace=True)

final = pd.DataFrame(columns = np.append(data.columns.values, "AMSstd"))

for trigger in triggers:
    for option in np.unique(data.options.values[trigger == data.trigger.values]):
        print(option)
        mean =  data[np.array([option == i for i in data.options]) & np.array([trigger == i for i in data.trigger])].mean(axis = 0)
        std  =  data[np.array([option == i for i in data.options]) & np.array([trigger == i for i in data.trigger])].std(axis = 0)
        mean["options"] = option
        mean["trigger"] = trigger
        #print(pd.Series([round(std["AMS"],3)], index = ["AMSstd"]))
        # mean = mean.append(pd.Series([round(std["AMS"],3)], index = ["AMSstd"]))
        mean = mean.append(pd.Series([round(std["AMSrecon"],3)], index = ["AMSreconstd"]))
        mean = mean.append(pd.Series([round(std["sigevents"],3)], index = ["sigeventsstd"]))

        #print(mean)
        final = final.append(mean, ignore_index = True)
        #print(final)
        
final.drop("trigger", axis = 1, inplace = True)

data = final


for i in range(len(data["options"].values)):
    data.set_value(i, "AMSrecon", round(data["AMSrecon"][i], 3))
    data.set_value(i, "sigevents", round(data["sigevents"][i], 0))
    data.set_value(i, "sigeventsstd", round(data["sigeventsstd"][i], 0))


    data.set_value(i, "options", data["options"][i])#[34:])
    data.set_value(i, "ROCIntegral", round(data["ROCIntegral"][i]*0.01, 3))
    data.set_value(i, "overtrain", round(data["overtrain"][i], 3))

data = data.reindex(["options", "AMSrecon", "AMSreconstd", "sigevents", "sigeventsstd", "bckevents", "ROCIntegral", "overtrain"], axis = 1)
data.drop("bckevents", axis =1, inplace = True)
data.columns = ["Option", "AMS", "AMS std", "Nsignal", "Nsignalstd", "ROCIntegral", "Overtrain"]

print(data)

pd.set_option('display.max_colwidth', -1)
data.to_latex("latex.txt", index=False)

data.to_csv("mean_"+sys.argv[1], index = False)

