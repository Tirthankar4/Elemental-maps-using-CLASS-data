import pandas as pd
import csv 

dir_path = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios"

R2_thres=0.8
output = pd.read_csv(f"{dir_path}/output.csv")
N=output.shape[0]
elementList="Fe Ti Ca Si Al Mg Na O".split(" ")
i=0
flareNum=68
fileNames = [f"{dir_path}/results/{element}_{flareNum}.csv" for element in elementList]
writers = [csv.writer(open(_file, 'w')) for _file in fileNames]

header = "BORE_LAT BORE_LON Integration_time amplitude_ratio_vs_Si".split(" ")
[writer_.writerow(header) for writer_ in writers]

while i < N:
    elem_bool = [False]*8
    temp = i
    while i < temp+4:
        row = output.iloc[i].tolist()
        if row[2]==8 and i!=temp:
            break
        for elemId in range(len(elementList)):
            elemDetected = row[3+6*elemId+1]
            elemSigma = row[3+6*elemId+3]
            elemR2 = row[3+6*elemId+4]
            SiR2 = row[3+6*3+4]
            elemRatio = row[3+6*elemId+5]
            if (elemDetected and (0.03< elemSigma <0.1) and (elemR2>R2_thres) and (SiR2 > R2_thres) and elem_bool[elemId]==False):                
                writers[elemId].writerow([row[0], row[1], row[2], elemRatio])
                elem_bool[elemId]=True

        i+=1
        


