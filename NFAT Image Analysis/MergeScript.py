# -*- coding: utf-8 -*-
"""
Script to join data exported from Fiji_CellDataScript

"""

from PyQt5.QtWidgets import*
import pandas as pd
import numpy as np


dlg = QFileDialog()
dlg.setFileMode(QFileDialog.ExistingFiles)
paths= []

if dlg.exec_():
    paths = dlg.selectedFiles()
    
frames = []
for path in paths:
    frame = pd.read_csv(path, sep = ",")
    frame["mean fraction"] = frame["nucleus mean"]/frame["cytoplasma mean"]
    frame["sampleShort"] = " ".join(path.replace(".csv","").split("/")[-1].split("_"))
    frames.append(frame)

allDf = pd.concat(frames)
tempDict = dict([[name,group["mean fraction"]] for name,group in allDf.groupby("sampleShort")])
finalFrame = pd.DataFrame(tempDict)
finalFrame.to_csv("filename" + ".csv")
