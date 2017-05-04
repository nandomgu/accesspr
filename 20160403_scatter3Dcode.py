import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



h5PeakTime=df[df['strain']=='Hxt5']['FLAlignedPeakTime']
h5normpeak=df[df['strain']=='Hxt5']['FLPeak']/df[df['strain']=='Hxt5']['InitialFLperOD']
h5OD=df[df['strain']=='Hxt5']['FinalOD']


h2PeakTime=df[df['strain']=='Hxt2']['FLAlignedPeakTime']
h2normpeak=df[df['strain']=='Hxt2']['FLPeak']/df[df['strain']=='Hxt2']['InitialFLperOD']
h2OD=df[df['strain']=='Hxt2']['FinalOD']


h3PeakTime=df[df['strain']=='Hxt3']['FLAlignedPeakTime']
h3normpeak=df[df['strain']=='Hxt3']['FLPeak']/df[df['strain']=='Hxt3']['InitialFLperOD']
h3OD=df[df['strain']=='Hxt3']['FinalOD']


h4PeakTime=df[df['strain']=='Hxt4']['FLAlignedPeakTime']
h4normpeak=df[df['strain']=='Hxt4']['FLPeak']/df[df['strain']=='Hxt4']['InitialFLperOD']
h4OD=df[df['strain']=='Hxt4']['FinalOD']


fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
ax.scatter(h4OD, h4PeakTime, h4normpeak, c='blue', marker='o')
ax.scatter(h2OD, h2PeakTime, h2normpeak, c='green', marker='^')
ax.scatter(h3OD, h3PeakTime, color='red', marker='s')
ax.scatter(h3OD, h5PeakTime, color='purple', marker='+')


scatter(h4OD, h4PeakTime, c='blue', marker='o')
scatter(h2OD, h2PeakTime,  c='green', marker='^')
scatter(h3OD, h3PeakTime, color='red', marker='s')
scatter(h3OD, h5PeakTime, color='purple', marker='+')


ax.set_xlabel('Final OD')
ax.set_ylabel('Fluorescence peak time')
ax.set_zlabel('Max fluorescence')

plt.show()
