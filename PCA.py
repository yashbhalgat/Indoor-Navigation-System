from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
from pylab import plot,subplot,axis,stem,show,figure

def princomp(A):
 """ performs principal components analysis 
     (PCA) on the n-by-p data matrix A
     Rows of A correspond to observations, columns to variables. 

 Returns :  
  coeff :
    is a p-by-p matrix, each column containing coefficients 
    for one principal component.
  score : 
    the principal component scores; that is, the representation 
    of A in the principal component space. Rows of SCORE 
    correspond to observations, columns to components.
  latent : 
    a vector containing the eigenvalues 
    of the covariance matrix of A.
 """
 # computing eigenvalues and eigenvectors of covariance matrix
 M = (A-mean(A.T,axis=1)).T # subtract the mean (along columns)
 [latent,coeff] = linalg.eig(cov(M)) # attention:not always sorted
 score = dot(coeff.T,M) # projection of the data in the new space
 return coeff,score,latent

if __name__=="__main__":
	myfile = open('walking_mag_static23.txt','r')
	data = myfile.read().split('\n')
	myfile.close()

	for i in range(len(data)-1):
	 data[i] = data[i].split(',')
	 data[i][0] = float(data[i][0])
	 data[i][1] = float(data[i][1])
	 data[i][2] = float(data[i][2])

	data_array = [[],[],[]]
	for i in range(len(data)-1):
	 data_array[0].append(data[i][0])
	 data_array[1].append(data[i][1])
	 data_array[2].append(data[i][2])

	mean_data_array = [range(len(data)-1),range(len(data)-1),range(len(data)-1)]
	for i in range(len(data)-1):
	 mean_data_array[0][i] = data_array[0][i]-sum(data_array[0]) / float(len(data_array[0]))
	 mean_data_array[1][i] = data_array[1][i]-sum(data_array[1]) / float(len(data_array[1]))
	 mean_data_array[2][i] = data_array[2][i]-sum(data_array[2]) / float(len(data_array[2]))
	#print data_array[0]
	#print data_array[1]
	#print data_array[2]

	data_array1 = [data_array[0],data_array[1]]


	A = array(data_array1)
	coeff, score, latent = princomp(A.T)
	print coeff
	print score
	print latent

	from matplotlib import pyplot
	import pylab
	from mpl_toolkits.mplot3d import Axes3D
	import random


	fig = pylab.figure()
	ax = Axes3D(fig)

	sequence_containing_x_vals = mean_data_array[0]
	sequence_containing_y_vals = mean_data_array[1]
	sequence_containing_z_vals = mean_data_array[2]

	ax.scatter(sequence_containing_x_vals, sequence_containing_y_vals)#,sequence_containing_z_vals)
	pylab.show()


	#figure()
	#subplot(121)
	# every eigenvector describe the direction
	# of a principal component.
	#m = mean(A,axis=1)
	#plot([0, -coeff[0,0]*2]+m[0], [0, -coeff[0,1]*2]+m[1],'--k')
	#plot([0, coeff[1,0]*2]+m[0], [0, coeff[1,1]*2]+m[1],'--k')
	#plot(A[0,:],A[1,:],'ob') # the data
	#axis('equal')
	#subplot(122)
	# new data
	#plot(score[0,:],score[1,:],'*g')
	#axis('equal')
	#show()
