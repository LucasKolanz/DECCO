import numpy as np

def S(sigma):
	return np.sum(1/sigma**2)

def Si(i,sigma):
	return np.sum(i/(sigma**2))

def Sii(x,y,sigma):
	return np.sum(x*y/(sigma**2))

def main():
	headers = np.loadtxt("averageData.csv",delimiter=',',dtype=str)[0]
	data = np.loadtxt("averageData.csv",delimiter=',',dtype=np.float64,skiprows=1)
	data = np.transpose(data)
	# print(headers)
	# print(headers.shape)
	# print(np.transpose(data))
	# print(np.transpose(data).shape)

	x = np.log(data[0])
	y_data = data[1:]

	print("Slope and slope uncertainties")
	print('===========================================')
	for i in range(0,len(y_data),2):
		y = y_data[i]
		sigma = y_data[i+1]

		delta = S(sigma)*Sii(x,x,sigma)-(Si(x,sigma))**2
		slope = (S(sigma)*Sii(x,y,sigma)-Si(x,sigma)*Si(y,sigma))/delta

		slope_sigma_sq = S(sigma)/delta

		print(headers[i+1])
		print("\t    slope: {}".format(slope))
		print("slope uncertainty: {}".format(np.sqrt(slope_sigma_sq)))
		if i < len(y_data)-2:
			print()
	print('===========================================')

		

if __name__ == '__main__':
	main()