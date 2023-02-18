import time, subprocess
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import sys

now_dir = os.getcwd()
if not os.path.isdir("out"):
	os.mkdir("out")
notepad = subprocess.Popen(now_dir + '/program.exe')
time_max_sec = 20
waiting = 0
while 1:
	if notepad.poll() == None:
		pass
	elif notepad.poll() == 0:
		print ('Python: C++ worked good! C++ working time (in sec) was ' + str(waiting) + '\n')
		break
	else:
		print ('Python: C++ has a poblem! C++ returns exit code ' + str(notepad.poll()) + '\n')
		break
	time.sleep(1)
	waiting += 1
	if waiting <= time_max_sec:
		pass
	else:
		print ('Python: C++ timeout error! Time max (in sec) was ' + str(time_max_sec) + '\n')
		break


work_info = open("work_info.txt", "r")
lines_wi = work_info.readlines()
for name in lines_wi:
	name = name[:-1]

	plt.close('all')
	fig = plt.figure()
	# subplot
	#ax = fig.gca(projection='3d') #old version
	ax = fig.add_subplot(projection='3d')
	X = []
	Y = []
	Z = []

	try:  # numpy loadtxt
		with open(name) as f:
			for line in f:
				z, x, y = map(float, line.split())
				Z.append(z)
				X.append(x)
				Y.append(y)
		# matplotlib contour contourf		
		surf = ax.plot_trisurf(X, Y, Z, cmap=cm.coolwarm)
		ax.view_init(40, 220)
		
		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		 
		fig.colorbar(surf, shrink=0.5, aspect=5)
		
	#	plt.show()
		new_file_name = name[:-4] + '.png'
		plt.savefig(fname=new_file_name)
		print('Python: Everything is OK! Graph is ready here: ' + new_file_name)
	except FileNotFoundError:
		print('Python: No such file named!', repr(name))
work_info.close
