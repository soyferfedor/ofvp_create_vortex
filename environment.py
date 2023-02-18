import time, subprocess
notepad = subprocess.Popen('/home/fedor/project/program.exe')
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
