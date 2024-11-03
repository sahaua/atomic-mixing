## This program calculates atomic mixing (apart from the E_Dn in denominator),
## that is, R-square from the trjaectory_out.xyz files output from the 
## TurboGAP simulations.
## -- Uttiyoarnab Saha --##
##-------------------------------------------------------------------------------

import numpy, os, datetime

day_execution = datetime.date.today()
time_execution = datetime.datetime.now().strftime('%H:%M:%S')
knowcwd = os.getcwd()

## get some inputs from the user

E_pka = float(input('PKA energy (eV): '))
num_atoms = int(input('Number of atoms in box: '))
side_length_x = float(input('Side length along x: '))
side_length_y = float(input('Side length along y: '))
side_length_z = float(input('Side length along z: '))

if (E_pka <= 0 or num_atoms <=0 or side_length_x <= 0 or side_length_y <= 0 or side_length_z <= 0):
	raise Exception('None of Energy, # atoms, length of a side can be 0 or negative!')

## this should be N_MDsteps/dump_frequency + 1
num_frames = int(input('Total number of frame dumps including 0th dump: '))

if (num_frames < 2):
	raise Exception('Including 0th step, there should be atleast 2 frame dumps in output trajectory file!')

## total number of different random trajectory cases simulated 
num_cases = int(input('Total number of data set (count starting from 1): '))

if (num_cases < 1):
	raise Exception('There should be atleast 1 dataset!')

## folder path of trajectories relative to the working directory
folder_path_trajectories = input('Enter path to trajectories relative to working directory (name should be the common part of the folders in case of more than 1 case): ')

## name of model used in simulation
model_name = input('Name of case, SRIM-ES, EPH, etc.: ')

## folder path to save output files relative to the current working directory
folder_path_outputs = input('Enter path to output folder relative to the working directory: ')

## find the n0
n_atom_density = num_atoms/(side_length_x*side_length_y*side_length_z)

## arrays to hold poitions
x0 = numpy.zeros(num_atoms)
y0 = numpy.zeros(num_atoms)
z0 = numpy.zeros(num_atoms)

xt = numpy.zeros(num_atoms)
yt = numpy.zeros(num_atoms)
zt = numpy.zeros(num_atoms)

## array to store the R-square values  
R_square = numpy.zeros((num_cases,num_frames))

## array to store the times
times = numpy.zeros((num_cases,num_frames))

for i in range(num_cases):
	print('case ', i+1)

	## read the trajectory files
	## to calculate R-square

	if (num_cases > 1):
		os.chdir(folder_path_trajectories+str(i+1))
	if (num_cases == 1):
		os.chdir(folder_path_trajectories)

	## Find out positions
	time_prev = 0
	ifile = open ('trajectory_out.xyz', 'r')
	flag = 0
	for j in range(num_frames):
		ifile.readline()
		## get the value of time from the second line
		line = ifile.readline()
		data = line.split()
		for k in range(len(data)):
			if (data[k][0:5] == 'time='):
				times[i][j] = float(data[k][5:])/1000.0
		## to handle (skip) repetition in consecutive trajectory snapshots
		## also, TurboGAP starts printing from time value = 0 in trajectory outputs for restarted runs
		if (j > 0):
			if (times[i][j] == 0.0 or (times[i][j] > 0.0 and times[i][j] == time_prev)):
				print('..repeated run <--')
				## since TurboGAP starts printing from time value = 0 in trajectory outputs for restarted runs
				if (times[i][j] == 0.0):
					jbase = j - 1
					flag = 1

				for k in range(num_atoms):
					ifile.readline()
				ifile.readline()
				## get the value of time from the second line
				line = ifile.readline()
				data = line.split()
				for k in range(len(data)):
					if (data[k][0:5] == 'time='):
						times[i][j] = float(data[k][5:])/1000.0

		time_prev = times[i][j]
		if (flag == 1):
			times[i][j] += times[i][jbase] 
		print(times[i][j])

		sum_sq = 0
		if (j == 0):
			for k in range(num_atoms):
				line = ifile.readline()
				x0[k] = float(line.split()[1])
				y0[k] = float(line.split()[2])
				z0[k] = float(line.split()[3])

		if (j > 0):
			for k in range(num_atoms):
				line = ifile.readline()
				xt[k] = float(line.split()[1])
				yt[k] = float(line.split()[2])
				zt[k] = float(line.split()[3])

			for k in range(num_atoms):
				dx = xt[k] - x0[k]
				dy = yt[k] - y0[k]
				dz = zt[k] - z0[k]
				if (dx > side_length_x * 0.5):
					dx = dx - side_length_x
				if (dx <= -side_length_x * 0.5):
					dx = dx + side_length_x
				if (dy > side_length_y * 0.5):
					dy = dy - side_length_y
				if (dy <= -side_length_y * 0.5):
					dy = dy + side_length_y
				if (dz > side_length_z * 0.5):
					dz = dz - side_length_z
				if (dz <= -side_length_z * 0.5):
					dz = dz + side_length_z

				sum_sq = sum_sq + dx*dx + dy*dy + dz*dz

		R_square[i][j] = sum_sq / (6*n_atom_density)

	ifile.close()

os.chdir('../')

R_square_av = numpy.mean(R_square, axis = 0)
R_square_stde = numpy.std(R_square, axis = 0)
times_av = numpy.mean(times, axis = 0)

for j in range(num_frames):
	R_square_stde[j] /= numpy.sqrt(num_cases)

## save the calculated ouputs in the specified directory

os.chdir(knowcwd)
os.chdir(folder_path_outputs)

## output file containing calculated values

ofile = open('R-square-upon-6n0_'+str(E_pka)+'_'+model_name, 'w')
print(f'{E_pka} eV, model {model_name}, time (ps)/mean(Rsq_6n0)/stde(Rsq_6n0)', file = ofile)
for j in range(num_frames):
	print(times_av[j], R_square_av[j], R_square_stde[j], file = ofile)
print(file = ofile)
print('Value of 6n0 = ', 6*n_atom_density, file = ofile)
print(file = ofile)
print('average R-square-upon-6n0 = ', numpy.mean(R_square_av[1:]), file = ofile)
print('standard_error R-square-upon-6n0 = ', numpy.std(R_square_av[1:])/numpy.sqrt(len(R_square_av[1:])), file = ofile)

ofile.close()

## output file containing given input data

ofile = open('log-input-'+str(E_pka)+'eV-'+model_name+'.Rsq_6n0', 'w')

print('LOG file -- calcR-square-upon-6n0.py --', file = ofile)
print(day_execution, time_execution, file = ofile)
print(file = ofile)

print('Following inputs were provided: ', file = ofile)
print('-----------------------------------------------------------', file = ofile)
print(file = ofile)
print('PKA energy (eV): ', E_pka, file = ofile)
print(file = ofile)
print('Number of atoms in box: ', num_atoms, file = ofile)
print(file = ofile)
print('Side length along x: ', side_length_x, file = ofile)
print(file = ofile)
print('Side length along y: ', side_length_y, file = ofile)
print(file = ofile)
print('Side length along z: ', side_length_z, file = ofile)
print(file = ofile)

print('Total number of frame dumps including 0th dump: ', num_frames, file = ofile)
print(file = ofile)
print('Total number of data set (count starting from 1): ', num_cases, file = ofile)
print(file = ofile)
print('Path to trajectories (name should be the common part of the folders in case of more than 1 case): ', \
folder_path_trajectories, file = ofile)

print(file = ofile)
print('Model (SRIM-ES, EPH, etc.): ', model_name, file = ofile)
print(file = ofile)
print('Path to ouput folder: ', folder_path_outputs, file = ofile)

ofile.close()
