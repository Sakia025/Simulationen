import os
import math
import yaml

### define, how many simulations will run to evaulate the uncertainty and mean

with open('Anzahl_Simulationen.yaml') as file:
	 Anzahl_Simulationen = yaml.load(file, Loader=yaml.FullLoader)['Anzahl_Simulationen']
nummerierung =  range(Anzahl_Simulationen)

### geant4 singularity container, für interactive Maschine
#aktuelles_Geant4_image = '/ceph/Singularity/Geant4-10.6.simg'

# all paths will be relative to this one, snakemake will put the output here
workdir: '/Users/smller/Simulationen/Snakefile'
#interaktive Maschine
#localrules: targets, create_run_RS, create_run_COB, build_RS, build_COB, TDK, Oberflaeche_RS, zusammenfassung,
MyConfigFile = 'Parameter_allgemein.yaml'
configfile: MyConfigFile

### define the scr and include files of the simulations (rotational symmetry)
src_dateien_RS = 'B1ActionInitialization.cc B1DetectorConstruction.cc B1PrimaryGeneratorAction.cc BrachyPhysicsList.cc'.split()
include_dateien_RS = 'B1ActionInitialization.hh B1DetectorConstruction.hh B1PrimaryGeneratorAction.hh BrachyPhysicsList.hh'.split()



### output der simulationen
rule targets:
	input:
		input1 = expand('B1-build/TDK.pdf'),
		input2 = expand('B1-build/TDK-Fitparameter.txt'),
		input3 = expand('B1-build/TDK-Werte.txt'),
		input4 = expand('B1-build/Oberflaeche.pdf'),
		input5 = expand('B1-build/Oberflaeche-Werte.txt'),
		input6 = ('Ergebnisse/TDK-Werte.txt', 'Ergebnisse/TDK-Fitparameter.txt', 'Ergebnisse/Oberflaeche-Werte.txt')



### erstellt die run-files mit hilfe des Pythonskriptes.
rule create_run_RS:
	input: 'py_run_dateien_RS.ipynb', MyConfigFile
	output:'B1/run{zahl}.mac', 'B1-build/run{zahl}_.mac'
	params:
		run_n = Anzahl_Simulationen,
		Teilchenzahl = config['Teilchenzahl'],
		winkelschritt = config['Winkelschritt'],
	notebook: 'py_run_dateien_RS.ipynb'

rule build_RS:
	input:
		input1 = ('B1/exampleB1.cc','B1/exampleB1.in','B1/exampleB1.out'),
		input2 = 'B1/CMakeLists.txt',
		input3 = (src_dateien_RS, include_dateien_RS),
	output:
		output1 = ('B1-build/exampleB1','B1-build/exampleB1.in','B1-build/exampleB1.out'),
		output2 = ('B1-build/cmake_install.cmake', 'B1-build/CMakeCache.txt', 'B1-build/init_vis.mac', 'B1-build/Makefile', 'B1-build/vis.mac'),
	#wildcard_constraints:
	#	applikator = separator.join(rs_applikatoren)
	singularity: aktuelles_Geant4_image
	threads: 1
	shell: 'mkdir -p {wildcards.applikator}-build && cd {wildcards.applikator}-build && cmake -DGeant4_DIR=/opt/geant4/lib/Geant4-10.6.0 ../{wildcards.applikator} && make '

### submittiert die simulationen an das cluster.
rule simulate_RS:
	input:
		input1 = rules.build_RS.output,
		input2 = 'B1-build/run{zahl}.mac',
	output:
		output1 = expand('{B1}-build/TDK_{{zahl}}.txt')
	#wildcard_constraints:
	#	applikator = separator.join(rs_applikatoren),
	threads: 1
	singularity: aktuelles_Geant4_image
	shell: 'cd {wildcards.applikator}-build && ./exampleB1 run{wildcards.zahl}.mac'


rule TDK:	# die TDK kann in jeder Simulation abgerufen werden
	input:
		input1 = expand('{B1}-build/TDK_{zahl}_0_0.txt', zahl=nummerierung),
		input2 = 'py_TDK.ipynb',
		input3 = ('BEBIG/TDK-Werte.txt', 'BEBIG/TDK-z.txt','BEBIG/TDK-US.txt'),
		input4 = 'Anzahl_Simulationen.yaml'
	output:
		output1 = ('B1-build/TDK.pdf','B1-build/TDK-Werte.txt', 'B1-build/TDK-Fitparameter.txt'),
		output2 = ('B1-build/TDK-BEBIG.pdf','B1-build/TDK-BEBIG.txt')
	params:
		run_n = Anzahl_Simulationen
	notebook: 'py_TDK.ipynb'



rule Oberflaeche_RS:	# bezieht sich nur auf die Rotationssymmetrischen
	input:
		input1 = expand('{B1}-build/TDK_{zahl}.txt', zahl=nummerierung),
		input2 = MyConfigFile,
		input3 = 'py_Oberflaeche.ipynb',
		input4 = 'Anzahl_Simulationen.yaml'
	output:
		output1 = ('B1-build/Oberflaeche-Werte.txt','B1-build/Oberflaeche.pdf')
	params:
		run_n = Anzahl_Simulationen,
		winkelschritt = config['Winkelschritt']
	#wildcard_constraints:
	#	applikator = separator.join(rs_applikatoren)
	notebook:	'py_Oberflaeche.ipynb'



rule zusammenfassung:
	input:
		input1 = expand('B1-build/TDK-Fitparameter.txt'),
		input2 = expand('B1-build/TDK-Werte.txt'),
		input3 = expand('B1-build/Oberflaeche-Werte.txt'),
		input4 = 'py_Zusammenfassung.ipynb'
	output:
		output1 = ('Ergebnisse/TDK-Werte.txt', 'Ergebnisse/TDK-Fitparameter.txt', 'Ergebnisse/Oberflaeche-Werte.txt')
	notebook: 'py_Zusammenfassung.ipynb'
