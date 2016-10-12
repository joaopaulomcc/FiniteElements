import numpy as np 
import math

################################################################
#Classes

class material:
	def __init__(self, name, density, young, poisson, *args):
		self.name = name
		self.density = density
		self.young = young
		self.poisson = poisson

class propertie:
	def __init__(self, name, kind, material, parameters, *args):
		self.name = name
		self.kind = kind
		self.material = material
		self.parameters = parameters

class load:
	def __init__(self, name, node, FX, FY, FZ, MX, MY, MZ):
		self.name = name
		self.node = node
		self.FX = FX
		self.FY = FY
		self.FZ = FZ
		self.MX = MX
		self.MY = MY
		self.MZ = MZ

class constrain:
	def __init__(self, name, node, X, Y, Z, RX, RY, RZ):
		self.name = name
		self.node = node
		self.X = X
		self.Y = Y
		self.Z = Z
		self.RX = RX
		self.RY = RY
		self.RZ = RZ

class job_data:
	def __init__(self, materials, properties, loads, constrains, elements, nodes):
		self.materials = materials
		self.properties = properties
		self.loads = loads
		self.constrains = constrains
		self.elements = elements
		self.nodes = nodes

##################################################################################
#Functions
def read_input_file(input_file):
	materials = {}
	properties = {}
	loads = {}
	constrains = {}
	elements = []
	nodes = []

	flag = None

	for line in input_file:
	
		line_s = line.split()
	
		if line[0] == "#":
			if line_s[0] == "#":
				flag = None
			elif line_s[0] == "#TITLE":
				flag = "title"
			elif line_s[0] == "#ANALYSIS":
				flag = "analysis"
			elif line_s[0] == "#MATERIALS":
				flag = "materials"
			elif line_s[0] == "#PROPERTIES":
				flag = "properties"
			elif line_s[0] == "#LOADS":
				flag = "loads"
			elif line_s[0] == "#CONSTRAINS":
				flag = "constrains"
			elif line_s[0] == "#ELEMENTS":
				flag = "elements"
			elif line_s[0] == "#NODES":
				flag = "nodes"
	
		else:
		
			if flag == "title":
				title = line
	
			elif flag == "analysis":
				analysis = line
		
			elif flag == "materials":			
				materials[line_s[0]] = material(line_s[0], float(line_s[1]), float(line_s[2]), float(line_s[3]))
	
			elif flag == "properties":
				properties[line_s[0]] = propertie(line_s[0], line_s[1], line_s[2], line_s[3:])
	
			elif flag == "loads":
				loads[line_s[0]] = load(line_s[0], int(line_s[1]), float(line_s[2]), float(line_s[3]), float(line_s[4]), float(line_s[5]), float(line_s[6]), float(line_s[7]) )
	
			elif flag == "constrains":
				constrains[line_s[0]] = constrain(line_s[0], int(line_s[1]), line_s[2], line_s[3], line_s[4], line_s[5], line_s[6], line_s[7])
	
			elif flag == "elements":
				elements.append(line_s)
	
			elif flag == "nodes":
				nodes.append(line_s)

	return job_data(materials, properties, loads, constrains, elements, nodes)

def BAR(element, job):
	
	element_propertie_name = element[1]
	element_propertie = job.properties[element_propertie_name]
	element_material_name = element_propertie.material
	element_material = job.materials[element_material_name]

	E = element_material.young
	A = float(element_propertie.parameters[0])

	#Problens for elements with different number of nodes!!!!
	Node1 = job.nodes[int(element[2]) - 1]
	Node2 = job.nodes[int(element[3]) - 1]
	
	X1 = float(Node1[1])
	Y1 = float(Node1[2])
	Z1 = float(Node1[3])

	X2 = float(Node2[1])
	Y2 = float(Node2[2])
	Z2 = float(Node2[3])

	L = math.sqrt((X2 - X1)**2 + (Y2 - Y1)**2)

	c = (X2 - X1)/L
	s = (Y2 - Y1)/L
	
	K = np.array([[1, -1], [-1, 1]])*(E*A/L)
	T = np.array([[c, s, 0, 0], [0, 0, c, s]])

	K_Local = np.dot(np.dot(np.transpose(T), K), T)

	return K_Local
##################################################################################
#Read Input File

input_file = open("input.inp", 'r')
job = read_input_file(input_file)

K_Global = np.zeros((len(job.nodes), len(job.nodes)))
DOF_Global = np.zeros((len(job.nodes), 1))
F_Global = np.zeros((len(job.nodes), 1))

for element in job.elements:
	if job.properties[element[1]].kind == "BAR":
		K_Local = BAR(element, job)

		for i in range(len(K_Local)):
			for j in range(len(K_Local)):
				#i_g = int(element[i + 1]) - 1
				#j_g = int(element[j + 1]) - 1

				print(i)
				print(element[i + 2])
				print(j)
				print(element[j + 2])

				#K_Global[i_g][j_g] += K_Local[i][j]
				#Relate nodes with Degrees of Freedon


	