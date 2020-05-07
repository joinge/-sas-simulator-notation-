import numpy as np
import traceback

class ObjLoader:
   """Very simple object for loading Wavefront (obj) files.
   
   Only vertices and normals are preserved, the rest gets discarded.
   """
   def __init__(self):
      self.vert_coords = []
      self.text_coords = []
      self.norm_coords = []

      self.vertex_index = []
      self.texture_index = []
      self.normal_index = []

      self.vertices = []
      self.textures = []
      self.normals  = []

      self.materials = None

   def loadMaterials(self, file):
      name = None
      color = None
      alpha = None

      for line in open(file, 'r'):
            if line.startswith('#'):
               continue
            values = line.split()
            if not values:
               continue

            if values[0] == 'newmtl':  # New material
               name  = values[1:]
               color = None
               alpha = None

            if values[0] == 'Kd':      # Color
               color = values[1:3]
            if values[0] == 'd':       # Alpha
               alpha = values[1]


   def loadModel(self, file):

      try:
         current_material = None
         for line in open(file, 'r'):
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue

            if values[0] == 'v': # Vertex
               self.vert_coords.append(values[1:4])
            if values[0] == 'vt': # UV
               self.text_coords.append(values[1:3])
            if values[0] == 'vn': # Normal
               self.norm_coords.append(values[1:4])

            if values[0] == 'usemtl':
               current_material = values[1]

            if values[0] == 'f':
               vertex_index = []
               uv_index = []
               normal_index = []
               for v in values[1:4]:
                  w = v.split('/')
                  vertex_index.append(int(w[0])-1)
                  if w[1] != '':
                     uv_index.append(int(w[1])-1)
                  normal_index.append(int(w[2])-1)
               self.vertex_index.append(vertex_index)
               self.texture_index.append(uv_index)
               self.normal_index.append(normal_index)

         self.vertex_index  = [y for x in self.vertex_index  for y in x]
         self.texture_index = [y for x in self.texture_index for y in x]
         self.normal_index  = [y for x in self.normal_index  for y in x]

         for i in self.vertex_index:
            self.vertices.extend(self.vert_coords[i])

         for i in self.texture_index:
            self.textures.extend(self.text_coords[i])

         for i in self.normal_index:
            self.normals.extend(self.norm_coords[i])

         self.model = []
         self.model.extend(self.vertices)
         self.model.extend(self.normals)

         self.vertices = np.array(self.vertices, dtype='float32')
         self.normals  = np.array(self.normals,  dtype='float32')

      except Exception as e:
         traceback.print_exc()
         print(e)

