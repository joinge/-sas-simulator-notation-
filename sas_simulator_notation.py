#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@file
@brief Sample code for side-looking sonar simulation with OpenGL.   

Released as complementary source code to the journal article titled
   "Real-time Synthetic Aperture Sonar Simulation for Automatic Target Recognition: Notation and Application"
   (Buskenes J.I, Midelfart, H.)
   submitted to IEEE Journal of Oceanic Engineering, 2020.

No license restrictions applies. Use, alter and/or redistribute as you wish.
     Author:      Jo Inge Buskenes
     Affiliation: University of Oslo / The Norwegian Defense Research Establishment (FFI)
     Date:        May 2020
     
Credits/acknowledgements:
   Adopted from a template for pyopengl found at:
   https://github.com/lukecampbell/pyopengl-examples

Tested on various Ubuntu installations, with these package requirements:
   For Python3: python3-numpy python3-pyqt5 python3-pyqt5.qtopengl python3-opengl
   For Python2: python-numpy python-pyqt5 python-pyqt5.qtopengl python-opengl

GL version must be at least 4.30, force MESA to use it by setting these environment variables:  
   export MESA_GLSL_VERSION_OVERRIDE=430 
   export MESA_GL_VERSION_OVERRIDE=4.3
"""

import sys,time
from numpy import abs,append,arctan,array,asarray,block,bmat,column_stack,cos, \
                  diagflat,diff,dot,eye,ones,pi,sign,sin,sqrt,squeeze,vstack,zeros
from numpy.linalg import multi_dot

from PyQt5 import QtCore, QtWidgets
# from PyQt5.QtOpenGL import QGLWidget 

import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo

from ObjLoader import ObjLoader

# from SimCL import SimCL
# import pyopencl as cl


def computeGeometry(type, psi, theta, phi):
      # Reference frames:
      # V       - Some vertex in the loaded model
      # M       - Model frame
      # I       - Image frame (on the seafloor)
      # S       - Sonar frame (local body)
      # C       - Cartesian camera (OpenGL projection plane)
      # R       - Radial    camera (OpenGL projection surface, curved)
      # N       - Normalized device coordinates (NDC)
      #
      # Syntax:
      # p_AB_A  - 3x1 Position vector from frame A to frame B, decomposed in A
      # L_AB    - 3x3 Linear transformation matrix (scale, rotation, shear, reflection)
      # T_AB_A  - 4x4 Projective transformation matrix (adds translation and optionally projection) 

      def mkHomo(p):
         """Append 1 to one or more position vectors to make them homogeneous"""
         newshape = list(p.shape)
         newshape[0] += 1
         p_new = ones(newshape)
         p_new[:-1] = p
         return p_new
         
      ext_I_I  = array([4, 4, 0])  # The extent of the image (I)

      # Define translations between the reference frames
      p_IM_I = array([0, 0, 0])  # Image to model
      p_SI_S = array([0, 10,4])  # Sonar to image
      p_CS_C = array([0, 0, 0])  # Sonar to viewing plane
      
      #############
      # Rotations #
      #############
      
      # Our models load with Z-up, Y-forward (default Blender view).
      # Rotation by 180 about Y-axis aligns it with our I-frame,
      # so this is the first rotation we construct R_IM from.
      if type=='model':
         # For show, we vary the model's
         R_IM = multi_dot([
            Rx(psi),                # roll,
            Ry(theta),                # pitch, and
            Rz(phi),                # yaw.
            Ry(pi)
            # R([0,1,0],pi)
            ])
      elif type=='image_frame':
         # The reference frame should just be still
         R_IM = eye(3) #Ry(pi)
         
      # The S- and I-frames have same orientation.
      R_SI = eye(3,3)

      # Rotate about x-axis until z-axis and p_SI are aligned. Two alternatives:
      # 1. Compute the rotation angle and use some axis-angle formula
      R_CS = \
         Rx(arctan(p_SI_S[1]/p_SI_S[2]))

      # 2. Construct R_CS directly from the normalized version of p_SI_S 
      # n_SI_S = p_SI_S / sum(p_SI_S**2)**0.5
      # R_CS = array([
      #    [1,  0,          0],
      #    [0,  n_SI_S[2], -n_SI_S[1]],
      #    [0,  n_SI_S[1],  n_SI_S[2]]
      # ])
      
      # Define scaling
      s_IM = 1

      # Represent it all using homogeneous coordinates 
      T_IM_I = vstack([column_stack(
         [  R_IM,  p_IM_I ]),
         [0, 0, 0, 1/s_IM ]])
      
      T_SI_S = vstack([column_stack(
         [  R_SI,  p_SI_S ]),
         [0, 0, 0, 1      ]])
      
      T_CS_C = vstack([column_stack(
         [  R_CS,  p_CS_C ]),
         [0, 0, 0, 1      ]])

      # The conversion from cartesian (C) to radial (R) is a non-linear function 
      def T_RC(p_C):
         p_R = zeros(p_C.shape)
         p_R[0] = p_C[0]
         p_R[1] = arctan(p_C[1]/p_C[2])
         p_R[2] = sign(p_C[2])*sqrt(p_C[1]**2 + p_C[2]**2)
         return p_R
      
      # To compute the projection matrices, we first need the image frame edges
      p_IE_I = array([
         ext_I_I[0]/2, ext_I_I[1]/2, 0
         ])
      
      # Then represent it in C
      p_CE_C = multi_dot([T_CS_C, T_SI_S, mkHomo(p_IE_I)])[:-1]
      p_CI_C = multi_dot([T_CS_C, mkHomo(p_SI_S)])[:-1]

      # First crate T_NC: the orthographic projection matrix
      p_IE_C = p_CE_C-p_CI_C
      p_IE_C *= array([1.5,3,3])    # Space around image frame, for demonstration only)
      T_NC_N = vstack([column_stack(
         [diagflat(1/p_IE_C), -p_CI_C/p_IE_C]),
         [    0, 0, 0,      1           ]])
      # Flip x-axis for OpenGL complicance
      T_NC_N[0,0] *= -1
      
      # Then create T_NR: the radial projection matrix
      p_RE_R = T_RC(p_CE_C)
      p_RI_R = T_RC(p_CI_C)
      p_IE_R = p_RE_R-p_RI_R
      p_IE_R *= array([1.5,1,3]) # Space around image frame, for demonstration only
      p_IE_R[1] = 23*pi/180      # 23 degree opening angle
      T_NR_N = vstack([column_stack(
         [diagflat(1/p_IE_R),  -p_RI_R/p_IE_R]),
         [0, 0, 0,    1  ]])
      # Flip x-axis for OpenGL complicance
      T_NR_N[0,0] *= -1                 

      return T_IM_I, T_SI_S, T_CS_C, T_NC_N, T_NR_N



# Window creation function.
def createWindow(window_class):
   """Create a Qt window in Python, or interactively in IPython with Qt GUI
   event loop integration:
      # in ~/.ipython/ipython_config.py
      c.TerminalIPythonApp.gui = 'qt'
      c.TerminalIPythonApp.pylab = 'qt'
   See also:
      http://ipython.org/ipython-doc/dev/interactive/qtconsole.html#qt-and-the-qtconsole
   """
   app_created = False
   app = QtCore.QCoreApplication.instance()
   if app is None:
      app = QtWidgets.QApplication(sys.argv)
      app_created = True
   app.references = set()
   window = window_class()
   app.references.add(window)
   window.show()
   if app_created:
      app.exec_()
   return window


def compileShader(source, shader_type):
   shader = gl.glCreateShader(shader_type)
   gl.glShaderSource(shader, source)
   gl.glCompileShader(shader)
   # check compilation error
   result = gl.glGetShaderiv(shader, gl.GL_COMPILE_STATUS)
   if not(result):
      raise RuntimeError(gl.glGetShaderInfoLog(shader))
   return shader

def compileVertexShader(source):
   """Compile a vertex shader from source."""
   return compileShader(source, gl.GL_VERTEX_SHADER)

def compileFragmentShader(source):
   """Compile a fragment shader from source."""
   return compileShader(source, gl.GL_FRAGMENT_SHADER)

def compileComputeShader(source):
   """Compile a compute shader from source."""
   return compileShader(source, gl.GL_COMPUTE_SHADER)

def createProgram(*shaders):
   # if len(shaders) > 1:
   program = gl.glCreateProgram()
   for shader in shaders:
      gl.glAttachShader(program, shader)
   gl.glLinkProgram(program)

   # check linking error
   result = gl.glGetProgramiv(program, gl.GL_LINK_STATUS)

   if not(result):
      raise RuntimeError(gl.glGetProgramInfoLog(program))
   return program

def Rx(psi):
   return array([
   [1,            0,             0           ],
   [0,            cos(psi),      -sin(psi)   ],
   [0,            sin(psi),      cos(psi)    ],
   ])
   
def Ry(theta):
   return array([
   [cos(theta),   0,             sin(theta)  ],
   [0,            1,             0           ],
   [-sin(theta),  0,             cos(theta)  ],
   ])
   
def Rz(phi):
   return array([
   [cos(phi),     -sin(phi),     0           ],
   [sin(phi),     cos(phi),      0           ],
   [0,            0,             1           ],
   ])
   
def S(sx,sy,sz):
   return array([
   [sx,       0,        0  ],
   [0,        sy,       0  ],
   [0,        0,        sz ],
   ])

def R(ur,theta):
   I = eye(3)
   X = array([
      [ 0,      -ur[2],   ur[1] ],
      [ ur[2],   0,      -ur[0] ],
      [-ur[1],   ur[0],   0     ]])
   return I + X*sin(theta) + dot(X,X)*(1-cos(theta))

# Vertex shader
VS = """
#version 430
// Attribute variable that contains coordinates of the vertices.
layout(location = 0) in vec3 p_MV_M;
layout(location = 1) in vec3 n_MV_M_input;

out vec3 p_SV_S;
out vec3 n_MV_S;

uniform int is_radial;

uniform mat4 T_IM_I;
uniform mat4 T_SI_S;
uniform mat4 T_CS_C;
uniform mat4 T_NC_N;
uniform mat4 T_NR_N;

vec4 T_RC(in vec4 p_C)
{
   vec4 p_R;
   
   p_R.x = p_C.x;
   p_R.y = atan(p_C.y/p_C.z);
   p_R.z = sign(p_C.z)*sqrt(p_C.y*p_C.y + p_C.z*p_C.z);
   p_R.w = 1;
   
   return p_R;
}

// Main function, which needs to set `gl_Position`.
void main()
{
   vec4 p_CV_C, p_RV_R, t;

   // Compute distance and normals as seen from the sonar, fragment shader needs this
   n_MV_S = mat3(T_SI_S) * mat3(T_IM_I) *      n_MV_M_input; 
   p_SV_S =     (T_SI_S  *      T_IM_I  * vec4(p_MV_M,1.0)).xyz;
   
   p_CV_C = T_CS_C * T_SI_S * T_IM_I * vec4(p_MV_M, 1.0);
   
   if( is_radial == 1 ) {
      p_RV_R = T_RC(p_CV_C);
      gl_Position = T_NR_N * p_RV_R;
   } 
   // Projection
   else { 
      gl_Position = T_NC_N * p_CV_C;
      //gl_Position = T_CS_C * T_IM_I*vec4(0.5*p_MV_M, 1);
   }
   
}
"""

# Fragment shader
FS = """
#version 430

in vec3 p_SV_S; // From vertex shader
in vec3 n_MV_S; // From vertex shader

// Output variable of the fragment shader:
// RGBA components of the pixel color.
out vec4 out_color;

// Main fragment shader function.
void main()
{
   vec3  n; // Pixel normal (from vertex shader)
   float I; // Pixel intensity
   vec3  d; // Light direction (i.e. sonar origin)
   
   // Just using a Lambertian scattering model
   n = normalize(n_MV_S);
   d = normalize(-p_SV_S);
   I = dot(d,n);
   
   // ...with 30% diffused energy
   out_color = vec4(.1+.7*I, .1+.7*I, .1+.7*I, .1);
}
"""


class GLPlotWidget(QtWidgets.QOpenGLWidget):
   # default window size
   width, height = 1600, 800
   def __init__(self):
      QtWidgets.QOpenGLWidget.__init__(self)
      self.i = 0
#       
#    def initializeCL(self):
#       clinit()

   def initializeGL(self):
      """Initialize OpenGL, VBOs, upload data on the GPU, etc."""

      print('Loaded OpenGL {} with GLSL {}'.format(
         str(gl.glGetString(gl.GL_VERSION), 'utf-8'), 
         str(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION), 'utf-8')
         ))
      
      gl.glClearColor(0.2, 0.2, 0.2, 0)      # Background color
      gl.glEnable(gl.GL_DEPTH_TEST)          # Enable depth buffer
      # gl.glEnable(gl.GL_BLEND)               # Enable blending for alpha
#       gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
#       gl.glBlendEquationSeparate(gl.GL_FUNC_ADD, gl.GL_FUNC_ADD)
      
      # Compiling vertex & fragment shader
      vs = compileVertexShader(VS)
      fs = compileFragmentShader(FS)         
      self.shaders_program = createProgram(vs, fs)
      gl.glUseProgram(self.shaders_program)
      
      # Set up uniforms for sending data from host to the shaders 
      self.T_IM_I_device = gl.glGetUniformLocation(self.shaders_program,'T_IM_I')
      self.T_SI_S_device = gl.glGetUniformLocation(self.shaders_program,'T_SI_S')
      self.T_CS_C_device = gl.glGetUniformLocation(self.shaders_program,'T_CS_C')
      self.T_NC_N_device = gl.glGetUniformLocation(self.shaders_program,'T_NC_N')
      self.T_NR_N_device = gl.glGetUniformLocation(self.shaders_program,'T_NR_N')
      self.is_radial     = gl.glGetUniformLocation(self.shaders_program,'is_radial')

      # Use the custom (very simple) object loader.
      # This only handles the most basic of models
      self.model = ObjLoader()
      self.model.loadModel('manta_test_object_with_plane.obj')
      
      self.image_frame_model = ObjLoader()
      self.image_frame_model.loadModel('image_test_object.obj')

      # Put the model vertices into vertex buffer objects (VBOs)
      self.vbo_model_vertices = glvbo.VBO(self.model.vertices)
      self.vbo_image_frame_vertices = glvbo.VBO(self.image_frame_model.vertices)
      # usage=gl.DYNAMIC_DRAW,
      # target=gl.GL_ARRAY_BUFFER 
      self.vbo_model_normals = glvbo.VBO(self.model.normals)
      self.vbo_image_frame_normals = glvbo.VBO(self.image_frame_model.normals)

      # For the display, we are going to wiggle the model's roll and pitch.
      self.psi = 0  # Roll
      self.phi = 0  # Pitch
      self.theta = 0  # Yaw

      gl.glUseProgram(self.shaders_program)


   
   def drawModel(self, type, image_type, model, vbo_vertices, vbo_normals):
      
      # Bind (enable) the VBO for this model's vertices
      vbo_vertices.bind()
      
      # Connect VBO to 1st layout in shader, i.e. to vec3 p_MV_M
      gl.glEnableVertexAttribArray(0)
      
      # Describe vertex format: 3 single precision coordinates
      gl.glVertexAttribPointer(0, 3, gl.GL_FLOAT, gl.GL_FALSE, 0, None)

      # Repeat the same procedure for the normals too.
      # The normals are needed to compute lighting in the fragment shader
      vbo_normals.bind()
      gl.glEnableVertexAttribArray(1)
      gl.glVertexAttribPointer(1, 3, gl.GL_FLOAT, gl.GL_FALSE, 0, None)
      

      # Tell vertex shader whether we are rendering the sonar image or not,
      # so it can adjust its transforms accordingly
      if image_type == 'cartesian_view':
         gl.glUniform1i(self.is_radial, 0)
      elif image_type == 'radial_view':
         gl.glUniform1i(self.is_radial, 1)
      else:
         raise Exception("This should never happen")
         
      T_IM_I, T_SI_S, T_CS_C, T_NC_N, T_NR_N \
         = computeGeometry(type, self.psi, self.theta, self.phi)
      
      # Upload all the newly computed transforms to the device
      gl.glUniformMatrix4fv(self.T_IM_I_device, 1, gl.GL_TRUE, T_IM_I)
      gl.glUniformMatrix4fv(self.T_SI_S_device, 1, gl.GL_TRUE, T_SI_S)
      gl.glUniformMatrix4fv(self.T_CS_C_device, 1, gl.GL_TRUE, T_CS_C)
      gl.glUniformMatrix4fv(self.T_NC_N_device, 1, gl.GL_TRUE, T_NC_N)
      gl.glUniformMatrix4fv(self.T_NR_N_device, 1, gl.GL_TRUE, T_NR_N)

      # The grand finale: Draw all the points
      gl.glDrawArrays(gl.GL_TRIANGLES, 0, len(model.vertices))

   def paintGL(self):
      """This function gets called repeatedly. Each time we render two images:
      Left.  A normal camera view, from the sonar looking at the image center
      Right. A bird view of the image (perfect fit)
      """
      
      # Start fresh. Clear the screen buffer.
      gl.glClear(gl.GL_COLOR_BUFFER_BIT|gl.GL_DEPTH_BUFFER_BIT)
   
      for image_type,image_vals in {
         'radial_view':{
            'viewport': list(map(int,[0,0,self.width/2,self.height]))},
         'cartesian_view':{
            'viewport': list(map(int,[self.width/2,0,self.width/2,self.height]))}}.items():
         
         # Specify which part of the screen to draw in (left or right)
         gl.glViewport(*image_vals['viewport'])
               
         for type,vals in {'model':      [self.model,
                                          self.vbo_model_vertices,
                                          self.vbo_model_normals],
                          'image_frame': [self.image_frame_model,
                                          self.vbo_image_frame_vertices,
                                          self.vbo_image_frame_normals]}.items():
            
            # Let drawModel draw the current model
            self.drawModel(type, image_type, *vals)
         
      # For each iteration, we adjust the roll and pitch of the model a bit.
      # This makes the rendering more interesting to watch, but also helps
      # identify corner cases (it's easier to "think" 3D when there's movement)
      speed = 1   # Adjust at will
      self.psi     = -pi*22.5/180  * sin(self.i*pi*speed/180)  # Roll  (around x-axis)
      self.theta   = -pi*12.5/180  * sin(self.i*pi*speed/360)  # Pitch (around y-axis)
      self.phi     = -pi/180 * speed * 0                 # Yaw   (around z-axis)
      self.update()
      self.i += 1

   def resizeGL(self, width, height):
      """Called upon window resizing: Reinitialize the viewport."""
      self.width, self.height = width, height

# Set up a Qt5 window with an OpenGL widget inside it
class TestWindow(QtWidgets.QMainWindow):
   def __init__(self):
      # Initialize the QMainWindow base class
      super(TestWindow, self).__init__()
      
      # Create the GL wigdet
      self.widget = GLPlotWidget()
      
      # Put the window at screen position (100,100)
      self.setGeometry(100, 100, self.widget.width, self.widget.height)
      
      # Let GlPlotWidget occupy the entire QMainWindow 
      self.setCentralWidget(self.widget)
      
      # Pass control to the Qt5 event loop
      self.show()
     
if __name__ == '__main__':
   win = createWindow(TestWindow)
#    computeGeometry('model', 0, 0, 0)
