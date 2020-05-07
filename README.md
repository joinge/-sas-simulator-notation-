# Synthetic aperture sonar simulator for automatic target recognition: Basic concept and notation showcase.

![alt text](sas_simulator_notation_screenshot.gif "Screenshot")

| Author  | Jo Inge Buskenes  |
|:----|:----|
| Affiliation | University of Oslo / The Norwegian Defense Research Establishment (FFI) |
| License | None, but credits are appreciated :) |

### Documentation
| Article | Real-time Synthetic Aperture Sonar Simulation for Automatic Target Recognition: Notation and Application   |
|:----|:----|
| Authors | *Buskenes J.I, Midelfart, H.* |
| Journal | IEEE Journal of Oceanic Engineering |
| Date    | May 2020 |

### Files
|    |    |
|----|----|
| sas_simulation_notation.py                | Simulator code |
| sas_simulation_notation_screenshot.png    | Screenshot of a running simulation |
| ObjLoader.py                              | Helper class to load Wavefront files (*.obj) |
| image_test_object.obj                     | 3D model of the image frame |
| manta_test_object_with_plane.obj          | 3D model of a manta mine |
     
### Credits
&ensp;&ensp;&ensp;Adopted from a template for pyopengl found at:
   https://github.com/lukecampbell/pyopengl-examples

### Prerequisites

Tested on various Ubuntu installations, with these package requirements:\
&ensp;&ensp;&ensp;For Python3: python3-numpy python3-pyqt5 python3-pyqt5.qtopengl python3-opengl\
&ensp;&ensp;&ensp;For Python2: python-numpy python-pyqt5 python-pyqt5.qtopengl python-opengl
   
GL version must be at least 4.30, force MESA to use it by setting these environment variables:\
&ensp;&ensp;&ensp;export MESA_GLSL_VERSION_OVERRIDE=430\
&ensp;&ensp;&ensp;export MESA_GL_VERSION_OVERRIDE=4.3
