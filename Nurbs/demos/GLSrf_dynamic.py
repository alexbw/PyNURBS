from Nurbs import Srf, Crv
import numpy as np

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import sys

lastx = 0
lasty = 0
lastkey = ''
width = 640
height = 480
lasttime = 0.0

def maprange(val, source_range=(-100, 500), dest_range=(-5,5), clip=True):
	if clip:
		val = np.clip(val, source_range[0], source_range[1])

	# Normalize
	val = (val - source_range[0]) / (source_range[1] - source_range[0])
	# And remap
	val = val*(dest_range[1]-dest_range[0]) + dest_range[0]

	return val

def norm(val, minval=0, maxval=1, clip=True):

	val = np.clip(val)

def on_keypress(key, x, y):
	global lastkey
	lastkey = key

	if key == 'x':
		sys.exit(1)

def on_motion(x, y):
	# store the mouse coordinate
	global lastx, lasty, width, height
	lastx = float(x)
	lasty = float(y)

	# redisplay
	glutPostRedisplay()

def on_reshape(this_width, this_height):
	# setup the viewport
	global width, height
	width, height = this_width, this_height

	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	aspect = float(width)/float(height)
	if aspect < 1.:
		glOrtho(-10., 20., -10. / aspect, 20. / aspect, -25., 25.)
	else:
		glOrtho(-10. * aspect, 20. * aspect, -10., 20., -25., 25.)
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
		
def on_display():
	import time
	global lasttime

	thistime = time.time()
	print (1.0/(thistime-lasttime))
	lasttime = thistime

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()

	glCallList(1)


	global lastx, lasty, width, height, lastkey

	cam_range = (-15,15)
	camx = maprange(lastx, (0,width), cam_range)
	camy = maprange(lasty, (0,height), cam_range)
	
	glTranslate(5, 5, 0.0)
	gluLookAt(0.0, 5.0, -2.5, 0, 0, 0, 0, 1, 0)
	glRotatef(-60, 1.0, 0.0, 0.0)
	glScale(3.0, 3.0, 3.0)


	scale_a = maprange(lastx, (0, width), (1, 10), clip=False)

	cntrl = np.zeros((4,4,4), np.float)
	for u in range(4):
	    for v in range(4):
	        cntrl[0][u][v] = 2.*(u - 1.5)
	        cntrl[1][u][v] = 2.*(v - 1.5)
	        if (u == 1 or u == 2) and (v == 1 or v == 2):
	            cntrl[2][u][v] = scale_a # used to be 2.0
	        else:
	            cntrl[2][u][v] = -2.0 # used to be -2.0
	        cntrl[3][u][v] = 1.

	cntrl_color = np.zeros_like(cntrl)

	cntrl_color[0:3,:,:] = cntrl[2,:,:]
	cntrl_color[3,:,:] = 1.0
	knots = [0.,0.,0.,0.,1.,1.,1.,1.]

	srf = Srf.Srf(cntrl, knots, knots)

	nurb2 = gluNewNurbsRenderer()
	gluNurbsProperty(nurb2, GLU_SAMPLING_TOLERANCE, 50.)
	#gluNurbsProperty(nurb2, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON)
	gluNurbsProperty(nurb2, GLU_DISPLAY_MODE, GLU_FILL)
	
	gluBeginSurface(nurb2)
	# gluNurbsSurface(nurb2, srf.uknots, srf.vknots, np.transpose(cntrl_color, (1,2,0)), type=GL_MAP2_COLOR_4)
	gluNurbsSurface(nurb2, srf.uknots, srf.vknots, np.transpose(srf.cntrl, (1,2,0)), type=GL_MAP2_VERTEX_4)
	gluEndSurface(nurb2)

	# nurb1 = nurb2
	# # nurb1 = gluNewNurbsRenderer()
	# gluBeginSurface(nurb1)
	# gluEndSurface(nurb1)
	
	glutSwapBuffers()

	data = glReadPixels(0,0,width,height, GL_RGB, GL_FLOAT)
	print "Maxval: %f" % data.max()
	print data.shape

def main():

	global height, width
	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
	glutInitWindowSize(width, height)
	glutCreateWindow('GLSrf - Press MB and drag to rotate')

	glutKeyboardFunc(on_keypress)
	glutMotionFunc(on_motion)
	glutDisplayFunc(on_display)
	glutReshapeFunc(on_reshape)
	
	glMatrixMode(GL_MODELVIEW)
	
	glNewList(1, GL_COMPILE)
	
	#glMaterialfv(GL_FRONT, GL_SPECULAR, ( 0.5, 0.5, 0.5, 0.5 ))
	glMaterialfv(GL_FRONT, GL_SHININESS, 100.0)
	# glMaterialfv(GL_FRONT, GL_DIFFUSE, ( 0.7, 0.0, 0.1, 1.0 ))
	# glEnable(GL_LIGHTING)
	# glEnable(GL_LIGHT0)
	glEnable(GL_DEPTH_TEST)
	glEnable(GL_AUTO_NORMAL)
	glEnable(GL_NORMALIZE)

	glEndList()




if __name__ == '__main__':
	try:
		GLU_VERSION_1_2
	except:
		print "Need GLU 1.2 to run this demo"
		sys.exit(1)
	main()
	glutMainLoop()

