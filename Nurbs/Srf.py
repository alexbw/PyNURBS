import math
from _Bas import  bspkntins, bspdegelev, bspbezdecom, bspeval # Lowlevel functions
from Util import translate, scale, roty, rotx, rotz
import Crv

dependencies = '''This module requires:
	Numeric Python
'''

try:
    import numpy as np
except ImportError, value:
	print dependencies
	raise

NURBSError = 'NURBSError'

class Srf:
    '''Construct a NURB surface structure, and check the format.
    
 The NURB surface is represented by a 4 dimensional b-spline.

 INPUT:

    cntrl  - Control points, homogeneous coordinates (wx,wy,wz,w)
        For a surface [dim,nu,nv] matrix
            where nu  is along the u direction and
            nv  is along the v direction.
            dim is the dimension valid options are:
            2 .... (x,y)        2D cartesian coordinates
            3 .... (x,y,z)      3D cartesian coordinates   
            4 .... (wx,wy,wz,w) 4D homogeneous coordinates

    uknots - Knot sequence along the parametric u direction.
    vknots - Knot sequence along the paramteric v direction.

 NOTES:

    Its assumed that the input knot sequences span the
    interval [0.0,1.0] and are clamped to the control
    points at the end by a knot multiplicity equal to
    the spline order.'''

    def __init__(self, cntrl, uknots, vknots):
        self._bezier = None
        cntrl = np.asarray(cntrl, np.float)
        (dim, nu, nv) = cntrl.shape
        if dim < 2 or dim > 4:
            raise NURBSError, 'Illegal control point format'
        elif dim < 4:
            self.cntrl = np.zeros((4, nu, nv), np.float)
            self.cntrl[0:dim,:,:] = cntrl
            self.cntrl[-1,:,:] = np.ones((nu,nv), np.float)
        else:
            self.cntrl = cntrl
            
        # Force the u knot sequence to be a vector in ascending order
        # and normalise between [0.0,1.0]
        uknots = np.sort(np.asarray(uknots, np.float))
        nku = uknots.shape[0]
        uknots = (uknots - uknots[0])/(uknots[-1] - uknots[0])
        if uknots[0] == uknots[-1]:
            raise NURBSError, 'Illegal uknots sequence'
        self.uknots = uknots
        
        # Force the v knot sequence to be a vector in ascending order
        # and normalise between [0.0,1.0]  
        vknots = -np.sort(-np.asarray(vknots, np.float))
        nkv = vknots.shape[0]
        vknots = (vknots - vknots[0])/(vknots[-1] - vknots[0])
        if vknots[0] == vknots[-1]:
            raise NURBSError, 'Illegal vknots sequence'
        self.vknots = vknots
        
        # Spline Degree
        self.degree = [nku-nu-1, nkv-nv-1]
        if self.degree[0] < 0 or self.degree[1] < 0:
            raise NURBSError, 'NURBS order must be a positive integer'

    def trans(self, mat):
        "Apply the 4D transform matrix to the NURB control points."
        for v in range(self.cntrl.shape[2]):
            self.cntrl[:,:,v] = np.dot(mat, self.cntrl[:,:,v])
        
    def swapuv(self):
        "Swap u and v parameters."
        self.cntrl = np.transpose(self.cntrl,(0,2,1))
        temp = self.uknots[:]
        self.uknots = self.vknots[:]
        self.vknots = temp
        udegree, vdegree = self.degree
        self.degree[0] = vdegree
        self.degree[1] = udegree
        
    def reverse(self):
        "Reverse evaluation directions."
        coefs = self.cntrl[:,:,::-1]
        self.cntrl = coefs[:,::-1,:]
        self.uknots = 1. - self.uknots[::-1]
        self.vknots = 1. - self.vknots[::-1]

    def extractV(self, u):
        "Extract curve in v-direction at parameter u."
        if np.less(u, 0.) or np.greater(u, 1.):
                raise NURBSError, 'Out of parameter range [0,1]'
        if u == 0.:
            cntrl = self.cntrl[:,0,:]
            knots = self.vknots[:]
        elif u == 1.:
            cntrl = self.cntrl[:,-1,:]
            knots = self.vknots[:]
        else:
            uknots = np.repeat(np.asarray([u], np.float),[self.degree[1]*(self.cntrl.shape[2] + 1)])
            coefs = np.transpose(self.cntrl,(0, 2, 1))
            coefs = np.resize(coefs,(4*self.cntrl.shape[2], self.cntrl.shape[1]))
            coefs, knots = bspkntins(self.degree[0], coefs, self.uknots, uknots)
            coefs = np.resize(coefs, (4, self.cntrl.shape[2], coefs.shape[1]))
            cntrl = np.transpose(coefs,(0,2,1))
            i = 0
            j = knots[0]
            for k in knots[1:]:
                if k == u:
                    break
                elif k != j:
                    i += 1
                    j = k
        return Crv.Crv(cntrl[:,i,:], self.vknots[:])

    def extractU(self, v):
        "Extract curve in u-direction at parameter v."
        if np.less(v, 0.) or np.greater(v, 1.):
                raise NURBSError, 'Out of parameter range [0,1]'
        if v == 0.:
            cntrl = self.cntrl[:,:,0]
            knots = self.uknots[:]
        elif v == 1.:
            cntrl = self.cntrl[:,:,-1]
            knots = self.uknots[:]
        else:
            vknots = np.repeat(np.asarray([v], np.float),[self.degree[0]*(self.cntrl.shape[1] + 1)])
            coefs = np.resize(self.cntrl,(4*self.cntrl.shape[1], self.cntrl.shape[2]))
            coefs, knots = bspkntins(self.degree[1], coefs, self.vknots, vknots)
            cntrl = np.resize(coefs, (4, self.cntrl.shape[1], coefs.shape[1]))
            i = 0
            j = knots[0]
            for k in knots[1:]:
                if k == v:
                    break
                elif k != j:
                    i += 1
                    j = k
        return Crv.Crv(cntrl[:,:,i], self.uknots[:])
    
    def kntins(self, uknots, vknots = None):
        """Insert new knots into the surface
	uknots - knots to be inserted along u direction
	vknots - knots to be inserted along v direction
	NOTE: No knot multiplicity will be increased beyond the order of the spline"""
        if len(vknots):
            # Force the v knot sequence to be a vector in ascending order
            vknots = np.sort(np.asarray(vknots, np.float))
            if np.less(vknots, 0.) or np.greater(vknots, 1.):
                raise NURBSError, 'Illegal vknots sequence'
            coefs = np.resize(self.cntrl,(4*self.cntrl.shape[1], self.cntrl.shape[2]))
            coefs, self.vknots = bspkntins(self.degree[1], coefs, self.vknots, vknots)
            self.cntrl = np.resize(coefs, (4, self.cntrl.shape[1], coefs.shape[1]))
        if len(uknots):
            # Force the u knot sequence to be a vector in ascending order
            uknots = np.sort(np.asarray(uknots, np.float))
            if np.less(uknots, 0.) or np.greater(uknots, 1.):
                raise NURBSError, 'Illegal uknots sequence'
            coefs = np.transpose(self.cntrl,(0, 2, 1))
            coefs = np.resize(coefs,(4*self.cntrl.shape[2], self.cntrl.shape[1]))
            coefs, self.uknots = bspkntins(self.degree[0], coefs, self.uknots, uknots)
            coefs = np.resize(coefs, (4, self.cntrl.shape[2], coefs.shape[1]))
            self.cntrl = np.transpose(coefs,(0,2,1))
            
    def degelev(self, utimes, vtimes = None):
        """Degree elevate the surface.
	utimes - degree elevate utimes along u direction.
	vtimes - degree elevate vtimes along v direction."""
        if vtimes:
            if vtimes < 0:
                raise NURBSError, 'Degree must be positive'
            coefs = np.resize(self.cntrl,(4*self.cntrl.shape[1], self.cntrl.shape[2]))
            coefs, vknots, nh = bspdegelev(self.degree[1], coefs, self.vknots, vtimes)
            coefs = coefs[:,:nh + 1]
            self.vknots = vknots[:nh + self.degree[1] + vtimes + 2]
            self.degree[1] += vtimes
            self.cntrl = np.resize(coefs, (4, self.cntrl.shape[1], coefs.shape[1]))
        if utimes:
            if utimes < 0:
                raise NURBSError, 'Degree must be positive'
            coefs = np.transpose(self.cntrl,(0, 2, 1))
            coefs = np.resize(coefs,(4*self.cntrl.shape[2], self.cntrl.shape[1]))
            coefs, uknots, nh = bspdegelev(self.degree[0], coefs, self.uknots, utimes)
            coefs = coefs[:,:nh + 1]
            self.uknots = uknots[:nh + self.degree[0] + utimes + 2]
            self.degree[0] += utimes
            coefs = np.resize(coefs, (4, self.cntrl.shape[2], coefs.shape[1]))
            self.cntrl = np.transpose(coefs,(0,2,1))

    def bezier(self, update = None):
        "Decompose surface to bezier patches and return overlaping control points."
        if update or not self._bezier:
            cntrl = np.resize(self.cntrl,(4*self.cntrl.shape[1], self.cntrl.shape[2]))
            cntrl = bspbezdecom(self.degree[1], cntrl, self.vknots)
            cntrl = np.resize(cntrl, (4, self.cntrl.shape[1], cntrl.shape[1]))
            temp1 = cntrl.shape[1]
            temp2 = cntrl.shape[2]
            cntrl = np.transpose(cntrl,(0, 2, 1))
            cntrl = np.resize(cntrl,(4*temp2, temp1))
            cntrl = bspbezdecom(self.degree[0], cntrl, self.uknots)
            cntrl = np.resize(cntrl, (4, temp2, cntrl.shape[1]))
            self._bezier = np.transpose(cntrl,(0, 2, 1))
        return self._bezier
        
    def bounds(self):
        "Return the bounding box for the surface."
        w = self.cntrl[3,:,:]
        cx = np.sort(np.ravel(self.cntrl[0,:,:]/w))
        cy = np.sort(np.ravel(self.cntrl[1,:,:]/w))
        cz = np.sort(np.ravel(self.cntrl[2,:,:]/w))
        return np.asarray([cx[0], cy[0], cz[0],
                                cx[-1], cy[-1], cz[-1]], np.float)
        
    def pnt3D(self, ut, vt = None):
        """Evaluate parametric point[s] and return 3D cartesian coordinate[s]
	If only ut is given then we will evaluate at scattered points.
	ut(0,:) represents the u direction.
	ut(1,:) represents the v direction.
	If both parameters are given then we will evaluate over a [u,v] grid."""
        val = self.pnt4D(ut, vt)
        if len(val.shape) < 3:
            return val[0:3,:]/np.resize(val[-1,:], (3, val.shape[1]))
        else: #FIX!
            return val[0:3,:,:]/np.resize(val[-1,:,:], (3, val.shape[1], val.shape[2]))
        
    def pnt4D(self, ut, vt = None):
        """Evaluate parametric point[s] and return 4D homogeneous coordinates.
	If only ut is given then we will evaluate at scattered points.
	ut(0,:) represents the u direction.
	ut(1,:) represents the v direction.
	If both parameters are given then we will evaluate over a [u,v] grid."""
        ut = np.asarray(ut, np.float)
        if np.less(ut, 0.) or np.greater(ut, 1.):
            raise NURBSError, 'NURBS curve parameter out of range [0,1]'
        
        if vt: #FIX!
            # Evaluate over a [u,v] grid
            vt = np.asarray(vt, np.float)
            if np.less(vt, 0.) or np.greater(vt, 1.):
                raise NURBSError, 'NURBS curve parameter out of range [0,1]'
    
            val = np.resize(self.cntrl,(4*self.cntrl.shape[1],self.cntrl.shape[2]))
            val = bspeval(self.degree[1], val, self.vknots, vt)
            val = np.resize(val,(4, self.cntrl.shape[1], vt.shape[0]))
    
            val = np.transpose(val,(0,2,1))
    
            val = np.resize(self.cntrl,(4*vt.shape[0],self.cntrl.shape[1]))
            val = bspeval(self.degree[0], val, self.uknots, ut)
            val = np.resize(val,(4, vt.shape[0], ut.shape[0]))
       
            return np.transpose(val,(0,2,1)) 

        # Evaluate at scattered points
        nt = ut.shape[1]
        uval = np.resize(self.cntrl,(4*self.cntrl.shape[1],self.cntrl.shape[2]))
        uval = bspeval(self.degree[1],uval,self.vknots,ut[1,:])
        uval = np.resize(uval,(4, self.cntrl.shape[1], nt))

        val = np.zeros((4,nt), np.float)
        for v in range(nt):
            val[:,v] = bspeval(self.degree[0],np.resize(uval[:,:,v],(4,self.cntrl.shape[1])),
                                self.uknots, (ut[0,v],))[:,0]
        return val
            
    def plot(self, n = 50, iso = 8):
        """A simple plotting function based on dislin for debugging purpose.
	n = number of subdivisions. iso = number of iso line to plot in each dir.
	TODO: plot ctrl poins and knots."""
        try:
            import dislin
        except ImportError, value:
            print 'dislin plotting library not available'
            return

        maxminx = np.sort(np.ravel(self.cntrl[0,:,:])/np.ravel(self.cntrl[3,:,:]))
        minx = maxminx[0]
        maxx = maxminx[-1]
        if minx == maxx:
            minx -= 1.
            maxx += 1.

        maxminy = np.sort(np.ravel(self.cntrl[1,:,:])/np.ravel(self.cntrl[3,:,:]))
        miny = maxminy[0]
        maxy = maxminy[-1]
        if miny == maxy:
            miny -= 1.
            maxy += 1.

        maxminz = np.sort(np.ravel(self.cntrl[2,:,:])/np.ravel(self.cntrl[3,:,:]))
        minz = maxminz[0]
        maxz = maxminz[-1]
        if minz == maxz:
            minz -= 1.
            maxz += 1.
                
        dislin.metafl('cons')
        dislin.disini()
        dislin.hwfont()
        dislin.pagera()
        dislin.name('X-axis', 'X')
        dislin.name('Y-axis', 'Y')
        dislin.name('Z-axis', 'Z')
        dislin.graf3d(minx, maxx, 0 , abs((maxx-minx)/4.),
                      miny, maxy, 0 , abs((maxy-miny)/4.),
                      minz, maxz, 0 , abs((maxz-minz)/4.))
            
        dislin.color('yellow')
        pnts0 = self.pnt3D([np.arange(n + 1, typecode = np.float)/n,
                            np.zeros(n + 1,np.float)])
        pnts1 = self.pnt3D([np.arange(n + 1, typecode = np.float)/n,
                            np.ones(n + 1,np.float)])
        pnts2 = self.pnt3D([np.zeros(n + 1,np.float),
                            np.arange(n + 1, typecode = np.float)/n])
        pnts3 = self.pnt3D([np.ones(n + 1,np.float),
                            np.arange(n + 1, typecode = np.float)/n])
        dislin.curv3d(pnts0[0,:], pnts0[1,:], pnts0[2,:], n+1)
        dislin.curv3d(pnts1[0,:], pnts1[1,:], pnts1[2,:], n+1)
        dislin.curv3d(pnts2[0,:], pnts2[1,:], pnts2[2,:], n+1)
        dislin.curv3d(pnts3[0,:], pnts3[1,:], pnts3[2,:], n+1)
            
        dislin.color('red')
        step = 1./iso
        for uv in np.arange(step, 1., step):
            pnts = self.pnt3D([np.arange(n + 1, typecode = np.float)/n,
                               np.zeros(n + 1,np.float) + uv])
            dislin.curv3d(pnts[0,:], pnts[1,:], pnts[2,:], n+1)
            pnts = self.pnt3D([np.zeros(n + 1,np.float) + uv,
                               np.arange(n + 1, typecode = np.float)/n])
            dislin.curv3d(pnts[0,:], pnts[1,:], pnts[2,:], n+1)
            
        dislin.disfin()

    def __call__(self, *args):
        return self.pnt3D(args[0])

    def __repr__(self):
        return 'Nurbs surface:\ndegree: %s\ncntrl: %s\nuknots: %s\nvknots: %s' % (`self.degree`,`self.cntrl`,`self.uknots`,`self.vknots`)

class Bilinear(Srf):
    '''Constructs a NURB surface defined by 4 points to form a bilinear surface.
    Inputs:
        p11 - 3D coordinate of the corner point
        p12 - 3D coordinate of the corner point 
        p21 - 3D coordinate of the corner point 
        p22 - 3D coordinate of the corner point 

    The position of the corner points

        ^ v direction
        |
        ----------------
        |p21        p22|
        |              |
        |    SRF       |
        |              |
        |p11        p12|
        -------------------> u direction'''       
    def __init__(self, p00 = [-.5,-.5], p01 = [.5, -.5], p10 = [-.5, .5], p11 = [.5, .5]):
        coefs = np.zeros((4,2,2), np.float)
        coefs[3,:,:] = np.ones((2,2), np.float)
        coefs[0:len(p00),0,0] = p00
        coefs[0:len(p01),0,1] = p01
        coefs[0:len(p10),1,0] = p10
        coefs[0:len(p11),1,1] = p11
        Srf.__init__(self, coefs, [0.,0.,1.,1.], [0., 0., 1., 1.])

class Extrude(Srf):
    '''Construct a NURBS surface by extruding a NURBS curve along the vector (dx,dy,dz).
    INPUT:
        curve   - NURBS curve to be extruded.
        vector  - Extrusion vector.'''
    def __init__(self, crv, vector):
        if not isinstance(crv, Crv.Crv):
            raise NURBSError, 'Parameter crv not derived from Crv class!'
        coefs = np.zeros((4,crv.cntrl.shape[1],2), np.float)
        coefs[:,:,0] = crv.cntrl
        coefs[:,:,1] = np.dot(translate(vector), crv.cntrl)
        Srf.__init__(self, coefs, crv.uknots, [0., 0., 1., 1.])

class Revolve(Srf):
    '''Construct a surface by revolving the profile curve
    around an axis defined by a point and a unit vector.

    srf = nrbrevolve(curve,point,vector[,theta])

    INPUT:

    crv    - NURB profile curve.
    pnt    - coordinate of the point. (default: [0, 0, 0])
    vector - rotation axis. (default: [1, 0, 0])
    theta  - angle to revolve curve (default: 2*pi).'''
    def __init__(self, crv, pnt = [0., 0., 0.], vector = [1., 0., 0.], theta = 2.*math.pi):
        if not isinstance(crv, Crv.Crv):
            raise NURBSError, 'Parameter crv not derived from Crv class!'
        # Translate and rotate the curve into alignment with the z-axis
        T = translate(-np.asarray(pnt, np.float))
        # Normalize vector
        vector = np.asarray(vector, np.float)
        len = np.sqrt(np.add.reduce(vector*vector))
        if len == 0:
            raise ZeroDivisionError, "Can't normalize a zero-length vector"
        vector = vector/len
        if vector[0] == 0.:
            angx = 0.
        else:
            angx = math.atan2(vector[0], vector[2])
        RY = roty(-angx)
        vectmp = np.ones((4,), np.float)
        vectmp[0:3] = vector
        vectmp = np.dot(RY, vectmp)
        if vectmp[1] == 0.:
            angy = 0.
        else:
            angy = math.atan2(vector[1], vector[2])
        RX = rotx(angy)
        crv.trans(np.dot(RX, np.dot(RY, T)))
        arc = Crv.Arc(1., [0., 0., 0.], 0., theta)

        narc = arc.cntrl.shape[1]
        ncrv = crv.cntrl.shape[1]
        coefs = np.zeros((4, narc, ncrv), np.float)
        angle = np.arctan2(crv.cntrl[1,:], crv.cntrl[0,:])
        vectmp = crv.cntrl[0:2,:]
        radius = np.sqrt(np.add.reduce(vectmp*vectmp))

        for i in xrange(0, ncrv):
            coefs[:,:,i] = np.dot(rotz(angle[i]),
                                       np.dot(translate((0., 0., crv.cntrl[2,i])),
                                                   np.dot(scale((radius[i], radius[i])), arc.cntrl)))
            coefs[3,:,i] = coefs[3,:,i] * crv.cntrl[3,i]
        Srf.__init__(self, coefs, arc.uknots, crv.uknots)
        T = translate(pnt)
        RX = rotx(-angy)
        RY = roty(angx)
        self.trans(np.dot(T, np.dot(RY, RX)))

class Ruled(Srf):
    '''Constructs a ruled surface between two NURBS curves.
    INPUT:
        crv1 - first  NURBS curve
        crv2 - second NURBS curve
    NOTE:
    The ruled surface is ruled along the V direction.'''
    def __init__(self, crv1, crv2):
        if not isinstance(crv1, Crv.Crv):
                raise NURBSError, 'Parameter crv1 not derived from Crv class!'
        if not isinstance(crv2, Crv.Crv):
                raise NURBSError, 'Parameter crv2 not derived from Crv class!'
        # ensure both curves have a common degree
        d = max(crv1.degree, crv2.degree)
        crv1.degelev(d - crv1.degree)
        crv2.degelev(d - crv2.degree)
        # merge the knot vectors, to obtain a common knot vector
        k1 = crv1.uknots
        k2 = crv2.uknots
        ku = []
        for item in k1:
            if not np.sometrue(np.equal(k2, item)):
                if item not in ku:
                    ku.append(item)
        for item in k2:
            if not np.sometrue(np.equal(k1, item)):
                if item not in ku:
                    ku.append(item)
        ku = np.sort(np.asarray(ku, np.float))
        n = ku.shape[0]
        ka = np.array([], np.float)
        kb = np.array([], np.float)
        for i in range(0, n):
            i1 = np.compress(np.equal(k1, ku[i]), k1).shape[0]
            i2 = np.compress(np.equal(k2, ku[i]), k2).shape[0]
            m = max(i1, i2)
            ka = np.concatenate((ka , ku[i] * np.ones((m - i1,), np.float)))
            kb = np.concatenate((kb , ku[i] * np.ones((m - i2,), np.float)))
        crv1.kntins(ka)
        crv2.kntins(kb)
        coefs = np.zeros((4, crv1.cntrl.shape[1], 2), np.float)
        coefs[:,:,0] = crv1.cntrl
        coefs[:,:,1] = crv2.cntrl
        Srf.__init__(self, coefs, crv1.uknots, [0., 0., 1., 1.])

class Coons(Srf):
    '''Construction of a bilinearly blended Coons patch.
    INPUTS:
        u1, u2, v1 ,v2 - four NURBS curves defining boundaries.
    NOTES:
        The orientation of the four NURBS boundary curves.

         ^ v direction
         |
         |     u2
         ------->--------
         |              |
         |              |
      v1 ^   Surface    ^ v2
         |              |
         |              |
         ------->-----------> u direction
               u1'''
    def __init__(self, u1, u2, v1, v2):
        if not isinstance(u1, Crv.Crv):
                raise NURBSError, 'Parameter u1 not derived from Crv class!'
        if not isinstance(u2, Crv.Crv):
                raise NURBSError, 'Parameter u2 not derived from Crv class!'
        if not isinstance(v1, Crv.Crv):
                raise NURBSError, 'Parameter v1 not derived from Crv class!'
        if not isinstance(v2, Crv.Crv):
                raise NURBSError, 'Parameter v2 not derived from Crv class!'
        r1 = Ruled(u1, u2)
        r2 = Ruled(v1, v2)
        r2.swapuv()
        t = Bilinear(u1.cntrl[:,0], u1.cntrl[:,-1], u2.cntrl[:,0], u2.cntrl[:,-1])
        # Raise all surfaces to a common degree
        du = max(r1.degree[0], r2.degree[0], t.degree[0])
        dv = max(r1.degree[1], r2.degree[1], t.degree[1])
        r1.degelev(du - r1.degree[0], dv - r1.degree[1])
        r2.degelev(du - r2.degree[0], dv - r2.degree[1])
        t.degelev(du - t.degree[0], dv - t.degree[1])
        # Merge the knot vectors, to obtain a common knot vector
        # uknots:
        k1 = r1.uknots
        k2 = r2.uknots
        k3 = t.uknots
        k = []
        for item in k1:
            if not np.sometrue(np.equal(k2, item)):
                if not np.sometrue(np.equal(k3, item)):
                    if item not in k:
                        k.append(item)
        for item in k2:
            if not np.sometrue(np.equal(k1, item)):
                if not np.sometrue(np.equal(k3, item)):
                    if item not in k:
                        k.append(item)
        for item in k3:
            if not np.sometrue(np.equal(k1, item)):
                if not np.sometrue(np.equal(k2, item)):
                    if item not in k:
                        k.append(item)        
        k = np.sort(np.asarray(k, np.float))
        n = k.shape[0]
        kua = np.array([], np.float)
        kub = np.array([], np.float)
        kuc = np.array([], np.float)
        for i in range(0, n):
            i1 = np.compress(np.equal(k1, k[i]), k1).shape[0]
            i2 = np.compress(np.equal(k2, k[i]), k2).shape[0]
            i3 = np.compress(np.equal(k3, k[i]), k3).shape[0]
            m = max(i1, i2, i3)
            kua = np.concatenate((kua , k[i] * np.ones((m - i1,), np.float)))
            kub = np.concatenate((kub , k[i] * np.ones((m - i2,), np.float)))
            kuc = np.concatenate((kuc , k[i] * np.ones((m - i3,), np.float)))

        # vknots:
        k1 = r1.vknots
        k2 = r2.vknots
        k3 = t.vknots
        k = []
        for item in k1:
            if not np.sometrue(np.equal(k2, item)):
                if not np.sometrue(np.equal(k3, item)):
                    if item not in k:
                        k.append(item)
        for item in k2:
            if not np.sometrue(np.equal(k1, item)):
                if not np.sometrue(np.equal(k3, item)):
                    if item not in k:
                        k.append(item)
        for item in k3:
            if not np.sometrue(np.equal(k1, item)):
                if not np.sometrue(np.equal(k2, item)):
                    if item not in k:
                        k.append(item)        
        k = np.sort(np.asarray(k, np.float))
        n = k.shape[0]
        kva = np.array([], np.float)
        kvb = np.array([], np.float)
        kvc = np.array([], np.float)
        for i in range(0, n):
            i1 = np.compress(np.equal(k1, k[i]), k1).shape[0]
            i2 = np.compress(np.equal(k2, k[i]), k2).shape[0]
            i3 = np.compress(np.equal(k3, k[i]), k3).shape[0]
            m = max(i1, i2, i3)
            kva = np.concatenate((kva , k[i] * np.ones((m - i1,), np.float)))
            kvb = np.concatenate((kvb , k[i] * np.ones((m - i2,), np.float)))
            kvc = np.concatenate((kvc , k[i] * np.ones((m - i3,), np.float)))

        r1.kntins(kua, kva)
        r2.kntins(kub, kvb)
        t.kntins(kuc, kvc)
        coefs = np.zeros((4 , t.cntrl.shape[1], t.cntrl.shape[2]), np.float)
        coefs[0,:,:] = r1.cntrl[0,:,:] + r2.cntrl[0,:,:] - t.cntrl[0,:,:]
        coefs[1,:,:] = r1.cntrl[1,:,:] + r2.cntrl[1,:,:] - t.cntrl[1,:,:]
        coefs[2,:,:] = r1.cntrl[2,:,:] + r2.cntrl[2,:,:] - t.cntrl[2,:,:]
        coefs[3,:,:] = r1.cntrl[3,:,:] + r2.cntrl[3,:,:] - t.cntrl[3,:,:]
        Srf.__init__(self, coefs, r1.uknots, r1.vknots)
        
if __name__ == '__main__':
    '''cntrl = np.zeros((4,4,4), np.float)
    for u in range(4):
        for v in range(4):
            cntrl[0][u][v] = 2.*(u - 1.5)
            cntrl[1][u][v] = 2.*(v - 1.5)
            if (u == 1 or u == 2) and (v == 1 or v == 2):
                cntrl[2][u][v] = 2.
            else:
                cntrl[2][u][v] = -2.
            cntrl[3][u][v] = 1.
    knots = [0.,0.,0.,0.,1.,1.,1.,1.]
    s = Srf(cntrl,knots,knots)'''
    #s = Bilinear([0,0,0], [2,0,0], [0,1,0], [2,1,1.8])
    c = Crv.Crv([[0,30,6,90],[0,0,50,50]],[0,0,0,0,1,1,1,1])
    s = Revolve(c)
    #s = Extrude(c,[0,0,5])
    '''cntrl = [[-50., -75., 25., 0., -25., 75., 50.],
             [25., 50., 50., 0., -50., -50., 25.]]
    knots = [0., 0., 0., .2, .4, .6, .8, 1., 1., 1.]
    c = Crv.Crv(cntrl, knots)
    s = Extrude(c,[0,0,5])'''
    
    '''pnts = [[0., 3., 4.5, 6.5, 8., 10.],
            [0., 0., 0., 0., 0., 0.],
            [2., 2., 7., 4., 7., 9.]]   
    crv1 = Crv.Crv(pnts, [0., 0., 0., 1./3., 0.5, 2./3., 1., 1., 1.])
    
    pnts= [[0., 3., 5., 8., 10.],
           [10., 10., 10., 10., 10.],
           [3., 5., 8., 6., 10.]]
    crv2 = Crv.Crv(pnts, [0., 0., 0., 1./3., 2./3., 1., 1., 1.])
    
    pnts= [[0.,0., 0., 0.],
           [0., 3., 8., 10.],
           [2., 0., 5., 3.]]
    crv3 = Crv.Crv(pnts, [0., 0., 0., 0.5, 1., 1., 1.])
    
    pnts= [[10., 10., 10., 10., 10.],
           [0., 3., 5., 8., 10.],
           [9., 7., 7., 10., 10.]]
    crv4 = Crv.Crv(pnts, [0., 0., 0., 0.25, 0.75, 1., 1., 1.])
    
    s = Coons(crv1, crv2, crv3, crv4)'''
    s.plot()
